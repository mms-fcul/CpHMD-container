#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

IMAGE=""
REPO_ROOT="${DEFAULT_REPO_ROOT}"
GPU_DEVICE=0
RUN_APP_SMOKE=1
CHECK_RUNTIME_HASHES=1
RUN_IMAGE_TEST=1

usage() {
    cat <<USAGE
Usage: $(basename "$0") --image IMAGE [options]

Validate the x86-64 CpHMD image containing isolated CPU and CUDA GROMACS
2026.3 installations and the public helper applications.

Options:
  --image FILE          SIF image to validate. Required.
  --repo-root DIR       Source repository used for embedded-runtime comparison.
                        Default: ${DEFAULT_REPO_ROOT}
  --gpu-device          Add --nv, require a visible NVIDIA device and inspect it.
  --skip-app-smoke      Verify app registration/backing files but do not copy
                        sample files into a temporary directory.
  --skip-runtime-hash   Do not compare embedded runtime hashes with the source.
  --skip-image-test     Do not run the image's built-in Apptainer test.
  -h, --help            Show this help.

This validator does not run a scientific CpHMD integration simulation. CPU,
GPU, restart and reduced-titration integration tests remain separate gates.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --image) IMAGE="${2:?Missing image path}"; shift 2 ;;
        --repo-root) REPO_ROOT="${2:?Missing repository root}"; shift 2 ;;
        --gpu-device) GPU_DEVICE=1; shift ;;
        --skip-app-smoke) RUN_APP_SMOKE=0; shift ;;
        --skip-runtime-hash) CHECK_RUNTIME_HASHES=0; shift ;;
        --skip-image-test) RUN_IMAGE_TEST=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

for command_name in apptainer sha256sum awk grep sort comm; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: ${command_name}" >&2
        exit 1
    }
done

[[ -n "${IMAGE}" ]] || {
    echo "ERROR: --image is required." >&2
    exit 2
}
IMAGE="$(readlink -f "${IMAGE}")"
[[ -s "${IMAGE}" ]] || {
    echo "ERROR: image missing or empty: ${IMAGE}" >&2
    exit 1
}

REPO_ROOT="$(cd "${REPO_ROOT}" && pwd)"

runtime=(apptainer exec --cleanenv)
run_app=(apptainer run --cleanenv)
if (( GPU_DEVICE )); then
    runtime+=(--nv)
    run_app+=(--nv)
fi
runtime+=("${IMAGE}")

cpu=/gromacs/bin/gmx
gpu=/gromacs/bin/gmx-gpu

expected_apps=(
    add-mol-db
    extract-databases
    extract-force-fields
    extract-protonation-script
    extract-tautomers
    get-settings
    get-tutorials
    make-sites
    pdb2cphmd
    protonation
)

if (( RUN_IMAGE_TEST )); then
    echo "=== Apptainer image test ==="
    apptainer test "${IMAGE}"
    echo
fi

echo "=== Image labels ==="
apptainer inspect "${IMAGE}"

echo
echo "=== Registered applications ==="
mapfile -t observed_apps < <(
    apptainer inspect --list-apps "${IMAGE}" \
        | awk 'NF {print $1}' \
        | sort -u
)
printf '%s\n' "${observed_apps[@]}"

expected_file="$(mktemp)"
observed_file="$(mktemp)"
cleanup_files=("${expected_file}" "${observed_file}")
cleanup() {
    local path
    for path in "${cleanup_files[@]}"; do
        [[ -e "${path}" ]] && rm -rf -- "${path}"
    done
}
trap cleanup EXIT

printf '%s\n' "${expected_apps[@]}" | sort -u > "${expected_file}"
printf '%s\n' "${observed_apps[@]}" | sort -u > "${observed_file}"

if ! diff_output="$(comm -3 "${expected_file}" "${observed_file}")"; then
    :
fi
if [[ -n "${diff_output}" ]]; then
    echo "ERROR: registered app set differs from the expected interface." >&2
    echo "Lines without indentation are missing; indented lines are unexpected:" >&2
    printf '%s\n' "${diff_output}" >&2
    exit 1
fi

echo
echo "=== Image-wide environment ==="
"${runtime[@]}" bash -lc '
    printf "PATH=%s\n" "${PATH}"
    printf "LD_LIBRARY_PATH=%s\n" "${LD_LIBRARY_PATH:-unset}"
    printf "PLUMED_KERNEL=%s\n" "${PLUMED_KERNEL:-unset}"
    printf "GMX_CPU_BIN=%s\n" "${GMX_CPU_BIN:-unset}"
    printf "GMX_CUDA_BIN=%s\n" "${GMX_CUDA_BIN:-unset}"
    printf "GMXLIB=%s\n" "${GMXLIB:-unset}"
'

echo
echo "=== CPU launcher ==="
cpu_output="$("${runtime[@]}" "${cpu}" --version 2>&1)"
printf '%s\n' "${cpu_output}"
grep -q 'GROMACS version:.*2026.3' <<< "${cpu_output}"
grep -q 'GPU support:.*disabled' <<< "${cpu_output}"
grep -q 'SIMD instructions:.*AVX2_256' <<< "${cpu_output}"
grep -q 'Plumed support:.*enabled' <<< "${cpu_output}"

echo
echo "=== CUDA launcher ==="
gpu_output="$("${runtime[@]}" "${gpu}" --version 2>&1)"
printf '%s\n' "${gpu_output}"
grep -q 'GROMACS version:.*2026.3' <<< "${gpu_output}"
grep -q 'GPU support:.*CUDA' <<< "${gpu_output}"
grep -q 'SIMD instructions:.*AVX2_256' <<< "${gpu_output}"
grep -q 'Plumed support:.*enabled' <<< "${gpu_output}"

echo
echo "=== Runtime library isolation ==="
"${runtime[@]}" bash -lc '
    set -eu
    cpu_root=/gromacs/gromacs-2026.3
    gpu_root=/gromacs/gromacs-2026.3-GPU

    cpu_lib="$(
        LD_LIBRARY_PATH="${cpu_root}/lib:${cpu_root}/lib64:/opt/plumed/lib:/opt/plumed/lib64" \
            ldd "${cpu_root}/bin/gmx" \
            | while read -r library arrow path remainder; do
                case "${library}" in
                    libgromacs*) printf "%s\n" "${path}"; break ;;
                esac
              done
    )"
    gpu_lib="$(
        LD_LIBRARY_PATH="${gpu_root}/lib:${gpu_root}/lib64:/opt/plumed/lib:/opt/plumed/lib64:/usr/local/cuda/lib64" \
            ldd "${gpu_root}/bin/gmx" \
            | while read -r library arrow path remainder; do
                case "${library}" in
                    libgromacs*) printf "%s\n" "${path}"; break ;;
                esac
              done
    )"

    cpu_lib="$(readlink -f "${cpu_lib}" 2>/dev/null || printf "%s\n" "${cpu_lib}")"
    gpu_lib="$(readlink -f "${gpu_lib}" 2>/dev/null || printf "%s\n" "${gpu_lib}")"

    printf "CPU libgromacs:  %s\n" "${cpu_lib}"
    printf "CUDA libgromacs: %s\n" "${gpu_lib}"

    case "${cpu_lib}" in
        /gromacs/gromacs-2026.3/lib/*|/gromacs/gromacs-2026.3/lib64/*) ;;
        *) echo "CPU executable resolved the wrong libgromacs: ${cpu_lib}" >&2; exit 1 ;;
    esac
    case "${gpu_lib}" in
        /gromacs/gromacs-2026.3-GPU/lib/*|/gromacs/gromacs-2026.3-GPU/lib64/*) ;;
        *) echo "CUDA executable resolved the wrong libgromacs: ${gpu_lib}" >&2; exit 1 ;;
    esac
'

echo
echo "=== PLUMED command interfaces ==="
"${runtime[@]}" "${cpu}" mdrun -h 2>&1 | grep -- '-plumed' | head -n 1
"${runtime[@]}" "${gpu}" mdrun -h 2>&1 | grep -- '-plumed' | head -n 1

echo
echo "=== App backing files ==="
"${runtime[@]}" bash -lc '
    set -eu
    test -s /etc/cphmd/apps.txt
    test -s /Templates/CpHMD-basic.settings
    test -s /Templates/CpHMD-advanced.settings
    test -s /analysis/extract_protonation.sh
    test -s /analysis/prepare-pdb.sh
    test -s /analysis/add-mol-database.sh
    test -s /CpHMD/scripts/make_sites
    test -d /Tutorials
    test -d /FFs/G54a7pH.ff
    test -d /FFs/CHARMM36pH.ff
    test -d /FFs/Amber14SBpH.ff
    test -d /STs/St-G54a7pH
    test -d /STs/St-CHARMM36pH
    test -d /STs/St-Amber14SBpH
    test -s /Databases/DataBaseT_G54a7pH.crg
    test -s /Databases/DataBaseT_G54a7pH.siz
    test -s /Databases/DataBaseT_CHARMM36pH.crg
    test -s /Databases/DataBaseT_CHARMM36pH.siz
    test -s /Databases/DataBaseT_Amber14SBpH.crg
    test -s /Databases/DataBaseT_Amber14SBpH.siz
    echo "App backing files are present."
'

if (( CHECK_RUNTIME_HASHES )); then
    echo
    echo "=== Embedded runtime hashes versus repository ==="

    host_files=(
        CpHMD/CpHMD-default.settings
        CpHMD/scripts/CpHMD.sh
        CpHMD/scripts/functions.sh
    )
    image_files=(
        /CpHMD/CpHMD-default.settings
        /CpHMD/scripts/CpHMD.sh
        /CpHMD/scripts/functions.sh
    )

    for index in "${!host_files[@]}"; do
        host_path="${REPO_ROOT}/${host_files[$index]}"
        image_path="${image_files[$index]}"
        [[ -f "${host_path}" ]] || {
            echo "ERROR: source file missing: ${host_path}" >&2
            exit 1
        }

        host_hash="$(sha256sum "${host_path}" | awk '{print $1}')"
        image_hash="$("${runtime[@]}" sha256sum "${image_path}" | awk '{print $1}')"
        printf '%s\n  host:  %s\n  image: %s\n' \
            "${host_files[$index]}" "${host_hash}" "${image_hash}"

        [[ "${host_hash}" == "${image_hash}" ]] || {
            echo "ERROR: embedded runtime does not match the repository source." >&2
            exit 1
        }
    done
fi

if (( RUN_APP_SMOKE )); then
    echo
    echo "=== Non-destructive application smoke tests ==="
    smoke_dir="$(mktemp -d)"
    cleanup_files+=("${smoke_dir}")

    (
        cd "${smoke_dir}"

        apptainer run --cleanenv --app get-settings "${IMAGE}" basic
        apptainer run --cleanenv --app get-settings "${IMAGE}" advanced
        apptainer run --cleanenv --app extract-protonation-script "${IMAGE}"
        apptainer run --cleanenv --app extract-databases "${IMAGE}" GROMOS
        apptainer run --cleanenv --app extract-tautomers "${IMAGE}" GROMOS
        apptainer run --cleanenv --app extract-force-fields "${IMAGE}" GROMOS

        test -s CpHMD-basic.settings
        test -s CpHMD-advanced.settings
        test -s extract_protonation.sh
        test -s Databases/DataBaseT_G54a7pH.crg
        test -s Databases/DataBaseT_G54a7pH.siz
        test -d St-G54a7pH
        test -d G54a7pH.ff
    )

    echo "Application smoke tests passed in a disposable directory."
fi

if (( GPU_DEVICE )); then
    echo
    echo "=== NVIDIA device visibility ==="
    "${runtime[@]}" nvidia-smi \
        --query-gpu=index,name,driver_version,compute_cap,memory.total \
        --format=csv,noheader

    # This confirms that the CUDA launcher starts with host driver libraries.
    # A real mdrun integration test remains necessary for release qualification.
    "${runtime[@]}" "${gpu}" mdrun -h >/dev/null
fi

echo
echo "PASS: image structure, dual GROMACS isolation, PLUMED support, apps and embedded runtime are consistent."
