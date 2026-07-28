#!/bin/bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: diagnose_cphmd_image_test.sh IMAGE.sif

Runs the same named checks as the image %test section, but from a host script.
This is useful when diagnosing an already-built candidate image.
USAGE
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 2
fi

image="$1"
[[ -f "$image" ]] || {
    echo "ERROR: image not found: $image" >&2
    exit 1
}

if command -v apptainer >/dev/null 2>&1; then
    runtime=apptainer
elif command -v singularity >/dev/null 2>&1; then
    runtime=singularity
else
    echo "ERROR: neither apptainer nor singularity is available." >&2
    exit 1
fi

"$runtime" exec --cleanenv "$image" /bin/sh -s <<'INNER'
set -eu

test_log="$(mktemp "${TMPDIR:-/tmp}/cphmd-image-test.XXXXXX")"
trap 'rm -f "${test_log}"' EXIT HUP INT TERM

require_output () {
    test_label="$1"
    test_pattern="$2"
    shift 2
    printf '[TEST] %s\n' "${test_label}"
    if "$@" >"${test_log}" 2>&1; then :; else
        test_rc=$?
        printf '[FAIL] %s: command exited non-zero\n' "${test_label}" >&2
        cat "${test_log}" >&2 || true
        exit "${test_rc}"
    fi
    if ! grep -Eq -- "${test_pattern}" "${test_log}"; then
        printf '[FAIL] %s: expected pattern not found: %s\n' "${test_label}" "${test_pattern}" >&2
        cat "${test_log}" >&2 || true
        exit 1
    fi
    printf '[PASS] %s\n' "${test_label}"
}

require_command () {
    test_label="$1"
    shift
    printf '[TEST] %s\n' "${test_label}"
    if "$@" >"${test_log}" 2>&1; then :; else
        test_rc=$?
        printf '[FAIL] %s: command exited non-zero\n' "${test_label}" >&2
        cat "${test_log}" >&2 || true
        exit "${test_rc}"
    fi
    printf '[PASS] %s\n' "${test_label}"
}

require_nonempty_file () {
    test_label="$1"; test_path="$2"
    printf '[TEST] %s\n' "${test_label}"
    [ -s "${test_path}" ] || { printf '[FAIL] %s: missing or empty file: %s\n' "${test_label}" "${test_path}" >&2; exit 1; }
    printf '[PASS] %s\n' "${test_label}"
}

require_executable () {
    test_label="$1"; test_path="$2"
    printf '[TEST] %s\n' "${test_label}"
    [ -x "${test_path}" ] || { printf '[FAIL] %s: missing or non-executable file: %s\n' "${test_label}" "${test_path}" >&2; ls -l "${test_path}" >&2 2>/dev/null || true; exit 1; }
    printf '[PASS] %s\n' "${test_label}"
}

require_directory () {
    test_label="$1"; test_path="$2"
    printf '[TEST] %s\n' "${test_label}"
    [ -d "${test_path}" ] || { printf '[FAIL] %s: missing directory: %s\n' "${test_label}" "${test_path}" >&2; exit 1; }
    printf '[PASS] %s\n' "${test_label}"
}

require_libgromacs_root () {
    test_label="$1"
    executable="$2"
    expected_root="$3"
    library_path="${expected_root}/lib:${expected_root}/lib64:/opt/plumed/lib:/opt/plumed/lib64"
    case "${expected_root}" in
        *-GPU) library_path="${library_path}:/usr/local/cuda/lib64" ;;
    esac

    printf '[TEST] %s\n' "${test_label}"
    if env LD_LIBRARY_PATH="${library_path}" ldd "${executable}" >"${test_log}" 2>&1; then
        :
    else
        test_rc=$?
        printf '[FAIL] %s: ldd exited non-zero\n' "${test_label}" >&2
        cat "${test_log}" >&2 || true
        exit "${test_rc}"
    fi

    lib_path="$(awk '$1 ~ /^libgromacs([.]so)?/ { print $3; exit }' "${test_log}")"
    [ -n "${lib_path}" ] || {
        printf '[FAIL] %s: libgromacs was not found in ldd output\n' "${test_label}" >&2
        cat "${test_log}" >&2 || true
        exit 1
    }

    canonical_path="$(readlink -f "${lib_path}" 2>/dev/null || printf '%s\n' "${lib_path}")"
    printf '       resolved libgromacs: %s\n' "${canonical_path}"
    case "${canonical_path}" in
        "${expected_root}"/lib/*|"${expected_root}"/lib64/*) ;;
        *)
            printf '[FAIL] %s: expected libgromacs below %s, observed %s\n' \
                "${test_label}" "${expected_root}" "${canonical_path}" >&2
            cat "${test_log}" >&2 || true
            exit 1
            ;;
    esac
    printf '[PASS] %s\n' "${test_label}"
}

cpu=/gromacs/bin/gmx
gpu=/gromacs/bin/gmx-gpu
require_output 'CPU GROMACS version' 'GROMACS version:[[:space:]]*2026\.3' "${cpu}" --version
require_output 'CPU build has GPU support disabled' 'GPU support:[[:space:]]*disabled' "${cpu}" --version
require_output 'CPU build uses AVX2_256' 'SIMD instructions:[[:space:]]*AVX2_256' "${cpu}" --version
require_output 'CUDA GROMACS version' 'GROMACS version:[[:space:]]*2026\.3' "${gpu}" --version
require_output 'CUDA build reports CUDA support' 'GPU support:[[:space:]]*CUDA' "${gpu}" --version
require_output 'CUDA build uses AVX2_256' 'SIMD instructions:[[:space:]]*AVX2_256' "${gpu}" --version
require_output 'CPU mdrun exposes PLUMED option' '-plumed' "${cpu}" mdrun -h
require_output 'CUDA mdrun exposes PLUMED option' '-plumed' "${gpu}" mdrun -h
require_libgromacs_root 'CPU executable resolves CPU libgromacs' /gromacs/gromacs-2026.3/bin/gmx /gromacs/gromacs-2026.3
require_libgromacs_root 'CUDA executable resolves CUDA libgromacs' /gromacs/gromacs-2026.3-GPU/bin/gmx /gromacs/gromacs-2026.3-GPU
require_command 'CpHMD.sh syntax' bash -n /CpHMD/scripts/CpHMD.sh
require_command 'functions.sh syntax' bash -n /CpHMD/scripts/functions.sh
require_command 'extract_protonation.sh syntax' bash -n /analysis/extract_protonation.sh
require_command 'prepare-pdb.sh syntax' bash -n /analysis/prepare-pdb.sh
require_command 'add-mol-database.sh syntax' bash -n /analysis/add-mol-database.sh
require_nonempty_file 'default CpHMD settings' /CpHMD/CpHMD-default.settings
require_nonempty_file 'basic settings template' /Templates/CpHMD-basic.settings
require_nonempty_file 'advanced settings template' /Templates/CpHMD-advanced.settings
require_executable 'PETIT executable' /programs/petit1.6.1/petit
require_executable 'DelPhi executable' /Delphi/delphiT
require_nonempty_file 'GROMOS charge database' /Databases/DataBaseT_G54a7pH.crg
require_nonempty_file 'CHARMM charge database' /Databases/DataBaseT_CHARMM36pH.crg
require_nonempty_file 'AMBER charge database' /Databases/DataBaseT_Amber14SBpH.crg
require_directory 'GROMOS tautomer directory' /STs/St-G54a7pH
require_directory 'CHARMM tautomer directory' /STs/St-CHARMM36pH
require_directory 'AMBER tautomer directory' /STs/St-Amber14SBpH
require_nonempty_file 'app registry' /etc/cphmd/apps.txt
require_command 'gzip availability' gzip --version
printf '[PASS] all diagnostic image tests completed successfully\n'
INNER
