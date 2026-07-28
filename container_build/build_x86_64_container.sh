#!/bin/bash
set -euo pipefail

# Build the versioned x86-64 CpHMD release image from the current repository.
#
# Safety and reproducibility:
#   - dry-run is the default;
#   - an existing SIF is never overwritten;
#   - a dirty Git worktree is rejected unless --allow-dirty is explicit;
#   - the source state is recorded before and after the build;
#   - the SIF, checksum, definition copy and build metadata are written outside
#     the tracked repository by default.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

GROMACS_VERSION="2026.3"
CUDA_VERSION="12.6.3"
IMAGE_VERSION="0.7.0-rc1"

REPO_ROOT="${DEFAULT_REPO_ROOT}"
DEFINITION="${SCRIPT_DIR}/CpHMD-suite.def"
OUTPUT_DIR=""
OUTPUT=""
WORK_ROOT=""
TMP_DIR=""
CACHE_DIR=""
MIN_FREE_GB=40
APPLY=0
USE_SUDO=0
ALLOW_DIRTY=0
RUN_VALIDATOR=1

usage() {
    cat <<USAGE
Usage: $(basename "$0") [options]

Build the x86-64 CpHMD image with CPU-only and NVIDIA CUDA GROMACS 2026.3.
The CUDA base image is fixed by the definition to CUDA 12.6.3.

Options:
  --repo-root DIR       Repository/build-context root.
                        Default: ${DEFAULT_REPO_ROOT}
  --definition FILE     Apptainer definition file.
                        Default: ${DEFINITION}
  --image-version VER   Version used in the SIF filename.
                        Default: ${IMAGE_VERSION}
  --output-dir DIR      Output directory. Default: a sibling directory named
                        CpHMD-container-builds next to the repository.
  --output FILE         Explicit SIF path. Overrides --output-dir.
  --work-root DIR       Large build-work root; creates DIR/tmp and DIR/cache.
  --tmp-dir DIR         Explicit Apptainer temporary-build directory.
  --cache-dir DIR       Explicit user-private Apptainer OCI cache directory.
  --min-free-gb N       Required free space in the temporary filesystem.
                        Default: ${MIN_FREE_GB} GiB.
  --allow-dirty         Permit an uncommitted source tree and record that fact.
                        Release builds should normally use a clean commit.
  --skip-validator      Do not run tests/validate_dual_gromacs_image.sh.
  --apply               Perform the build. Default: validate and print only.
  --sudo                Use sudo instead of apptainer --fakeroot.
  -h, --help            Show this help.

No existing SIF is modified or replaced.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-root) REPO_ROOT="${2:?Missing repository root}"; shift 2 ;;
        --definition) DEFINITION="${2:?Missing definition}"; shift 2 ;;
        --image-version) IMAGE_VERSION="${2:?Missing image version}"; shift 2 ;;
        --output-dir) OUTPUT_DIR="${2:?Missing output directory}"; shift 2 ;;
        --output) OUTPUT="${2:?Missing output file}"; shift 2 ;;
        --work-root) WORK_ROOT="${2:?Missing work root}"; shift 2 ;;
        --tmp-dir) TMP_DIR="${2:?Missing temporary directory}"; shift 2 ;;
        --cache-dir) CACHE_DIR="${2:?Missing cache directory}"; shift 2 ;;
        --min-free-gb) MIN_FREE_GB="${2:?Missing minimum free-space value}"; shift 2 ;;
        --allow-dirty) ALLOW_DIRTY=1; shift ;;
        --skip-validator) RUN_VALIDATOR=0; shift ;;
        --apply) APPLY=1; shift ;;
        --sudo) USE_SUDO=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ "${MIN_FREE_GB}" =~ ^[0-9]+$ ]] || {
    echo "ERROR: --min-free-gb must be a non-negative integer." >&2
    exit 2
}

[[ "${IMAGE_VERSION}" =~ ^[0-9A-Za-z][0-9A-Za-z._-]*$ ]] || {
    echo "ERROR: invalid --image-version value: ${IMAGE_VERSION}" >&2
    exit 2
}

for command_name in apptainer git sha256sum awk df stat; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: ${command_name}" >&2
        exit 1
    }
done

REPO_ROOT="$(cd "${REPO_ROOT}" && pwd)"
DEFINITION="$(readlink -f "${DEFINITION}")"

[[ -d "${REPO_ROOT}/.git" ]] || {
    echo "ERROR: not a Git repository: ${REPO_ROOT}" >&2
    exit 1
}
[[ -f "${DEFINITION}" ]] || {
    echo "ERROR: definition file not found: ${DEFINITION}" >&2
    exit 1
}

required_paths=(
    CpHMD
    DelphiTools
    FFs
    STs
    Databases
    tutorials
    template
    petit1.6.1
    analysis
    plumed
)

for path in "${required_paths[@]}"; do
    [[ -e "${REPO_ROOT}/${path}" ]] || {
        echo "ERROR: build-context path is missing: ${REPO_ROOT}/${path}" >&2
        exit 1
    }
done

# Ensure the selected definition is the intended fixed toolchain definition.
grep -q '^From: nvidia/cuda:12\.6\.3-devel-ubuntu24\.04$' "${DEFINITION}" || {
    echo "ERROR: definition does not use the expected CUDA 12.6.3 build image." >&2
    exit 1
}
grep -q '^From: nvidia/cuda:12\.6\.3-runtime-ubuntu24\.04$' "${DEFINITION}" || {
    echo "ERROR: definition does not use the expected CUDA 12.6.3 runtime image." >&2
    exit 1
}
grep -q 'GROMACS version:.*2026\.3' "${DEFINITION}" || {
    echo "ERROR: definition does not contain the GROMACS 2026.3 validation." >&2
    exit 1
}

branch="$(git -C "${REPO_ROOT}" branch --show-current)"
commit="$(git -C "${REPO_ROOT}" rev-parse HEAD)"
dirty_status="$(git -C "${REPO_ROOT}" status --porcelain --untracked-files=all)"

if [[ -n "${dirty_status}" && "${ALLOW_DIRTY}" -ne 1 ]]; then
    echo "ERROR: the repository contains uncommitted or untracked changes." >&2
    echo "       Commit the exact release source, or use --allow-dirty only for a" >&2
    echo "       deliberately non-release candidate build." >&2
    git -C "${REPO_ROOT}" status --short --branch --untracked-files=all >&2
    exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
    OUTPUT_DIR="$(dirname "${REPO_ROOT}")/CpHMD-container-builds"
fi
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"

if [[ -z "${OUTPUT}" ]]; then
    OUTPUT="${OUTPUT_DIR}/CpHMD-container-v${IMAGE_VERSION}-x86_64-gromacs-${GROMACS_VERSION}-cuda-${CUDA_VERSION}.sif"
elif [[ "${OUTPUT}" != /* ]]; then
    OUTPUT="$(pwd)/${OUTPUT}"
fi
OUTPUT="$(readlink -m "${OUTPUT}")"

if [[ -e "${OUTPUT}" ]]; then
    echo "ERROR: refusing to overwrite existing output: ${OUTPUT}" >&2
    exit 1
fi

if [[ -n "${WORK_ROOT}" ]]; then
    [[ -z "${TMP_DIR}" ]] && TMP_DIR="${WORK_ROOT}/tmp"
    [[ -z "${CACHE_DIR}" ]] && CACHE_DIR="${WORK_ROOT}/cache"
fi

TMP_DIR="${TMP_DIR:-${APPTAINER_TMPDIR:-${TMPDIR:-/tmp}}}"
CACHE_DIR="${CACHE_DIR:-${APPTAINER_CACHEDIR:-${HOME}/.apptainer/cache}}"

mkdir -p "${TMP_DIR}" "${CACHE_DIR}" "$(dirname "${OUTPUT}")"
chmod 0700 "${CACHE_DIR}" 2>/dev/null || true

TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
CACHE_DIR="$(cd "${CACHE_DIR}" && pwd)"

filesystem_type() {
    local path="$1"
    if command -v findmnt >/dev/null 2>&1; then
        findmnt -n -o FSTYPE -T "${path}" 2>/dev/null || stat -f -c %T "${path}"
    else
        stat -f -c %T "${path}"
    fi
}

free_gib() {
    local path="$1"
    df -Pk "${path}" | awk 'NR==2 {printf "%d\n", $4/1024/1024}'
}

hash_runtime_sources() {
    (
        cd "${REPO_ROOT}"
        sha256sum \
            CpHMD/CpHMD-default.settings \
            CpHMD/scripts/CpHMD.sh \
            CpHMD/scripts/functions.sh \
            template/CpHMD-basic.settings \
            template/CpHMD-advanced.settings \
            analysis/extract_protonation.sh \
            analysis/prepare-pdb.sh \
            analysis/add-mol-database.sh
    )
}

tmp_fs="$(filesystem_type "${TMP_DIR}")"
tmp_free_gb="$(free_gib "${TMP_DIR}")"
output_free_gb="$(free_gib "$(dirname "${OUTPUT}")")"

case "${tmp_fs}" in
    nfs|nfs4|lustre|cifs|smb*|gpfs|panfs|ceph|glusterfs|fuse.*)
        if (( USE_SUDO == 0 )); then
            echo "WARNING: build temporary directory uses ${tmp_fs}." >&2
            echo "         Fakeroot builds are generally safer on local ext4/XFS storage." >&2
        fi
        ;;
esac

if (( APPLY == 1 && tmp_free_gb < MIN_FREE_GB )); then
    echo "ERROR: only ${tmp_free_gb} GiB are free in ${TMP_DIR}." >&2
    echo "       Select a larger filesystem with --work-root or --tmp-dir." >&2
    exit 1
fi

metadata="${OUTPUT}.build-info.txt"
checksum="${OUTPUT}.sha256"
definition_copy="${OUTPUT%.sif}.definition.def"
source_hashes_before="$(hash_runtime_sources)"
definition_hash="$(sha256sum "${DEFINITION}" | awk '{print $1}')"
diff_hash="clean"
if [[ -n "${dirty_status}" ]]; then
    diff_hash="$({
        git -C "${REPO_ROOT}" diff --binary
        while IFS= read -r -d '' untracked; do
            sha256sum "${REPO_ROOT}/${untracked}"
        done < <(git -C "${REPO_ROOT}" ls-files --others --exclude-standard -z)
    } | sha256sum | awk '{print $1}')"
fi

build_cmd=(
    env
    "APPTAINER_TMPDIR=${TMP_DIR}"
    "APPTAINER_CACHEDIR=${CACHE_DIR}"
    apptainer build --fakeroot --notest
    "${OUTPUT}" "${DEFINITION}"
)
if (( USE_SUDO )); then
    build_cmd=(
        sudo -E env
        "APPTAINER_TMPDIR=${TMP_DIR}"
        "APPTAINER_CACHEDIR=${CACHE_DIR}"
        apptainer build
        "${OUTPUT}" "${DEFINITION}"
    )
fi

printf 'Repository     : %s\n' "${REPO_ROOT}"
printf 'Branch/commit  : %s / %s\n' "${branch:-detached}" "${commit}"
printf 'Working tree   : %s\n' "$([[ -n "${dirty_status}" ]] && echo dirty || echo clean)"
printf 'Definition     : %s\n' "${DEFINITION}"
printf 'Output         : %s\n' "${OUTPUT}"
printf 'Image version  : %s\n' "${IMAGE_VERSION}"
printf 'GROMACS / CUDA : %s / %s\n' "${GROMACS_VERSION}" "${CUDA_VERSION}"
printf 'Build tmp      : %s (%s, %s GiB free)\n' "${TMP_DIR}" "${tmp_fs}" "${tmp_free_gb}"
printf 'OCI cache      : %s\n' "${CACHE_DIR}"
printf 'Output free    : %s GiB\n' "${output_free_gb}"
printf 'Mode           : %s\n' "$([[ ${APPLY} -eq 1 ]] && echo APPLY || echo DRY-RUN)"
printf 'Build command  :'
printf ' %q' "${build_cmd[@]}"
printf '\n'

if (( APPLY == 0 )); then
    echo
    echo "Dry run complete. No image, checksum, metadata or definition copy was written."
    exit 0
fi

cat > "${metadata}.tmp.$$" <<META
CpHMD container build metadata
==============================
Prepared:          $(date --iso-8601=seconds)
Repository:        ${REPO_ROOT}
Branch:            ${branch:-detached}
Commit:            ${commit}
Working tree:      $([[ -n "${dirty_status}" ]] && echo dirty || echo clean)
Dirty-state hash:  ${diff_hash}
Definition:        ${DEFINITION}
Definition SHA256: ${definition_hash}
Image version:     ${IMAGE_VERSION}
GROMACS:           ${GROMACS_VERSION}
CUDA base/runtime: ${CUDA_VERSION}
Architecture:      x86_64 / AVX2_256
Apptainer:         $(apptainer version 2>/dev/null || apptainer --version)
Build host:        $(hostname)
Build temporary:   ${TMP_DIR}
Build filesystem:  ${tmp_fs}
OCI cache:         ${CACHE_DIR}
Output:            ${OUTPUT}

Runtime source SHA256 before build
----------------------------------
${source_hashes_before}
META
mv -f "${metadata}.tmp.$$" "${metadata}"
cp -a "${DEFINITION}" "${definition_copy}"

(
    cd "${REPO_ROOT}"
    "${build_cmd[@]}"
)

source_hashes_after="$(hash_runtime_sources)"
if [[ "${source_hashes_before}" != "${source_hashes_after}" ]]; then
    echo "ERROR: release source files changed while the image was being built." >&2
    echo "       The image is retained for diagnosis but must not be released." >&2
    exit 42
fi

echo
printf 'Running the embedded %%test suite against the completed SIF:\n'
if apptainer test --cleanenv "${OUTPUT}"; then
    :
else
    test_rc=$?
    echo "ERROR: the completed SIF failed its embedded %test suite." >&2
    echo "       The image has been retained for diagnosis:" >&2
    echo "       ${OUTPUT}" >&2
    echo >&2
    echo "Re-run the test directly with:" >&2
    printf '  apptainer test --cleanenv %q\n' "${OUTPUT}" >&2
    exit "${test_rc}"
fi

sha256sum "${OUTPUT}" | tee "${checksum}"

{
    echo
    echo "Completed: $(date --iso-8601=seconds)"
    echo "Image SHA256: $(awk '{print $1}' "${checksum}")"
    echo
    echo "Runtime source SHA256 after build"
    echo "---------------------------------"
    printf '%s\n' "${source_hashes_after}"
} >> "${metadata}"

validator="${REPO_ROOT}/tests/validate_dual_gromacs_image.sh"
if (( RUN_VALIDATOR == 1 )); then
    [[ -x "${validator}" ]] || {
        echo "ERROR: validator is missing or not executable: ${validator}" >&2
        exit 1
    }
    "${validator}" --image "${OUTPUT}" --repo-root "${REPO_ROOT}"
fi

echo
echo "Image built and CPU-side validation completed successfully."
echo "A real GPU integration test with --nv is still required before release."
printf 'Image:      %s\n' "${OUTPUT}"
printf 'Checksum:   %s\n' "${checksum}"
printf 'Metadata:   %s\n' "${metadata}"
printf 'Definition: %s\n' "${definition_copy}"
