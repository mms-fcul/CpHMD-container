#!/bin/bash
#
# Generate a portable SLURM wrapper for one persistent CpHMD run.
#
# The generated job:
#   - selects CPU or CUDA from export GPU=0/1 in 01_CpHMD.settings;
#   - adds --nv only for GPU=1;
#   - uses the released /CpHMD runtime embedded in the SIF by default;
#   - stages compact restart, RT history/debug, HILLS and COLVAR files;
#   - synchronizes only finalized output-contract data;
#   - automatically resubmits only after clean-boundary exit code 10.
#
# Usage:
#   generate_CpHMD_slurm.sh [scheduler-settings] [--submit|--no-submit]

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
config="${script_dir}/CpHMD-slurm.settings"
submit_override=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [scheduler-settings] [--submit|--no-submit]

Generate a SLURM job from a scheduler-variable file and the run's
01_CpHMD.settings. Generation is the default; submission is explicit.
EOF
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --submit) submit_override="1"; shift ;;
        --no-submit) submit_override="0"; shift ;;
        -h|--help) usage; exit 0 ;;
        -*) echo "ERROR: unknown option: $1" >&2; exit 1 ;;
        *) config="$1"; shift ;;
    esac
done

[ -f "${config}" ] || {
    echo "ERROR: scheduler settings file not found: ${config}" >&2
    exit 1
}

# shellcheck source=/dev/null
source "${config}" || {
    echo "ERROR: could not source scheduler settings: ${config}" >&2
    exit 1
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

canonical_path() {
    path="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath -m "${path}"
        return
    fi
    parent="$(cd "$(dirname "${path}")" 2>/dev/null && pwd)"
    [ -n "${parent}" ] || return 1
    printf '%s/%s\n' "${parent}" "$(basename "${path}")"
}

settings_export_value() {
    file="$1"
    variable="$2"
    awk -v variable="${variable}" '
        $0 ~ "^[[:space:]]*export[[:space:]]+" variable "=" {
            value=$0
            sub("^[[:space:]]*export[[:space:]]+" variable "=", "", value)
            sub(/[[:space:]]*#.*/, "", value)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            gsub(/^"/, "", value)
            gsub(/"$/, "", value)
            print value
            exit
        }
    ' "${file}"
}

time_to_seconds() {
    value="$1"
    days=0
    hours=0
    minutes=0
    seconds=0
    rest="${value}"

    case "${rest}" in
        *-*) days="${rest%%-*}"; rest="${rest#*-}" ;;
    esac

    old_ifs="${IFS}"
    IFS=':'
    set -- ${rest}
    IFS="${old_ifs}"

    case "$#" in
        1) minutes="$1" ;;
        2) minutes="$1"; seconds="$2" ;;
        3) hours="$1"; minutes="$2"; seconds="$3" ;;
        *) return 1 ;;
    esac

    case "${days}:${hours}:${minutes}:${seconds}" in
        *[!0-9:]*) return 1 ;;
    esac

    echo $((10#${days}*86400 + 10#${hours}*3600 + 10#${minutes}*60 + 10#${seconds}))
}

write_array_assignment() {
    name="$1"
    shift
    printf '%s=(\n' "${name}"
    for item in "$@"; do
        printf '    %q\n' "${item}"
    done
    printf ')\n'
}

validate_directives() {
    label="$1"
    shift
    for directive in "$@"; do
        case "${directive}" in
            --*) ;;
            *) fail "${label} entry must begin with --: ${directive}" ;;
        esac
    done
}

run_dir="${CPHMD_RUN_DIR:-}"
[ -n "${run_dir}" ] || fail "CPHMD_RUN_DIR is required."
run_dir="$(canonical_path "${run_dir}")" || fail "could not resolve CPHMD_RUN_DIR."
[ -d "${run_dir}" ] || fail "run directory not found: ${run_dir}"

settings_file="${CPHMD_SETTINGS_FILE:-}"
[ -n "${settings_file}" ] || settings_file="${run_dir}/01_CpHMD.settings"
settings_file="$(canonical_path "${settings_file}")" || fail "could not resolve CPHMD_SETTINGS_FILE."
[ -f "${settings_file}" ] || fail "CpHMD settings file not found: ${settings_file}"

image="${CPHMD_IMAGE:-}"
[ -n "${image}" ] || fail "CPHMD_IMAGE is required."
image="$(canonical_path "${image}")" || fail "could not resolve CPHMD_IMAGE."
[ -f "${image}" ] || fail "container image not found: ${image}"

gpu_value="$(settings_export_value "${settings_file}" GPU)"
ncpu_value="$(settings_export_value "${settings_file}" nCPU)"
sys_name="$(settings_export_value "${settings_file}" SysName)"

case "${gpu_value}" in
    0) mode="cpu"; partition="${CPHMD_CPU_PARTITION:-}" ;;
    1) mode="gpu"; partition="${CPHMD_GPU_PARTITION:-}" ;;
    *) fail "${settings_file} must define export GPU=0 or export GPU=1." ;;
esac

[ -n "${partition}" ] || fail "no ${mode^^} partition was configured."

case "${ncpu_value}" in
    ''|*[!0-9]*|0) fail "${settings_file} must define export nCPU as a positive integer." ;;
esac

cpus="${CPHMD_CPUS:-${ncpu_value}}"
case "${cpus}" in
    ''|*[!0-9]*|0) fail "CPHMD_CPUS must be a positive integer." ;;
esac

[ "${cpus}" = "${ncpu_value}" ] || \
    fail "SLURM CPU count (${cpus}) differs from export nCPU=${ncpu_value}."

[ -n "${sys_name}" ] || fail "${settings_file} must define export SysName."
job_name="${CPHMD_JOB_NAME:-${sys_name}}"
job_name="$(printf '%s' "${job_name}" | tr ' /' '__')"
[ -n "${job_name}" ] || fail "generated job name is empty."

walltime="${CPHMD_WALLTIME:-24:00:00}"
walltime_seconds="$(time_to_seconds "${walltime}")" || fail "invalid CPHMD_WALLTIME: ${walltime}"

signal_seconds="${CPHMD_SIGNAL_SECONDS:-900}"
stop_buffer="${CPHMD_STOP_BUFFER_SECONDS:-4500}"
max_resub="${CPHMD_MAX_RESUB:-500}"

for pair in \
    "CPHMD_SIGNAL_SECONDS:${signal_seconds}" \
    "CPHMD_STOP_BUFFER_SECONDS:${stop_buffer}" \
    "CPHMD_MAX_RESUB:${max_resub}"
do
    label="${pair%%:*}"
    value="${pair#*:}"
    case "${value}" in
        ''|*[!0-9]*) fail "${label} must be a non-negative integer." ;;
    esac
done

[ "${stop_buffer}" -gt "${signal_seconds}" ] || \
    fail "CPHMD_STOP_BUFFER_SECONDS must exceed CPHMD_SIGNAL_SECONDS."
[ "${walltime_seconds}" -gt $((stop_buffer + 60)) ] || \
    fail "walltime must exceed the stop buffer by at least 60 seconds."

max_override="${CPHMD_MAX_SECONDS_OVERRIDE:-}"
case "${max_override}" in
    ''|*[!0-9]*) [ -z "${max_override}" ] || fail "CPHMD_MAX_SECONDS_OVERRIDE must be an integer." ;;
esac

scratch_policy="${CPHMD_SCRATCH_POLICY:-prefer-local}"
case "${scratch_policy}" in
    require-local|prefer-local|allow-network) ;;
    *) fail "CPHMD_SCRATCH_POLICY must be require-local, prefer-local or allow-network." ;;
esac

slurm_dir="${CPHMD_SLURM_DIR:-${run_dir}/slurm_files}"
slurm_dir="$(canonical_path "${slurm_dir}")" || fail "could not resolve CPHMD_SLURM_DIR."
mkdir -p "${slurm_dir}" || fail "could not create SLURM directory: ${slurm_dir}"

slurm_file="${CPHMD_SLURM_FILE:-${slurm_dir}/${job_name}.slurm}"
slurm_file="$(canonical_path "${slurm_file}")" || fail "could not resolve CPHMD_SLURM_FILE."
mkdir -p "$(dirname "${slurm_file}")" || fail "could not create generated-script directory."

memory="${CPHMD_MEMORY:-}"
account="${CPHMD_ACCOUNT:-}"
scratch_root="${CPHMD_SCRATCH_ROOT:-}"
host_code_dir="${CPHMD_HOST_CODE_DIR:-}"

if [ -n "${host_code_dir}" ]; then
    host_code_dir="$(canonical_path "${host_code_dir}")" || fail "could not resolve CPHMD_HOST_CODE_DIR."
    [ -d "${host_code_dir}" ] || fail "host-code override not found: ${host_code_dir}"
fi

declare -p CPHMD_BIND_PATHS >/dev/null 2>&1 || CPHMD_BIND_PATHS=()
declare -p CPHMD_CONTAINER_ENV >/dev/null 2>&1 || CPHMD_CONTAINER_ENV=()
declare -p CPHMD_COMMON_SBATCH_DIRECTIVES >/dev/null 2>&1 || CPHMD_COMMON_SBATCH_DIRECTIVES=()
declare -p CPHMD_CPU_SBATCH_DIRECTIVES >/dev/null 2>&1 || CPHMD_CPU_SBATCH_DIRECTIVES=()
declare -p CPHMD_GPU_SBATCH_DIRECTIVES >/dev/null 2>&1 || CPHMD_GPU_SBATCH_DIRECTIVES=()

validate_directives CPHMD_COMMON_SBATCH_DIRECTIVES "${CPHMD_COMMON_SBATCH_DIRECTIVES[@]}"
validate_directives CPHMD_CPU_SBATCH_DIRECTIVES "${CPHMD_CPU_SBATCH_DIRECTIVES[@]}"
validate_directives CPHMD_GPU_SBATCH_DIRECTIVES "${CPHMD_GPU_SBATCH_DIRECTIVES[@]}"

for bind_spec in "${CPHMD_BIND_PATHS[@]}"; do
    host_path="${bind_spec%%:*}"
    [ -n "${host_path}" ] && [ "${host_path}" != "${bind_spec}" ] || \
        fail "invalid bind specification: ${bind_spec}"
    [ -e "${host_path}" ] || fail "bind source does not exist: ${host_path}"
done

for env_spec in "${CPHMD_CONTAINER_ENV[@]}"; do
    case "${env_spec}" in
        [A-Za-z_][A-Za-z0-9_]*=*) ;;
        *) fail "invalid container environment entry: ${env_spec}" ;;
    esac
done

log_file="${slurm_dir}/${job_name}_%j.log"
blockname="${sys_name}_CpHrun"

{
    printf '#!/bin/bash\n'
    printf '#SBATCH --job-name=%s\n' "${job_name}"
    [ -z "${account}" ] || printf '#SBATCH --account=%s\n' "${account}"
    printf '#SBATCH --partition=%s\n' "${partition}"
    printf '#SBATCH --ntasks=1\n'
    printf '#SBATCH --cpus-per-task=%s\n' "${cpus}"
    [ -z "${memory}" ] || printf '#SBATCH --mem=%s\n' "${memory}"
    printf '#SBATCH --time=%s\n' "${walltime}"
    printf '#SBATCH --signal=B:TERM@%s\n' "${signal_seconds}"
    printf '#SBATCH --output=%s\n' "${log_file}"

    for directive in "${CPHMD_COMMON_SBATCH_DIRECTIVES[@]}"; do
        printf '#SBATCH %s\n' "${directive}"
    done
    if [ "${mode}" = "cpu" ]; then
        for directive in "${CPHMD_CPU_SBATCH_DIRECTIVES[@]}"; do
            printf '#SBATCH %s\n' "${directive}"
        done
    else
        for directive in "${CPHMD_GPU_SBATCH_DIRECTIVES[@]}"; do
            printf '#SBATCH %s\n' "${directive}"
        done
    fi

    printf '\n'
    printf 'HOME_RUN_DIR=%q\n' "${run_dir}"
    printf 'SETTINGS_SOURCE=%q\n' "${settings_file}"
    printf 'IMAGE=%q\n' "${image}"
    printf 'MODE=%q\n' "${mode}"
    printf 'EXPECTED_GPU=%q\n' "${gpu_value}"
    printf 'EXPECTED_CPUS=%q\n' "${cpus}"
    printf 'BLOCKNAME=%q\n' "${blockname}"
    printf 'JOB_LABEL=%q\n' "${job_name}"
    printf 'REQUESTED_WALLTIME=%q\n' "${walltime}"
    printf 'STOP_BUFFER_SECONDS=%q\n' "${stop_buffer}"
    printf 'MAX_SECONDS_OVERRIDE=%q\n' "${max_override}"
    printf 'MAX_RESUB=%q\n' "${max_resub}"
    printf 'SCRATCH_ROOT_CONFIG=%q\n' "${scratch_root}"
    printf 'SCRATCH_POLICY=%q\n' "${scratch_policy}"
    printf 'HOST_CODE_DIR=%q\n' "${host_code_dir}"
    write_array_assignment BIND_SPECS "${CPHMD_BIND_PATHS[@]}"
    write_array_assignment CONTAINER_ENV_SPECS "${CPHMD_CONTAINER_ENV[@]}"

    cat <<'JOB_BODY'

settings_export_value() {
    file="$1"
    variable="$2"
    awk -v variable="${variable}" '
        $0 ~ "^[[:space:]]*export[[:space:]]+" variable "=" {
            value=$0
            sub("^[[:space:]]*export[[:space:]]+" variable "=", "", value)
            sub(/[[:space:]]*#.*/, "", value)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            gsub(/^"/, "", value)
            gsub(/"$/, "", value)
            print value
            exit
        }
    ' "${file}"
}

time_to_seconds() {
    value="$1"
    days=0; hours=0; minutes=0; seconds=0; rest="${value}"
    case "${rest}" in *-*) days="${rest%%-*}"; rest="${rest#*-}" ;; esac
    old_ifs="${IFS}"; IFS=':'; set -- ${rest}; IFS="${old_ifs}"
    case "$#" in
        1) minutes="$1" ;;
        2) minutes="$1"; seconds="$2" ;;
        3) hours="$1"; minutes="$2"; seconds="$3" ;;
        *) return 1 ;;
    esac
    case "${days}:${hours}:${minutes}:${seconds}" in *[!0-9:]*) return 1 ;; esac
    echo $((10#${days}*86400 + 10#${hours}*3600 + 10#${minutes}*60 + 10#${seconds}))
}

actual_slurm_timelimit_seconds() {
    timelimit=""
    if [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
        timelimit="$(
            scontrol show job "${SLURM_JOB_ID}" 2>/dev/null \
                | awk 'match($0,/TimeLimit=[^ ]+/){print substr($0,RSTART+10,RLENGTH-10); exit}'
        )"
    fi
    if [ -n "${timelimit}" ] && [ "${timelimit}" != "UNLIMITED" ] && \
       [ "${timelimit}" != "NOT_SET" ]; then
        time_to_seconds "${timelimit}"
    else
        time_to_seconds "${REQUESTED_WALLTIME}"
    fi
}

write_status_marker() {
    marker="$1"; shift
    marker_dir="${HOME_RUN_DIR}/resub_maint"
    marker_path="${marker_dir}/${marker}"
    marker_tmp="${marker_path}.tmp.$$"
    mkdir -p "${marker_dir}" || return 1
    {
        echo "TIMESTAMP=$(date --iso-8601=seconds 2>/dev/null || date)"
        echo "JOB_ID=${SLURM_JOB_ID:-manual}"
        echo "JOB_NAME=${SLURM_JOB_NAME:-unknown}"
        echo "HOST=$(hostname 2>/dev/null || echo unknown)"
        printf '%s\n' "$@"
    } > "${marker_tmp}" || return 1
    mv -f "${marker_tmp}" "${marker_path}"
}

persistent_checkpoint_is_complete() {
    state="${HOME_RUN_DIR}/resub_maint/checkpoint.state"
    failures=0
    [ -f "${state}" ] || return 0

    start="$(awk -F= '$1=="CPHMD_LAST_SEGMENT_START"{gsub(/[[:space:]]/,"",$2); print $2; exit}' "${state}")"
    end="$(awk -F= '$1=="CPHMD_LAST_SEGMENT_END"{gsub(/[[:space:]]/,"",$2); print $2; exit}' "${state}")"
    finalized="$(awk -F= '$1=="CPHMD_SEGMENT_FINALIZED"{gsub(/[[:space:]]/,"",$2); print $2; exit}' "${state}")"
    rt_enabled="$(awk -F= '$1=="CPHMD_REDUCED_TITRATION"{gsub(/[[:space:]]/,"",$2); print $2; exit}' "${state}")"
    rt_active_count="$(awk -F= '$1=="CPHMD_RT_ACTIVE_COUNT"{gsub(/[[:space:]]/,"",$2); print $2; exit}' "${state}")"

    case "${start}:${end}:${finalized}:${rt_enabled}" in
        *[!0-9:]*) echo "ERROR: malformed checkpoint: ${state}" >&2; return 1 ;;
    esac
    [ "${finalized}" = "1" ] || {
        echo "ERROR: checkpoint is not finalized: ${state}" >&2
        return 1
    }

    segment_prefix="${HOME_RUN_DIR}/segments/${BLOCKNAME}_cycles${start}-${end}"
    log_prefix="${HOME_RUN_DIR}/log-files/${BLOCKNAME}_cycles${start}-${end}"

    for ext in edr xtc occ mocc; do
        [ -s "${segment_prefix}.${ext}" ] || {
            echo "MISSING: ${segment_prefix}.${ext}" >&2
            failures=$((failures + 1))
        }
    done

    for ext in gro top tpr; do
        if [ -s "${segment_prefix}.${ext}" ]; then
            :
        elif [ -s "${segment_prefix}.${ext}.gz" ] && gzip -t "${segment_prefix}.${ext}.gz"; then
            :
        else
            echo "MISSING: ${segment_prefix}.${ext}[.gz]" >&2
            failures=$((failures + 1))
        fi
    done

    if [ -s "${log_prefix}.gromacs.log" ]; then
        :
    elif [ -s "${log_prefix}.gromacs.log.gz" ] && gzip -t "${log_prefix}.gromacs.log.gz"; then
        :
    else
        echo "MISSING: ${log_prefix}.gromacs.log[.gz]" >&2
        failures=$((failures + 1))
    fi
    [ -s "${log_prefix}.info" ] || {
        echo "MISSING: ${log_prefix}.info" >&2
        failures=$((failures + 1))
    }

    for ext in gro top tpr; do
        [ -s "${HOME_RUN_DIR}/restart/${BLOCKNAME}.${ext}" ] || {
            echo "MISSING: ${HOME_RUN_DIR}/restart/${BLOCKNAME}.${ext}" >&2
            failures=$((failures + 1))
        }
    done

    if [ "${rt_enabled}" = "1" ]; then
        active="${HOME_RUN_DIR}/restart/${BLOCKNAME}.RT-active.sites"
        template_occ="${HOME_RUN_DIR}/restart/${BLOCKNAME}.RT-template.occ"
        template_mocc="${HOME_RUN_DIR}/restart/${BLOCKNAME}.RT-template.mocc"

        [ -e "${active}" ] || { echo "MISSING: ${active}" >&2; failures=$((failures + 1)); }
        [ -s "${template_occ}" ] || { echo "MISSING: ${template_occ}" >&2; failures=$((failures + 1)); }
        [ -s "${template_mocc}" ] || { echo "MISSING: ${template_mocc}" >&2; failures=$((failures + 1)); }

        if [ -e "${active}" ]; then
            observed_active="$(awk 'NF{n++} END{print n+0}' "${active}")"
            [ "${observed_active}" = "${rt_active_count}" ] || {
                echo "ERROR: RT active-count mismatch: checkpoint=${rt_active_count}, file=${observed_active}" >&2
                failures=$((failures + 1))
            }
        fi
    fi

    [ "${failures}" -eq 0 ]
}

scratch_checkpoint_is_finalized() {
    state="${SCRATCH_DIR}/resub_maint/checkpoint.state"
    [ -s "${state}" ] || return 1
    finalized="$(awk -F= '$1=="CPHMD_SEGMENT_FINALIZED"{gsub(/[[:space:]]/,"",$2); print $2; exit}' "${state}")"
    [ "${finalized}" = "1" ]
}

sync_from_scratch() {
    echo "Synchronizing finalized CpHMD data to persistent storage..."
    rsync -a \
        --include='segments/***' \
        --include='restart/***' \
        --include='log-files/***' \
        --include='resub_maint/' \
        --include='resub_maint/checkpoint.state' \
        --include='*-all.sites' \
        --include='*_RT-sites.dat' \
        --include='HILLS*' \
        --include='COLVAR*' \
        --exclude='*' \
        "${SCRATCH_DIR}/" "${HOME_RUN_DIR}/"
}

sync_debug_from_scratch() {
    debug_dir="${HOME_RUN_DIR}/log-files"
    debug_log="${debug_dir}/failure_job-${SLURM_JOB_ID:-manual}.log"
    debug_tmp="${debug_log}.tmp.$$"
    mkdir -p "${debug_dir}" || return 0
    {
        echo "=== Scheduler failure summary ==="
        echo "TIMESTAMP=$(date --iso-8601=seconds 2>/dev/null || date)"
        echo "JOB_ID=${SLURM_JOB_ID:-unset}"
        echo "JOB_NAME=${SLURM_JOB_NAME:-unset}"
        echo "HOST=$(hostname 2>/dev/null || echo unknown)"
        echo "SCRATCH_DIR=${SCRATCH_DIR}"
        echo
        echo "=== Scratch inventory ==="
        find "${SCRATCH_DIR}" -maxdepth 4 -type f \
            -printf '%P\t%s bytes\n' 2>/dev/null | sort
    } > "${debug_tmp}" 2>/dev/null || {
        rm -f "${debug_tmp}" 2>/dev/null
        return 0
    }
    mv -f "${debug_tmp}" "${debug_log}" 2>/dev/null || rm -f "${debug_tmp}" 2>/dev/null
}

cleanup_scratch() {
    rm -rf "${SCRATCH_DIR}"
}

on_term() {
    echo "Caught termination signal at $(date)."
    sync_debug_from_scratch
    rm -f "${HOME_RUN_DIR}/resub_maint/RUNNING"
    write_status_marker FAILED \
        "STATUS=interrupted" \
        "EXIT_CODE=143" \
        "SCRATCH_DIR=${SCRATCH_DIR}"
    echo "Scratch preserved at ${SCRATCH_DIR}"
    exit 143
}

trap on_term TERM INT

echo "=== CpHMD SLURM job start ==="
echo "Date            : $(date)"
echo "Host            : $(hostname)"
echo "Persistent run  : ${HOME_RUN_DIR}"
echo "Container image : ${IMAGE}"
echo "Execution mode  : ${MODE}"

[ -d "${HOME_RUN_DIR}" ] || { echo "ERROR: missing run directory: ${HOME_RUN_DIR}" >&2; exit 1; }
[ -f "${SETTINGS_SOURCE}" ] || { echo "ERROR: missing settings: ${SETTINGS_SOURCE}" >&2; exit 1; }
[ -f "${IMAGE}" ] || { echo "ERROR: missing image: ${IMAGE}" >&2; exit 1; }

current_gpu="$(settings_export_value "${SETTINGS_SOURCE}" GPU)"
current_cpus="$(settings_export_value "${SETTINGS_SOURCE}" nCPU)"
[ "${current_gpu}" = "${EXPECTED_GPU}" ] || {
    echo "ERROR: GPU changed after generation: expected ${EXPECTED_GPU}, found ${current_gpu}" >&2
    exit 1
}
[ "${current_cpus}" = "${EXPECTED_CPUS}" ] || {
    echo "ERROR: nCPU changed after generation: expected ${EXPECTED_CPUS}, found ${current_cpus}" >&2
    exit 1
}

mkdir -p "${HOME_RUN_DIR}/resub_maint" "${HOME_RUN_DIR}/slurm_files" || exit 1
cd "${HOME_RUN_DIR}" || exit 1
rm -f "${HOME_RUN_DIR}/resub_maint/DONE" "${HOME_RUN_DIR}/resub_maint/FAILED"
write_status_marker RUNNING "STATUS=running" || exit 1

if ! persistent_checkpoint_is_complete; then
    rm -f "${HOME_RUN_DIR}/resub_maint/RUNNING"
    write_status_marker FAILED "STATUS=persistent_checkpoint_invalid" "EXIT_CODE=42"
    exit 42
fi

actual_limit_seconds="$(actual_slurm_timelimit_seconds)" || {
    echo "ERROR: could not determine SLURM time limit." >&2
    exit 1
}

if [ -n "${MAX_SECONDS_OVERRIDE}" ]; then
    CPHMD_MAX_SECONDS="${MAX_SECONDS_OVERRIDE}"
else
    CPHMD_MAX_SECONDS=$((actual_limit_seconds - STOP_BUFFER_SECONDS))
    [ "${CPHMD_MAX_SECONDS}" -ge 60 ] || CPHMD_MAX_SECONDS=60
fi
export CPHMD_MAX_SECONDS
export CPHMD_WALL_START_EPOCH="$(date +%s)"

counter_file="${HOME_RUN_DIR}/resub_maint/resub_count"
resubmission_count=0
[ ! -f "${counter_file}" ] || resubmission_count="$(tr -d '[:space:]' < "${counter_file}")"
case "${resubmission_count}" in
    ''|*[!0-9]*) echo "ERROR: invalid resubmission count: ${resubmission_count}" >&2; exit 9 ;;
esac

echo "Actual SLURM limit : ${actual_limit_seconds} seconds"
echo "Clean-stop guard   : ${CPHMD_MAX_SECONDS} seconds"
echo "Stop buffer        : ${STOP_BUFFER_SECONDS} seconds"
echo "Resubmissions      : ${resubmission_count}/${MAX_RESUB}"

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_RUNTIME="singularity"
else
    echo "ERROR: neither apptainer nor singularity is available." >&2
    exit 1
fi

if [ -n "${SCRATCH_ROOT_CONFIG}" ]; then
    CPHMD_SCRATCH_ROOT="${SCRATCH_ROOT_CONFIG}"
else
    CPHMD_SCRATCH_ROOT="${TMPDIR:-/tmp}"
fi
mkdir -p "${CPHMD_SCRATCH_ROOT}" || exit 1

if command -v findmnt >/dev/null 2>&1; then
    scratch_fs="$(findmnt -n -o FSTYPE -T "${CPHMD_SCRATCH_ROOT}" 2>/dev/null)"
else
    scratch_fs="$(stat -f -c %T "${CPHMD_SCRATCH_ROOT}" 2>/dev/null)"
fi
scratch_is_network=0
case "${scratch_fs}" in
    nfs|nfs4|cifs|smb*|lustre|efs|gpfs|panfs|ceph|glusterfs|fuse.sshfs) scratch_is_network=1 ;;
esac

case "${SCRATCH_POLICY}" in
    require-local)
        if [ -z "${scratch_fs}" ] || [ "${scratch_is_network}" -eq 1 ]; then
            echo "ERROR: local scratch required; detected ${scratch_fs:-unknown}" >&2
            exit 40
        fi
        ;;
    prefer-local)
        if [ -z "${scratch_fs}" ]; then
            echo "WARNING: scratch filesystem type is unknown." >&2
        elif [ "${scratch_is_network}" -eq 1 ]; then
            echo "WARNING: scratch is network-backed: ${scratch_fs}" >&2
        fi
        ;;
    allow-network) ;;
esac

SCRATCH_DIR="${CPHMD_SCRATCH_ROOT}/${USER}/CpHMD/${JOB_LABEL}_${SLURM_JOB_ID:-manual}"
export SCRATCH_DIR
mkdir -p "${SCRATCH_DIR}/restart" "${SCRATCH_DIR}/resub_maint" "${SCRATCH_DIR}/log-files" || exit 1
rsync -a "${SETTINGS_SOURCE}" "${SCRATCH_DIR}/01_CpHMD.settings" || exit 1

if [ -f "${HOME_RUN_DIR}/resub_maint/checkpoint.state" ]; then
    rsync -a "${HOME_RUN_DIR}/resub_maint/checkpoint.state" "${SCRATCH_DIR}/resub_maint/" || exit 1
fi

shopt -s nullglob
for source_file in "${HOME_RUN_DIR}/restart/"*; do
    [ -e "${source_file}" ] || continue
    rsync -a "${source_file}" "${SCRATCH_DIR}/restart/" || exit 1
done
for source_file in \
    "${HOME_RUN_DIR}/"*-all.sites \
    "${HOME_RUN_DIR}/"*_RT-sites.dat \
    "${HOME_RUN_DIR}/"HILLS* \
    "${HOME_RUN_DIR}/"COLVAR*
do
    [ -e "${source_file}" ] || continue
    rsync -a "${source_file}" "${SCRATCH_DIR}/" || exit 1
done
for source_file in "${HOME_RUN_DIR}/log-files/"*_RT-debug.pocc_RT; do
    [ -e "${source_file}" ] || continue
    rsync -a "${source_file}" "${SCRATCH_DIR}/log-files/" || exit 1
done
shopt -u nullglob

run_args=(exec --cleanenv)
[ "${MODE}" != "gpu" ] || run_args+=(--nv)
for bind_spec in "${BIND_SPECS[@]}"; do run_args+=(--bind "${bind_spec}"); done
[ -z "${HOST_CODE_DIR}" ] || run_args+=(--bind "${HOST_CODE_DIR}:/CpHMD:ro")
for env_spec in "${CONTAINER_ENV_SPECS[@]}"; do run_args+=(--env "${env_spec}"); done
run_args+=(--bind "${SCRATCH_DIR}:/out")

"${CONTAINER_RUNTIME}" "${run_args[@]}" "${IMAGE}" \
    bash -lc 'cd /out || exit 1; /CpHMD/scripts/CpHMD.sh ./01_CpHMD.settings'
cphmd_rc=$?

echo "CpHMD.sh exit code: ${cphmd_rc}"
echo "Date: $(date)"

if [ "${cphmd_rc}" -eq 0 ] || [ "${cphmd_rc}" -eq 10 ]; then
    if ! scratch_checkpoint_is_finalized; then
        echo "ERROR: CpHMD returned ${cphmd_rc} without a finalized checkpoint." >&2
        cphmd_rc=42
    fi
fi

if [ "${cphmd_rc}" -eq 0 ] || [ "${cphmd_rc}" -eq 10 ]; then
    if ! sync_from_scratch; then
        echo "ERROR: finalized synchronization failed; preserving ${SCRATCH_DIR}" >&2
        exit 20
    fi
fi

if [ "${cphmd_rc}" -eq 0 ]; then
    rm -f "${HOME_RUN_DIR}/resub_maint/RUNNING" "${HOME_RUN_DIR}/resub_maint/FAILED"
    write_status_marker DONE "STATUS=complete" "EXIT_CODE=0"
    cleanup_scratch
    echo "CpHMD run completed. Done."
    exit 0
fi

if [ "${cphmd_rc}" -eq 10 ]; then
    if [ "${resubmission_count}" -ge "${MAX_RESUB}" ]; then
        rm -f "${HOME_RUN_DIR}/resub_maint/RUNNING"
        write_status_marker FAILED \
            "STATUS=resubmission_cap_reached" \
            "EXIT_CODE=9" \
            "RESUBMISSION_COUNT=${resubmission_count}"
        exit 9
    fi

    resubmission_count=$((resubmission_count + 1))
    echo "${resubmission_count}" > "${counter_file}"
    write_status_marker RUNNING \
        "STATUS=resubmitting" \
        "RESUBMISSION_COUNT=${resubmission_count}"
    cleanup_scratch

    echo "Stopped cleanly at a validated CpHMD segment boundary."
    echo "Submitting successor job from ${BASH_SOURCE[0]}."
    if sbatch "${BASH_SOURCE[0]}"; then
        exit 0
    fi

    rm -f "${HOME_RUN_DIR}/resub_maint/RUNNING"
    write_status_marker FAILED \
        "STATUS=resubmission_failed" \
        "EXIT_CODE=21" \
        "RESUBMISSION_COUNT=${resubmission_count}"
    exit 21
fi

echo "ERROR: CpHMD failed with exit code ${cphmd_rc}; not resubmitting." >&2
sync_debug_from_scratch
rm -f "${HOME_RUN_DIR}/resub_maint/RUNNING"
write_status_marker FAILED \
    "STATUS=failed" \
    "EXIT_CODE=${cphmd_rc}" \
    "SCRATCH_DIR=${SCRATCH_DIR}"
exit "${cphmd_rc}"
JOB_BODY
} > "${slurm_file}" || fail "could not write generated SLURM file."

chmod 0755 "${slurm_file}" || fail "could not make generated SLURM file executable."
bash -n "${slurm_file}" || fail "generated SLURM script failed Bash syntax validation."

echo "Generated SLURM script:"
echo "${slurm_file}"
echo
echo "Mode      : ${mode}"
echo "Partition : ${partition}"
echo "CPUs      : ${cpus}"
echo "Walltime  : ${walltime}"
echo "Run       : ${run_dir}"
echo "Settings  : ${settings_file}"
echo "Image     : ${image}"
echo
echo "No host code is bound over /CpHMD unless CPHMD_HOST_CODE_DIR is set."

submit="${CPHMD_SUBMIT:-0}"
[ -z "${submit_override}" ] || submit="${submit_override}"

case "${submit}" in
    0)
        echo
        echo "Generation only. Review, then submit with:"
        printf 'cd %q && sbatch %q\n' "${run_dir}" "${slurm_file}"
        ;;
    1)
        echo
        echo "Submitting generated job..."
        cd "${run_dir}" || exit 1
        sbatch "${slurm_file}"
        ;;
    *)
        fail "CPHMD_SUBMIT must be 0 or 1."
        ;;
esac
