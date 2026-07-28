# CpHMD function library
#
# This file is sourced by CpHMD.sh after the user settings have been read.
# It contains the complete implementation of the CpHMD workflow:
#
#   1. runtime and GROMACS/GPU preflight checks;
#   2. input, topology and force-field validation;
#   3. auxiliary coordinate/site preparation;
#   4. Delphi PB and Monte Carlo titration;
#   5. topology/protonation-state updates;
#   6. relaxation and effective molecular dynamics;
#   7. transactional trajectory/energy accumulation;
#   8. checkpoint, validation, compression and cleanup.
#
# The functions rely on variables exported by CpHMD.sh and the settings file.
# Functions that form the persistent-output contract return exit code 42 on
# missing, incomplete or invalid MD/output files. PB/MC validation failures
# return exit code 43. The scheduler wrapper resubmits neither failure code.
#
# Do not execute this file directly; source it from CpHMD.sh.

# -----------------------------------------------------------------------------
# CpHMD runtime/output safety helpers
# -----------------------------------------------------------------------------
# These helpers deliberately use inexpensive file-size and command-status checks
# during each CpHMD cycle. Full binary scans with `gmx check` are reserved for a
# scheduler-segment boundary, where their cost is negligible compared with the
# MD/PB-MC work completed in the segment.

# cphmd_output_error MESSAGE...
#   Print an output-contract error and return CPHMD_EXIT_OUTPUT_FAILURE (42).
#   Unlike `message E`, this function returns to its caller so the failure code
#   can propagate through functions.sh, CpHMD.sh and the scheduler wrapper.
cphmd_output_error ()
{
    # Return a distinctive application error without calling message E, because
    # message E always exits with status 1. CpHMD.sh propagates this status to the
    # scheduler wrapper, which must not resubmit it.
    echo "${prog:-CpHMD.sh}: Error: $*" >&2
    return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
}

# cphmd_require_nonempty_files LABEL FILE...
#   Verify that every listed file exists and has non-zero size. The label is
#   included in the diagnostic. Used after GROMACS/PBMC commands and before any
#   checkpoint is written, preventing stale or partial files from advancing the
#   simulation state.
cphmd_require_nonempty_files ()
{
    local label="$1"
    shift
    local missing=()
    local file

    for file in "$@"; do
        [[ -s "${file}" ]] || missing+=("${file}")
    done

    if (( ${#missing[@]} > 0 )); then
        echo "${prog:-CpHMD.sh}: Error: ${label}; missing or empty file(s):" >&2
        printf '  %s\n' "${missing[@]}" >&2
        return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
    fi
}

# cphmd_copy_atomic SOURCE DESTINATION
#   Copy one persistent output transactionally: write a temporary file in the
#   destination directory, then rename it into place. The final name is never
#   exposed until the copy succeeds completely.
cphmd_copy_atomic ()
{
    # Copy to a temporary file in the destination directory and rename only
    # after cp succeeds. The destination parent is created first so that the
    # organized output layout can be used from a fresh run directory.
    local source="$1"
    local destination="$2"
    local destination_dir
    local temporary

    destination_dir="$(dirname -- "${destination}")"
    mkdir -p -- "${destination_dir}" || {
        cphmd_output_error "Could not create output directory ${destination_dir}."
        return $?
    }

    temporary="${destination}.tmp.$$"
    rm -f -- "${temporary}"

    cp -f -- "${source}" "${temporary}" || {
        rm -f -- "${temporary}"
        cphmd_output_error "Could not copy ${source} to ${destination}."
        return $?
    }

    mv -f -- "${temporary}" "${destination}" || {
        rm -f -- "${temporary}"
        cphmd_output_error "Could not finalize ${destination}."
        return $?
    }
}

# cphmd_write_state_atomic FINISHED_CYCLE FINALIZED SOURCE_LABEL
#   Write resub_maint/checkpoint.state through a temporary file and atomic
#   rename. In reduced-titration mode the state also records the last full-site
#   refresh and the size of the active RT subset, allowing strict restart checks.
cphmd_write_state_atomic ()
{
    local finished_cycle="$1"
    local finalized="$2"
    local source_label="$3"
    local segment_start="${CPHMD_SEGMENT_START_CYCLE:-1}"
    local state_path="../${CPHMD_CHECKPOINT_FILE:-resub_maint/checkpoint.state}"
    local state_dir
    local state_tmp
    local rt_enabled="${ReduceTitration:-0}"
    local rt_last_full="${CPHMD_RT_LAST_FULL_CYCLE:-0}"
    local rt_active_count=0

    state_dir="$(dirname -- "${state_path}")"
    mkdir -p -- "${state_dir}" || {
        cphmd_output_error "Could not create checkpoint directory ${state_dir}."
        return $?
    }

    if (( rt_enabled == 1 )); then
        [[ -e "${CPHMD_RT_ACTIVE_SITES_FILE}" ]] || {
            cphmd_output_error "Reduced titration is enabled, but ${CPHMD_RT_ACTIVE_SITES_FILE} is missing."
            return $?
        }
        rt_active_count="$(awk 'NF{n++} END{print n+0}' "${CPHMD_RT_ACTIVE_SITES_FILE}")"
    fi

    state_tmp="${state_path}.tmp.$$"

    cat > "${state_tmp}" <<EOF
CPHMD_DONE_CYCLES=${finished_cycle}
CPHMD_LAST_SEGMENT_START=${segment_start}
CPHMD_LAST_SEGMENT_END=${finished_cycle}
CPHMD_SEGMENT_FINALIZED=${finalized}
CPHMD_STATE_SOURCE=${source_label}
CPHMD_REDUCED_TITRATION=${rt_enabled}
CPHMD_RT_LAST_FULL_CYCLE=${rt_last_full}
CPHMD_RT_ACTIVE_COUNT=${rt_active_count}
EOF

    mv -f -- "${state_tmp}" "${state_path}" || {
        rm -f -- "${state_tmp}"
        cphmd_output_error "Could not finalize checkpoint state ${state_path}."
        return $?
    }
}

# cphmd_gzip_atomic FILE
#   Compress an immutable finalized file using CPHMD_GZIP_LEVEL. The .gz archive
#   is created under a temporary name, verified with `gzip -t`, and only then
#   replaces the uncompressed source.
cphmd_gzip_atomic ()
{
    # Compress one finalized file without exposing a partial .gz. The original is
    # removed only after gzip -t verifies the completed archive.
    local source="$1"
    local level="${CPHMD_GZIP_LEVEL:-6}"
    local destination="${source}.gz"
    local temporary="${destination}.tmp.$$"

    [[ -s "${source}" ]] || {
        cphmd_output_error "Cannot compress missing or empty file ${source}."
        return $?
    }

    if ! [[ "${level}" =~ ^[1-9]$ ]]; then
        cphmd_output_error "CPHMD_GZIP_LEVEL must be an integer from 1 to 9; got '${level}'."
        return $?
    fi

    rm -f "${temporary}"
    if ! gzip -c -"${level}" -- "${source}" > "${temporary}"; then
        rm -f "${temporary}"
        cphmd_output_error "gzip failed for ${source}."
        return $?
    fi

    if ! gzip -t -- "${temporary}"; then
        rm -f "${temporary}"
        cphmd_output_error "gzip verification failed for ${source}."
        return $?
    fi

    mv -f -- "${temporary}" "${destination}" || {
        rm -f "${temporary}"
        cphmd_output_error "Could not finalize ${destination}."
        return $?
    }
    rm -f -- "${source}"
}

# cphmd_preserve_failure_diagnostics [CYCLE]
#   Write one flat text diagnostic under log-files/. The persistent log
#   directory must contain files only, never recursive copies of the scratch
#   layout. The node-local scratch tree is preserved separately by the scheduler
#   wrapper and remains the authoritative source for deeper debugging.
cphmd_preserve_failure_diagnostics ()
{
    # This function is best-effort and must never hide the original simulation
    # error. All diagnostic reads/writes are therefore allowed to fail quietly.
    local cycle="${1:-unknown}"
    local log_dir="../${CPHMD_LOG_DIR:-log-files}"
    local diagnostic_log="${log_dir}/${blockname}_failure_cycle${cycle}_pid$$.log"
    local temporary="${diagnostic_log}.tmp.$$"

    mkdir -p -- "${log_dir}" 2>/dev/null || return 0

    {
        echo "=== CpHMD failure diagnostic ==="
        echo "TIMESTAMP=$(date -Is 2>/dev/null || date)"
        echo "CYCLE=${cycle}"
        echo "PID=$$"
        echo "JOB_ID=${SLURM_JOB_ID:-unset}"
        echo "JOB_NAME=${SLURM_JOB_NAME:-unset}"
        echo "HOST=$(hostname 2>/dev/null || echo unknown)"
        echo "WORKING_DIRECTORY=$(pwd -P 2>/dev/null || pwd)"
        echo

        if [[ -f "${blockname}_MD.log" ]]; then
            echo "=== Driver log: ${blockname}_MD.log ==="
            cat -- "${blockname}_MD.log" 2>/dev/null || true
            echo
        fi

        if [[ -f "${blockname}.info" ]]; then
            echo "=== Run information: ${blockname}.info ==="
            cat -- "${blockname}.info" 2>/dev/null || true
            echo
        fi

        echo "=== Relevant temporary files ==="
        find . -maxdepth 1 -type f \
            \( -name 'TMP_*.log' -o -name 'TMP_*.out' \
               -o -name 'TMP_*.err' -o -name 'TMP_*.mdp' \
               -o -name 'TMP_*.tpr' \) \
            -printf '%f\t%s bytes\n' 2>/dev/null \
          | sort || true
    } > "${temporary}" 2>/dev/null || {
        rm -f -- "${temporary}" 2>/dev/null || true
        return 0
    }

    mv -f -- "${temporary}" "${diagnostic_log}" 2>/dev/null || \
        rm -f -- "${temporary}" 2>/dev/null || true
}


# cphmd_pbmc_error STAGE MESSAGE...
#   Preserve Delphi/convert/PETIT diagnostics in one flat log file and return
#   the dedicated PB/MC failure code (43). The scheduler will preserve scratch
#   and will not resubmit the failed segment.
cphmd_pbmc_error ()
{
    local stage="$1"
    shift
    local log_dir="../${CPHMD_LOG_DIR:-log-files}"
    local diagnostic="${log_dir}/${blockname}_PBMC-failure_cycle${Cycle:-pre}_stage-${stage}_pid$$.log"
    local temporary="${diagnostic}.tmp.$$"
    local file

    mkdir -p -- "${log_dir}" 2>/dev/null || true

    {
        echo "=== CpHMD PB/MC failure ==="
        echo "TIMESTAMP=$(date -Is 2>/dev/null || date)"
        echo "CYCLE=${Cycle:-unset}"
        echo "STAGE=${stage}"
        echo "JOB_ID=${SLURM_JOB_ID:-unset}"
        echo "JOB_NAME=${SLURM_JOB_NAME:-unset}"
        echo "HOST=$(hostname 2>/dev/null || echo unknown)"
        echo "MESSAGE=$*"
        echo

        for file in TMP_delphi.out TMP_delphi.err TMP_convert.err TMP_MCarlo.out TMP_MCarlo.err; do
            if [[ -f "${file}" ]]; then
                echo "=== ${file} ==="
                cat -- "${file}" 2>/dev/null || true
                echo
            fi
        done

        echo "=== PB/MC file inventory ==="
        find . -maxdepth 1 -type f \
            \( -name 'PBMC*' -o -name 'TMP_PBMC*' -o -name 'TMP_MCarlo*' \
               -o -name 'TMP_delphi*' -o -name 'TMP_convert*' \) \
            -printf '%f\t%s bytes\n' 2>/dev/null | sort || true
    } > "${temporary}" 2>/dev/null || true

    if [[ -s "${temporary}" ]]; then
        mv -f -- "${temporary}" "${diagnostic}" 2>/dev/null || true
    else
        rm -f -- "${temporary}" 2>/dev/null || true
    fi

    echo "${prog:-CpHMD.sh}: Error: PB/MC ${stage} failed: $*" >&2
    echo "${prog:-CpHMD.sh}: Error: PB/MC diagnostic: ${diagnostic}" >&2
    return "${CPHMD_EXIT_PBMC_FAILURE:-43}"
}

# cphmd_count_nonblank FILE
#   Print the number of nonblank lines in FILE. Missing files are errors.
cphmd_count_nonblank ()
{
    local file="$1"
    [[ -e "${file}" ]] || return 1
    awk 'NF{n++} END{print n+0}' "${file}"
}

# cphmd_rt_full_refresh_cycle CYCLE
#   Return success for cycles 1, 1+RTInterval, 1+2*RTInterval, ... . Using
#   (cycle-1) modulo interval also handles RTInterval=1 correctly.
cphmd_rt_full_refresh_cycle ()
{
    local cycle="$1"
    (( (cycle - 1) % RTInterval == 0 ))
}

# cphmd_validate_petit_output FILE EXPECTED_SITES
#   Require one complete PETIT state line plus one header and population record
#   for every site. This prevents an empty/partial Monte Carlo result from being
#   mistaken for a legitimate zero-site reduced-titration selection.
cphmd_validate_petit_output ()
{
    local file="$1"
    local expected="$2"

    [[ -s "${file}" ]] || {
        cphmd_pbmc_error petit-output "PETIT produced no nonempty output file ${file}."
        return $?
    }

    awk -v expected="${expected}" '
        function numeric(x) {
            return x ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
        }
        $1 == "f" {
            f_lines++
            f_sites = NF - 1
            for (i=2; i<=NF; i++) if ($i !~ /^[0-9]+$/) bad_state=1
        }
        $1 == ">" {
            if ($3 !~ /^[0-9]+$/) bad_header=1
            else header[$3+1]++
        }
        substr($0,1,1) == "." && $0 !~ /tot/ {
            if ($4 !~ /^[0-9]+$/ || !numeric($5)) bad_population=1
            else population[$4+1]++
            for (i=6; i<=NF; i++) if (!numeric($i)) bad_population=1
        }
        END {
            if (f_lines != 1 || f_sites != expected || bad_state || bad_header || bad_population) exit 1
            for (i=1; i<=expected; i++) {
                if (header[i] != 1 || population[i] != 1) exit 1
            }
        }
    ' "${file}" || {
        cphmd_pbmc_error petit-output \
            "PETIT output ${file} does not contain a complete ${expected}-site state/population result."
        return $?
    }
}

# cphmd_refresh_reduced_titration_state
#   After a validated all-site PETIT calculation, create four files atomically:
#     * current active RT subset;
#     * full-length fixed-state template;
#     * full-length fixed-population template;
#     * full-length states/populations actually applied in this cycle.
#   Fixed sites use their most populated state. Active sites use the sampled
#   PETIT state. NTR/CTR sites are always retained in the active subset.
cphmd_refresh_reduced_titration_state ()
{
    local all_sites="${CPHMD_ALL_SITES_FILE}"
    local active_sites="${CPHMD_RT_ACTIVE_SITES_FILE}"
    local template_occ="${CPHMD_RT_TEMPLATE_OCC_FILE}"
    local template_mocc="${CPHMD_RT_TEMPLATE_MOCC_FILE}"
    local current_occ="${CPHMD_RT_CURRENT_OCC_FILE}"
    local current_mocc="${CPHMD_RT_CURRENT_MOCC_FILE}"
    local total
    local active_tmp="${active_sites}.tmp.$$"
    local template_occ_tmp="${template_occ}.tmp.$$"
    local template_mocc_tmp="${template_mocc}.tmp.$$"
    local current_occ_tmp="${current_occ}.tmp.$$"
    local current_mocc_tmp="${current_mocc}.tmp.$$"
    local summary
    local rc
    local active_count

    total="$(cphmd_count_nonblank "${all_sites}")" || {
        cphmd_pbmc_error reduced-selection "All-site definition ${all_sites} is missing."
        return $?
    }

    cphmd_validate_petit_output TMP_MCarlo.out "${total}" || return $?

    rm -f -- "${active_tmp}" "${template_occ_tmp}" "${template_mocc_tmp}" \
        "${current_occ_tmp}" "${current_mocc_tmp}"
    : > "${active_tmp}"
    : > "${template_occ_tmp}"
    : > "${template_mocc_tmp}"
    : > "${current_occ_tmp}"
    : > "${current_mocc_tmp}"

    summary="$(gawk \
        -v threshold="${RTThreshold}" \
        -v allsites="${all_sites}" \
        -v activefile="${active_tmp}" \
        -v templateocc="${template_occ_tmp}" \
        -v templatemocc="${template_mocc_tmp}" \
        -v currentocc="${current_occ_tmp}" \
        -v currentmocc="${current_mocc_tmp}" '
        function numeric(x) {
            return x ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
        }
        BEGIN {
            while ((getline line < allsites) > 0) {
                if (line ~ /[^[:space:]]/) site[++total] = line
            }
            close(allsites)

            while ((getline < "TMP_MCarlo.out") > 0) {
                if ($1 == ">") {
                    idx = $3 + 1
                    if ($2 ~ /NTR/ || $2 ~ /CTR/) terminal[idx] = 1
                }
                if ($1 == "f") {
                    fcount++
                    nstates = NF - 1
                    for (i=2; i<=NF; i++) sampled[i-1] = $i
                }
                if (substr($0,1,1) == "." && $0 !~ /tot/) {
                    idx = $4 + 1
                    mocc[idx] = $5
                    maxp = -1
                    maxi = -1
                    for (i=6; i<=NF; i++) {
                        if (!numeric($i)) exit 20
                        if ($i > maxp) {
                            maxp = $i
                            maxi = i - 6
                        }
                    }
                    maxprob[idx] = maxp
                    maxstate[idx] = maxi
                }
            }
            close("TMP_MCarlo.out")

            if (fcount != 1 || nstates != total) exit 21

            active = 0
            for (i=1; i<=total; i++) {
                if (!(i in sampled) || !(i in mocc) || !(i in maxstate)) exit 22
                keep = terminal[i] || maxprob[i] <= (1-threshold)

                if (keep) {
                    print site[i] > activefile
                    printf "- " > templateocc
                    printf "- " > templatemocc
                    printf "%d ", sampled[i] > currentocc
                    active++
                } else {
                    printf "%d ", maxstate[i] > templateocc
                    printf "%.10g ", mocc[i] > templatemocc
                    printf "%d ", maxstate[i] > currentocc
                }
                printf "%.10g ", mocc[i] > currentmocc
            }
            printf "\n" > templateocc
            printf "\n" > templatemocc
            printf "\n" > currentocc
            printf "\n" > currentmocc

            close(activefile); close(templateocc); close(templatemocc)
            close(currentocc); close(currentmocc)
            printf "%d %d\n", active, total
        }
    ')"
    rc=$?

    if (( rc != 0 )); then
        rm -f -- "${active_tmp}" "${template_occ_tmp}" "${template_mocc_tmp}" \
            "${current_occ_tmp}" "${current_mocc_tmp}"
        cphmd_pbmc_error reduced-selection \
            "Could not derive a complete reduced-titration state from validated PETIT output (gawk exit ${rc})."
        return $?
    fi

    read -r active_count parsed_total <<< "${summary}"
    if [[ "${parsed_total}" != "${total}" || ! "${active_count}" =~ ^[0-9]+$ ]]; then
        rm -f -- "${active_tmp}" "${template_occ_tmp}" "${template_mocc_tmp}" \
            "${current_occ_tmp}" "${current_mocc_tmp}"
        cphmd_pbmc_error reduced-selection "Reduced-titration selection returned an invalid summary: ${summary}."
        return $?
    fi

    mv -f -- "${active_tmp}" "${active_sites}" || return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
    mv -f -- "${template_occ_tmp}" "${template_occ}" || return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
    mv -f -- "${template_mocc_tmp}" "${template_mocc}" || return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
    mv -f -- "${current_occ_tmp}" "${current_occ}" || return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
    mv -f -- "${current_mocc_tmp}" "${current_mocc}" || return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"

    export CPHMD_RT_LAST_FULL_CYCLE="${Cycle}"
    export CPHMD_RT_ACTIVE_COUNT="${active_count}"

    {
        echo "Active reduced-titration sites after full update at Cycle = ${Cycle}"
        echo "Active sites: ${active_count}/${total}"
        cat -- "${active_sites}"
        echo
    } >> "${CPHMD_RT_HISTORY_FILE}"

    message W "Reduced-titration refresh at cycle ${Cycle}: ${active_count}/${total} sites remain active."
}

# cphmd_merge_reduced_titration_cycle
#   Merge a validated PETIT result for only the active subset into the full-size
#   fixed templates. The generated current OCC/MOCC lines always contain one
#   entry per all-site definition and are used for both topology and output.
cphmd_merge_reduced_titration_cycle ()
{
    local active_count
    local all_count
    local current_occ_tmp="${CPHMD_RT_CURRENT_OCC_FILE}.tmp.$$"
    local current_mocc_tmp="${CPHMD_RT_CURRENT_MOCC_FILE}.tmp.$$"
    local rc

    active_count="$(cphmd_count_nonblank "${CPHMD_RT_ACTIVE_SITES_FILE}")" || return 42
    all_count="$(cphmd_count_nonblank "${CPHMD_ALL_SITES_FILE}")" || return 42

    (( active_count > 0 )) || {
        cphmd_pbmc_error reduced-merge "Reduced PETIT merge requested with zero active sites."
        return $?
    }

    cphmd_validate_petit_output TMP_MCarlo.out "${active_count}" || return $?
    cphmd_require_nonempty_files "reduced-titration templates are incomplete" \
        "${CPHMD_RT_TEMPLATE_OCC_FILE}" "${CPHMD_RT_TEMPLATE_MOCC_FILE}" || return $?

    rm -f -- "${current_occ_tmp}" "${current_mocc_tmp}"

    gawk \
        -v expected_all="${all_count}" \
        -v expected_active="${active_count}" \
        -v templateocc="${CPHMD_RT_TEMPLATE_OCC_FILE}" \
        -v templatemocc="${CPHMD_RT_TEMPLATE_MOCC_FILE}" \
        -v outocc="${current_occ_tmp}" \
        -v outmocc="${current_mocc_tmp}" '
        BEGIN {
            if ((getline line < templateocc) <= 0) exit 30
            nt = split(line, tocc)
            close(templateocc)
            if ((getline line < templatemocc) <= 0) exit 31
            nm = split(line, tmocc)
            close(templatemocc)
            if (nt != expected_all || nm != expected_all) exit 32

            while ((getline < "TMP_MCarlo.out") > 0) {
                if ($1 == "f") {
                    nf = NF - 1
                    for (i=2; i<=NF; i++) activeocc[i-1] = $i
                }
                if (substr($0,1,1) == "." && $0 !~ /tot/) activemocc[$4+1] = $5
            }
            close("TMP_MCarlo.out")
            if (nf != expected_active) exit 33

            c = 0
            for (i=1; i<=expected_all; i++) {
                if (tocc[i] == "-") {
                    c++
                    if (!(c in activeocc) || !(c in activemocc)) exit 34
                    printf "%d ", activeocc[c] > outocc
                    printf "%.10g ", activemocc[c] > outmocc
                } else {
                    printf "%d ", tocc[i] > outocc
                    printf "%.10g ", tmocc[i] > outmocc
                }
            }
            if (c != expected_active) exit 35
            printf "\n" > outocc
            printf "\n" > outmocc
            close(outocc); close(outmocc)
        }
    '
    rc=$?

    if (( rc != 0 )); then
        rm -f -- "${current_occ_tmp}" "${current_mocc_tmp}"
        cphmd_pbmc_error reduced-merge \
            "Could not merge active PETIT states into the full reduced-titration templates (gawk exit ${rc})."
        return $?
    fi

    mv -f -- "${current_occ_tmp}" "${CPHMD_RT_CURRENT_OCC_FILE}" || return 42
    mv -f -- "${current_mocc_tmp}" "${CPHMD_RT_CURRENT_MOCC_FILE}" || return 42
}

# cphmd_materialize_fixed_rt_cycle
#   Build the current full state/population lines when the active RT subset is
#   legitimately empty. No PB/MC command is needed; all sites remain fixed in
#   the states selected at the last full-site refresh.
cphmd_materialize_fixed_rt_cycle ()
{
    local all_count
    local active_count
    local occ_tokens
    local mocc_tokens

    all_count="$(cphmd_count_nonblank "${CPHMD_ALL_SITES_FILE}")" || return 42
    active_count="$(cphmd_count_nonblank "${CPHMD_RT_ACTIVE_SITES_FILE}")" || return 42
    (( active_count == 0 )) || return 42

    cphmd_require_nonempty_files "fixed reduced-titration templates are incomplete" \
        "${CPHMD_RT_TEMPLATE_OCC_FILE}" "${CPHMD_RT_TEMPLATE_MOCC_FILE}" || return $?

    if grep -qw -- '-' "${CPHMD_RT_TEMPLATE_OCC_FILE}" || \
       grep -qw -- '-' "${CPHMD_RT_TEMPLATE_MOCC_FILE}"; then
        cphmd_output_error "Active RT set is empty, but reduced-titration templates still contain active-site placeholders."
        return $?
    fi

    occ_tokens="$(awk '{print NF; exit}' "${CPHMD_RT_TEMPLATE_OCC_FILE}")"
    mocc_tokens="$(awk '{print NF; exit}' "${CPHMD_RT_TEMPLATE_MOCC_FILE}")"
    [[ "${occ_tokens}" == "${all_count}" && "${mocc_tokens}" == "${all_count}" ]] || {
        cphmd_output_error "Fixed RT templates do not match the ${all_count}-site all-site definition."
        return $?
    }

    cp -f -- "${CPHMD_RT_TEMPLATE_OCC_FILE}" "${CPHMD_RT_CURRENT_OCC_FILE}" || return 42
    cp -f -- "${CPHMD_RT_TEMPLATE_MOCC_FILE}" "${CPHMD_RT_CURRENT_MOCC_FILE}" || return 42
}

# cphmd_materialize_nonrt_cycle
#   Convert a validated all-site PETIT output into full OCC/MOCC lines when
#   reduced titration is disabled.
cphmd_materialize_nonrt_cycle ()
{
    local all_count
    local current_occ_tmp="${CPHMD_RT_CURRENT_OCC_FILE}.tmp.$$"
    local current_mocc_tmp="${CPHMD_RT_CURRENT_MOCC_FILE}.tmp.$$"
    local rc

    all_count="$(cphmd_count_nonblank "${CPHMD_ALL_SITES_FILE}")" || return 42
    cphmd_validate_petit_output TMP_MCarlo.out "${all_count}" || return $?

    gawk -v expected="${all_count}" -v outocc="${current_occ_tmp}" -v outmocc="${current_mocc_tmp}" '
        $1 == "f" {
            for (i=2; i<=NF; i++) printf "%d ", $i > outocc
            printf "\n" > outocc
            n = NF - 1
        }
        substr($0,1,1) == "." && $0 !~ /tot/ {mocc[$4+1] = $5}
        END {
            if (n != expected) exit 40
            for (i=1; i<=expected; i++) {
                if (!(i in mocc)) exit 41
                printf "%.10g ", mocc[i] > outmocc
            }
            printf "\n" > outmocc
            close(outocc); close(outmocc)
        }
    ' TMP_MCarlo.out
    rc=$?
    if (( rc != 0 )); then
        rm -f -- "${current_occ_tmp}" "${current_mocc_tmp}"
        cphmd_pbmc_error petit-output "Could not materialize full non-RT fractions (gawk exit ${rc})."
        return $?
    fi

    mv -f -- "${current_occ_tmp}" "${CPHMD_RT_CURRENT_OCC_FILE}" || return 42
    mv -f -- "${current_mocc_tmp}" "${CPHMD_RT_CURRENT_MOCC_FILE}" || return 42
}

# cphmd_record_current_fractions
#   Append the full applied state and population lines for the current cycle.
cphmd_record_current_fractions ()
{
    cphmd_require_nonempty_files "current protonation-state record is incomplete" \
        "${CPHMD_RT_CURRENT_OCC_FILE}" "${CPHMD_RT_CURRENT_MOCC_FILE}" || return $?
    cat -- "${CPHMD_RT_CURRENT_OCC_FILE}" >> TMP_CpHMD.occ || return 42
    cat -- "${CPHMD_RT_CURRENT_MOCC_FILE}" >> TMP_CpHMD.mocc || return 42
}

# cphmd_initialize_or_restore_rt_state
#   Remove the temporary legacy .sites name, initialize a fresh RT run, or
#   restore and validate active-site/templates from restart/ on continuation.
cphmd_initialize_or_restore_rt_state ()
{
    local all_count
    local expected_last_full
    local restart_dir="../${CPHMD_RESTART_DIR:-restart}"
    local restart_active="${restart_dir}/${blockname}.RT-active.sites"
    local restart_occ="${restart_dir}/${blockname}.RT-template.occ"
    local restart_mocc="${restart_dir}/${blockname}.RT-template.mocc"
    local active_count
    local dash_count
    local occ_tokens
    local mocc_tokens

    rm -f -- "${runname}.sites"
    all_count="$(cphmd_count_nonblank "${CPHMD_ALL_SITES_FILE}")" || {
        cphmd_output_error "All-site definition ${CPHMD_ALL_SITES_FILE} is missing."
        return $?
    }

    if (( ${CPHMD_DONE_CYCLES:-0} > 0 )); then
        if [[ "${CPHMD_REDUCED_TITRATION:-0}" != "${ReduceTitration:-0}" ]]; then
            cphmd_output_error \
                "Checkpoint reduced-titration mode (${CPHMD_REDUCED_TITRATION:-unset}) does not match settings (${ReduceTitration:-unset})."
            return $?
        fi

        if [[ -e "../${CPHMD_ALL_SITES_FILE}" ]] && ! cmp -s -- "${CPHMD_ALL_SITES_FILE}" "../${CPHMD_ALL_SITES_FILE}"; then
            cphmd_output_error "Regenerated all-site definition differs from the finalized checkpoint copy."
            return $?
        fi
    fi

    if (( ${ReduceTitration:-0} == 0 )); then
        rm -f -- "${CPHMD_RT_ACTIVE_SITES_FILE}" \
            "${CPHMD_RT_TEMPLATE_OCC_FILE}" "${CPHMD_RT_TEMPLATE_MOCC_FILE}" \
            "${CPHMD_RT_CURRENT_OCC_FILE}" "${CPHMD_RT_CURRENT_MOCC_FILE}"
        export CPHMD_RT_LAST_FULL_CYCLE=0
        export CPHMD_RT_ACTIVE_COUNT=0
        return 0
    fi

    if (( ${CPHMD_DONE_CYCLES:-0} == 0 )); then
        rm -f -- "${CPHMD_RT_ACTIVE_SITES_FILE}" \
            "${CPHMD_RT_TEMPLATE_OCC_FILE}" "${CPHMD_RT_TEMPLATE_MOCC_FILE}" \
            "${CPHMD_RT_CURRENT_OCC_FILE}" "${CPHMD_RT_CURRENT_MOCC_FILE}" \
            "${CPHMD_RT_HISTORY_FILE}"
        if (( all_count == 0 )); then
            : > "${CPHMD_RT_ACTIVE_SITES_FILE}"
            printf '\n' > "${CPHMD_RT_TEMPLATE_OCC_FILE}"
            printf '\n' > "${CPHMD_RT_TEMPLATE_MOCC_FILE}"
        fi
        export CPHMD_RT_LAST_FULL_CYCLE=0
        export CPHMD_RT_ACTIVE_COUNT=0
        return 0
    fi

    if [[ -f "../${SysName}_RT-sites.dat" ]]; then
        cp -f -- "../${SysName}_RT-sites.dat" "${CPHMD_RT_HISTORY_FILE}" || return 42
    fi
    if (( ${CPHMD_KEEP_RT_DEBUG:-0} == 1 )) && \
       [[ -f "../${CPHMD_LOG_DIR:-log-files}/${SysName}_RT-debug.pocc_RT" ]]; then
        cp -f -- "../${CPHMD_LOG_DIR:-log-files}/${SysName}_RT-debug.pocc_RT" \
            "${runname}.pocc_RT" || return 42
    fi

    if (( all_count == 0 )); then
        expected_last_full=0
    else
        expected_last_full=$(( 1 + ((CPHMD_DONE_CYCLES - 1) / RTInterval) * RTInterval ))
    fi
    if [[ "${CPHMD_RT_LAST_FULL_CYCLE:-0}" != "${expected_last_full}" ]]; then
        cphmd_output_error \
            "Checkpoint RT state reports last full cycle ${CPHMD_RT_LAST_FULL_CYCLE:-unset}; expected ${expected_last_full} after cycle ${CPHMD_DONE_CYCLES}."
        return $?
    fi

    [[ -e "${restart_active}" ]] || {
        cphmd_output_error "Missing checkpointed active RT site set ${restart_active}."
        return $?
    }
    cphmd_require_nonempty_files "missing checkpointed reduced-titration templates" \
        "${restart_occ}" "${restart_mocc}" || return $?

    cp -f -- "${restart_active}" "${CPHMD_RT_ACTIVE_SITES_FILE}" || return 42
    cp -f -- "${restart_occ}" "${CPHMD_RT_TEMPLATE_OCC_FILE}" || return 42
    cp -f -- "${restart_mocc}" "${CPHMD_RT_TEMPLATE_MOCC_FILE}" || return 42

    active_count="$(cphmd_count_nonblank "${CPHMD_RT_ACTIVE_SITES_FILE}")" || return 42
    occ_tokens="$(awk '{print NF; exit}' "${CPHMD_RT_TEMPLATE_OCC_FILE}")"
    mocc_tokens="$(awk '{print NF; exit}' "${CPHMD_RT_TEMPLATE_MOCC_FILE}")"
    dash_count="$(awk '{for(i=1;i<=NF;i++) if($i=="-") n++} END{print n+0}' "${CPHMD_RT_TEMPLATE_OCC_FILE}")"

    if [[ "${active_count}" != "${CPHMD_RT_ACTIVE_COUNT:-unset}" || \
          "${occ_tokens}" != "${all_count}" || "${mocc_tokens}" != "${all_count}" || \
          "${dash_count}" != "${active_count}" ]]; then
        cphmd_output_error \
            "Checkpointed reduced-titration state is inconsistent (all=${all_count}, active=${active_count}, checkpoint-active=${CPHMD_RT_ACTIVE_COUNT:-unset}, dashes=${dash_count})."
        return $?
    fi

    message W \
        "Restored reduced-titration state from cycle ${CPHMD_RT_LAST_FULL_CYCLE}: ${active_count}/${all_count} active sites."
}

# -----------------------------------------------------------------------------
# Optional settings-block preparation and validation
# -----------------------------------------------------------------------------
# The settings template may always contain #fixgro# and #plumed# example blocks.
# They are inert text until the corresponding feature is enabled. These helpers
# decide whether a block is needed and reject missing or untouched template
# placeholders before PB/MC begins.

# cphmd_fixgro_required
#   Return success when the historical fixgro PBC path is selected. The current
#   workflow uses fixgro only for non-membrane systems that either contain
#   multiple interacting solutes or include the solute through an external ITP.
#   Membrane systems instead use MembraneCenteringProtocol.
cphmd_fixgro_required ()
{
    (( ${memb:-0} == 0 \
       && ( ${multiple_solutes:-0} == 1 || ${include_itp:-0} == 1 ) ))
}

# cphmd_validate_fixgro_file FILE
#   Require a non-empty, user-adapted fixgro instruction file and the fixgro
#   executable. Comment-only/blank blocks and obvious template placeholders are
#   rejected. This function is called only when cphmd_fixgro_required succeeds.
cphmd_validate_fixgro_file ()
{
    local file="$1"
    local fixgro_bin="${CpHDIR}/scripts/fixgro"

    [[ -s "${file}" ]] || {
        message E "fixgro PBC treatment is enabled, but ${file} is missing or empty. Adapt the #fixgro# block in the settings file."
    }

    if ! grep -Eq '^[[:space:]]*[^#[:space:]]' "${file}"; then
        message E "fixgro PBC treatment is enabled, but ${file} contains only comments or blank lines."
    fi

    if grep -Eq 'TOTALATOM|REPLACE_[A-Z0-9_]+|/path/to/' "${file}"; then
        message E "The active #fixgro# block still contains template placeholders. Replace them with system-specific definitions before running."
    fi

    [[ -x "${fixgro_bin}" ]] || {
        message E "fixgro PBC treatment is enabled, but ${fixgro_bin} is missing or not executable."
    }
}

# cphmd_validate_plumed_file FILE
#   Require a non-empty PLUMED input when plumed=1. Obvious template paths and
#   replacement tokens are rejected so a copied example cannot silently reach
#   mdrun. PLUMED syntax itself is checked by the zero-step mdrun preflight.
cphmd_validate_plumed_file ()
{
    local file="$1"

    [[ -s "${file}" ]] || {
        message E "PLUMED is enabled, but ${file} is missing or empty. Adapt the #plumed# block in the settings file."
    }

    if ! grep -Eq '^[[:space:]]*[^#[:space:]]' "${file}"; then
        message E "PLUMED is enabled, but ${file} contains only comments or blank lines."
    fi

    if grep -Eq 'REPLACE_[A-Z0-9_]+|/path/to/' "${file}"; then
        message E "The active #plumed# block still contains template placeholders. Replace them with system-specific values before running."
    fi
}

# -----------------------------------------------------------------------------
# Runtime configuration and fail-fast GROMACS/GPU preflight
# -----------------------------------------------------------------------------
# These checks run before any PB/MC work. Their purpose is to stop immediately
# when the requested GROMACS mode cannot possibly run on the current host.
#
# Important distinction:
#   * `gmx --version` proves how GROMACS was compiled.
#   * `nvidia-smi` proves that an NVIDIA device and host driver are visible
#     inside the running container.
#
# A CUDA-enabled binary can print `GPU support: CUDA` even on a CPU-only host.
# Therefore, GPU=1 requires both checks unless CPHMD_GPU_PREFLIGHT=0 is set
# explicitly for troubleshooting.

# cphmd_require_boolean_setting NAME VALUE
#   Enforce the common 0/1 convention used by safety and output flags.
cphmd_require_boolean_setting ()
{
    local name="$1"
    local value="$2"

    if [[ "${value}" != "0" && "${value}" != "1" ]]; then
        message E "${name} must be 0 or 1; received '${value}'."
    fi
}

# cphmd_validate_runtime_settings
#   Validate GPU mode, safety flags, thread count and requested GPU ID before
#   constructing the mdrun command. Invalid settings stop the run immediately.
cphmd_validate_runtime_settings ()
{
    cphmd_require_boolean_setting "GPU" "${GPU:-}"
    cphmd_require_boolean_setting \
        "CPHMD_STRICT_GROMACS_MODE" \
        "${CPHMD_STRICT_GROMACS_MODE:-1}"
    cphmd_require_boolean_setting \
        "CPHMD_GPU_PREFLIGHT" \
        "${CPHMD_GPU_PREFLIGHT:-1}"
    cphmd_require_boolean_setting \
        "CPHMD_MDRUN_PREFLIGHT" \
        "${CPHMD_MDRUN_PREFLIGHT:-1}"

    if ! [[ "${nCPU:-}" =~ ^[1-9][0-9]*$ ]]; then
        message E "nCPU must be a positive integer; received '${nCPU:-unset}'."
    fi

    if [[ "${CPHMD_GPU_DEVICE_ID:-0}" != "auto" ]] \
       && ! [[ "${CPHMD_GPU_DEVICE_ID:-0}" =~ ^[0-9]+$ ]]
    then
        message E "CPHMD_GPU_DEVICE_ID must be 'auto' or a non-negative integer; received '${CPHMD_GPU_DEVICE_ID:-unset}'."
    fi
}

# cphmd_build_mdrun_command
#   Build CPHMD_MDRUN_CMD as a Bash array. Empty mdruncpu/mdrungpu values select
#   safe defaults based on GPU, nCPU and CPHMD_GPU_DEVICE_ID; non-empty values
#   are treated as expert overrides without using `eval`. Also records a quoted
#   command string for readable logging.
cphmd_build_mdrun_command ()
{
    local custom_command=""
    local -a command_tail=()
    local gpu_id="${CPHMD_GPU_DEVICE_ID:-0}"
    local i

    cphmd_validate_runtime_settings

    export OMP_NUM_THREADS="${nCPU}"

    if (( GPU == 1 )); then
        custom_command="${mdrungpu:-}"

        if [[ -n "${custom_command}" ]]; then
            # Advanced override. This is split into argv words without `eval`;
            # shell pipelines, redirects, command substitutions and embedded
            # quoting are intentionally unsupported.
            read -r -a command_tail <<< "${custom_command}"
        else
            command_tail=(
                mdrun
                -ntmpi 1
                -ntomp "${nCPU}"
                -pin on
                -nb gpu
                -pme auto
                -bonded auto
                -update auto
            )

            if [[ "${gpu_id}" != "auto" ]]; then
                command_tail+=( -gpu_id "${gpu_id}" )
            fi
        fi

        # Determine which device ID the actual command will request. A custom
        # mdrungpu value containing -gpu_id takes precedence over the general
        # CPHMD_GPU_DEVICE_ID setting.
        export CPHMD_EFFECTIVE_GPU_DEVICE_ID="${gpu_id}"
        for (( i=0; i<${#command_tail[@]}; i++ )); do
            if [[ "${command_tail[i]}" == "-gpu_id" ]]; then
                if (( i + 1 >= ${#command_tail[@]} )); then
                    message E "mdrungpu contains -gpu_id without a following device ID."
                fi
                export CPHMD_EFFECTIVE_GPU_DEVICE_ID="${command_tail[i+1]}"
                break
            fi
        done

        if [[ "${CPHMD_EFFECTIVE_GPU_DEVICE_ID}" != "auto" ]] \
           && ! [[ "${CPHMD_EFFECTIVE_GPU_DEVICE_ID}" =~ ^[0-9]+$ ]]
        then
            message E "The effective GPU device ID must be 'auto' or one non-negative integer. Received '${CPHMD_EFFECTIVE_GPU_DEVICE_ID}'."
        fi

        message W "GPU execution requested."
    else
        custom_command="${mdruncpu:-}"

        if [[ -n "${custom_command}" ]]; then
            read -r -a command_tail <<< "${custom_command}"
        else
            command_tail=(
                mdrun
                -nt "${nCPU}"
                -pin on
            )
        fi

        unset CPHMD_EFFECTIVE_GPU_DEVICE_ID
        message W "CPU execution requested."
    fi

    if (( ${#command_tail[@]} == 0 )) || [[ "${command_tail[0]}" != "mdrun" ]]; then
        message E "The selected mdrun command must begin with 'mdrun'. Received: ${custom_command:-<automatic command>}"
    fi

    CPHMD_MDRUN_CMD=( "${GroDIR}" "${command_tail[@]}" )

    printf -v CPHMD_MDRUN_DISPLAY '%q ' "${CPHMD_MDRUN_CMD[@]}"
    CPHMD_MDRUN_DISPLAY="${CPHMD_MDRUN_DISPLAY% }"
    export CPHMD_MDRUN_DISPLAY

    message W "mdrun command: ${CPHMD_MDRUN_DISPLAY}"
}

# cphmd_preflight_gromacs_runtime
#   Run `gmx --version`, verify CPU/CUDA compile mode, and—when GPU=1—confirm a
#   visible NVIDIA device with nvidia-smi. This check happens before PB/MC and
#   catches a wrong executable, missing --nv, hidden GPU or invalid GPU ID.
cphmd_preflight_gromacs_runtime ()
{
    local version_file
    local version=""
    local gpu_support=""
    local gpu_listing=""
    local gpu_count=0
    local requested_gpu="${CPHMD_EFFECTIVE_GPU_DEVICE_ID:-${CPHMD_GPU_DEVICE_ID:-0}}"
    local gpu_summary=""

    if [[ ! -x "${GroDIR}" ]]; then
        message E "Selected GROMACS executable is missing or not executable: ${GroDIR}"
    fi

    version_file="$(mktemp "${TMPDIR:-/tmp}/cphmd-gromacs-version.XXXXXX")" \
        || message E "Could not create a temporary GROMACS-version file."

    if ! "${GroDIR}" --version > "${version_file}" 2>&1; then
        cat "${version_file}" >&2 || true
        rm -f -- "${version_file}"
        message E "Unable to execute '${GroDIR} --version'."
    fi

    version="$(
        awk -F: '
            /^GROMACS version:/ {
                v=$2
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
                print v
                exit
            }
        ' "${version_file}"
    )"

    gpu_support="$(
        awk -F: '
            /^GPU support:/ {
                v=$2
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
                print v
                exit
            }
        ' "${version_file}"
    )"

    rm -f -- "${version_file}"

    if [[ -z "${version}" ]]; then
        message E "Could not parse the GROMACS version from '${GroDIR} --version'."
    fi

    if [[ -z "${gpu_support}" ]]; then
        message E "Could not parse GPU support from '${GroDIR} --version'."
    fi

    export CPHMD_GROMACS_VERSION="${version}"
    export CPHMD_GROMACS_GPU_SUPPORT="${gpu_support}"

    if (( GPU == 1 )); then
        if [[ "${gpu_support}" != "CUDA" ]]; then
            message E "GPU=1 was requested, but ${GroDIR} reports 'GPU support: ${gpu_support}', not CUDA. Select the CUDA GROMACS build or set GPU=0."
        fi
    elif (( ${CPHMD_STRICT_GROMACS_MODE:-1} == 1 )); then
        if [[ "${gpu_support}" != "disabled" ]]; then
            message E "GPU=0 requested the dedicated CPU build, but ${GroDIR} reports 'GPU support: ${gpu_support}'. Correct GroDIRCPU or set CPHMD_STRICT_GROMACS_MODE=0 to permit CPU execution with a GPU-enabled build."
        fi
    fi

    message W "Using GROMACS ${version}; GPU support: ${gpu_support}."

    if (( GPU == 0 )); then
        return 0
    fi

    if (( ${CPHMD_GPU_PREFLIGHT:-1} == 0 )); then
        message W "WARNING: CPHMD_GPU_PREFLIGHT=0. GPU device visibility will not be checked before PB/MC."
        return 0
    fi

    if [[ ${CUDA_VISIBLE_DEVICES+x} == x ]]; then
        case "${CUDA_VISIBLE_DEVICES}" in
            ""|-1|NoDevFiles|none|None)
                message E "GPU=1 was requested, but CUDA_VISIBLE_DEVICES='${CUDA_VISIBLE_DEVICES}' disables GPU visibility."
                ;;
            *)
                message W "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
                ;;
        esac
    else
        message W "CUDA_VISIBLE_DEVICES is unset; checking devices directly with nvidia-smi."
    fi

    if ! command -v nvidia-smi >/dev/null 2>&1; then
        message E "GPU=1 was requested, but nvidia-smi is unavailable inside the container. On Apptainer/Singularity, launch the container with --nv on an NVIDIA GPU node."
    fi

    if ! gpu_listing="$(nvidia-smi -L 2>&1)"; then
        printf '%s\n' "${gpu_listing}" >&2
        message E "GPU=1 was requested, but nvidia-smi could not access an NVIDIA device/driver. Refusing to enter PB/MC."
    fi

    gpu_count="$(grep -c '^GPU [0-9][0-9]*:' <<< "${gpu_listing}" || true)"
    if (( gpu_count < 1 )); then
        printf '%s\n' "${gpu_listing}" >&2
        message E "GPU=1 was requested, but no NVIDIA GPU was visible. Refusing to enter PB/MC."
    fi

    if [[ "${requested_gpu}" != "auto" ]]; then
        if ! [[ "${requested_gpu}" =~ ^[0-9]+$ ]]; then
            message E "Requested GPU ID '${requested_gpu}' is not a single non-negative integer."
        fi
        if (( requested_gpu >= gpu_count )); then
            printf '%s\n' "${gpu_listing}" >&2
            message E "Requested GPU ID ${requested_gpu}, but only ${gpu_count} visible NVIDIA GPU(s) were detected."
        fi
    fi

    gpu_summary="$(
        nvidia-smi \
            --query-gpu=index,name,driver_version \
            --format=csv,noheader 2>/dev/null \
        || true
    )"

    if [[ -n "${gpu_summary}" ]]; then
        message W "Visible NVIDIA GPU(s): ${gpu_summary//$'\n'/; }"
    else
        message W "Visible NVIDIA GPU count: ${gpu_count}"
    fi

    if [[ "${requested_gpu}" == "auto" ]]; then
        message W "GROMACS will select from the visible GPU devices automatically."
    else
        message W "Preflight accepted requested GPU device ID ${requested_gpu}."
    fi
}


# cphmd_preflight_mdrun_execution
#   Execute the real configured mdrun command against TMP_CpHMD.tpr with
#   `-nsteps 0`. This initializes the selected CPU/GPU path without advancing
#   the trajectory, detecting runtime/driver/option errors before PB/MC starts.
cphmd_preflight_mdrun_execution ()
{
    # The preflight runs in a new private directory once per scheduler
    # invocation. Because that directory is empty, normal append mode is safe
    # and avoids GROMACS' numbered `.partNNNN` output names.
    local preflight_dir=".cphmd-mdrun-preflight.$$"
    local preflight_output="${preflight_dir}/stdout-stderr.log"
    local preflight_log="${preflight_dir}/preflight.log"
    local tpr_path=""
    local rc=0
    local log_dir="../${CPHMD_LOG_DIR:-log-files}"
    local failure_log="${log_dir}/${blockname}_mdrun-preflight_job${SLURM_JOB_ID:-nojob}_pid$$.log"
    local failure_tmp="${failure_log}.tmp.$$"

    if (( ${CPHMD_MDRUN_PREFLIGHT:-1} == 0 )); then
        message W "WARNING: CPHMD_MDRUN_PREFLIGHT=0. The real mdrun command will not be tested before PB/MC."
        return 0
    fi

    cphmd_require_nonempty_files \
        "mdrun preflight requires a valid TMP_CpHMD.tpr" \
        TMP_CpHMD.tpr || return $?

    tpr_path="$(readlink -f TMP_CpHMD.tpr)" \
        || message E "Could not resolve the preflight TPR path."

    rm -rf -- "${preflight_dir}"
    mkdir -p -- "${preflight_dir}" \
        || message E "Could not create ${preflight_dir}."

    (
        cd "${preflight_dir}" || exit 125

        # State append mode explicitly after any expert override. This prevents
        # a custom command from generating numbered part files in the disposable
        # preflight directory. The log filename is also fixed explicitly.
        "${CPHMD_MDRUN_CMD[@]}" \
            -s "${tpr_path}" \
            -nsteps 0 \
            -deffnm preflight \
            -g preflight.log \
            -append
    ) > "${preflight_output}" 2>&1
    rc=$?

    if (( rc != 0 )) || [[ ! -s "${preflight_log}" ]]; then
        mkdir -p -- "${log_dir}" 2>/dev/null || true

        {
            echo "=== GROMACS mdrun preflight failure ==="
            echo "TIMESTAMP=$(date -Is 2>/dev/null || date)"
            echo "JOB_ID=${SLURM_JOB_ID:-unset}"
            echo "JOB_NAME=${SLURM_JOB_NAME:-unset}"
            echo "HOST=$(hostname 2>/dev/null || echo unknown)"
            echo "EXIT_CODE=${rc}"
            echo "EXPECTED_LOG=${preflight_log}"
            echo "LOG_PRESENT=$([[ -s "${preflight_log}" ]] && echo 1 || echo 0)"
            echo "COMMAND=${CPHMD_MDRUN_DISPLAY} -s ${tpr_path} -nsteps 0 -deffnm preflight -g preflight.log -append"
            echo

            echo "=== stdout and stderr ==="
            cat -- "${preflight_output}" 2>/dev/null || true
            echo

            if [[ -s "${preflight_log}" ]]; then
                echo "=== GROMACS log ==="
                cat -- "${preflight_log}" 2>/dev/null || true
                echo
            fi

            echo "=== Preflight-directory inventory ==="
            find "${preflight_dir}" -maxdepth 1 -type f \
                -printf '%f\t%s bytes\n' 2>/dev/null \
              | sort || true
        } > "${failure_tmp}" 2>/dev/null || true

        if [[ -s "${failure_tmp}" ]]; then
            mv -f -- "${failure_tmp}" "${failure_log}" 2>/dev/null || true
        else
            rm -f -- "${failure_tmp}" 2>/dev/null || true
        fi

        cat -- "${preflight_output}" >&2 2>/dev/null || true

        if (( rc != 0 )); then
            message E "The configured GROMACS mdrun command failed its zero-step preflight (exit code ${rc}). PB/MC was not started. Diagnostic: ${failure_log}"
        fi

        message E "The zero-step mdrun returned success but produced no readable GROMACS log at ${preflight_log}. PB/MC was not started. Diagnostic: ${failure_log}"
    fi

    rm -rf -- "${preflight_dir}"
    message W "Zero-step mdrun preflight passed before PB/MC."
}

# correct_variables
#   Normalize selected numeric settings after the settings file has been
#   sourced. Historically this prevents empty/zero values from being handled
#   inconsistently by later AWK expressions.
correct_variables ()
{                   
    #removed offset value since CpHMD sources this value
    GridSize=`awk -v i=$GridSize 'BEGIN{print (i==0) ? 0 : i}'`
    RTThreshold=`awk -v i=$RTThreshold 'BEGIN{print (i==0) ? 0 : i}'`
}

# build_forcefield
#   Copy the selected constant-pH force-field directory into the isolated
#   working directory so GROMACS can resolve its local force-field files.
build_forcefield ()
{
    if [ ! -d $ffDIR ]; then 
          message E  "Force field directory not defined in the .settings!!!... Program will crash"
    fi
    
    # Call the modified forcefield:
    cp -rdf "$ffDIR"/${ffID}.ff .;     
}

# check_files
#   Perform the main compatibility audit before the CpHMD loop: topology naming,
#   required programs/directories, PB dimensionality, GROMACS/cutoff support,
#   GPU-incompatible energy groups, force-field cutoffs and key MDP constraints.
#   This is validation only; it must not advance the simulation.
check_files ()
{

    # Check topology to see if the name of the molecule is set to protein, otherwise it will crash
    mol_name=`awk '/nrexcl/{getline; print $1}' $TOPin`
    if [[ ! $mol_name =~ "Protein" ]] ; then
	message E  "Topology does not have the molecule set as Protein. Program will crash"
    fi

    # Check directories
    for d in $PetitDIR $StDIR $DelphiDir
      do
      if [ ! -d $d ]; then 
          message E  "Directory $d does not exist!!!... Program will crash"
      fi
    done    

    if [[ $PBdim -ne 0 && $PBdim -ne 2 ]]; then
        message E  "Number of dimensions ($PBdim) is not supported in ${CpHModule} module!!!... Program will crash"
    fi
    #
    # GROMACS compile mode and GPU visibility were checked before the scratch
    # directory was created. Reuse those recorded values here for MDP-specific
    # compatibility checks.
    version="${CPHMD_GROMACS_VERSION:-}"
    gpu_support="${CPHMD_GROMACS_GPU_SUPPORT:-}"

    if [[ -z "${version}" || -z "${gpu_support}" ]]; then
        message E "Internal error: GROMACS preflight results are unavailable."
    fi

    # Detect cutoff-scheme to be used.
    scheme=`awk '/cutoff-scheme/ {print $3}' ${runname}.mdp `
    if [[ $scheme == "verlet" ]]; then
        message W "Scheme used is verlet."
    elif [[ $scheme == "group" ]]; then
        message W "Scheme used is Group based."
        if [[ $version != "5.1.5" ]]; then
            message E  "Cut-off scheme ${scheme} is not supported in ${GroDIR} module!!!... Program will crash"
        fi
    else
        message W "Scheme is not group nor verlet... Performing method as group based cut off"
    fi

    ###########################################################################
    ## Check if gromacs version is not 2020.2 since it has freezedims bugged ##
    ###########################################################################
    if [[ $version == "2020.2" ]]; then
        message E  "Gromacs ${version} is not supported for CpHMD since freezedims is bugged !!... Program will crash"
    fi

    #############################################
    ### Confirm the preflight-selected GROMACS mode ###
    #############################################
    # This is a defense-in-depth check. The primary executable/device preflight
    # has already run before any scratch setup or PB/MC preparation.
    if [[ ${GPU} -eq 1 && "${gpu_support}" != "CUDA" ]]; then
        message E "Internal error: GPU mode reached check_files without a CUDA GROMACS build."
    fi

    if [[ ${GPU} -eq 0 \
          && ${CPHMD_STRICT_GROMACS_MODE:-1} -eq 1 \
          && "${gpu_support}" != "disabled" ]]
    then
        message E "Internal error: strict CPU mode reached check_files with GPU support '${gpu_support}'."
    fi

    ##evaluate if mdp contains only one energy group ##
    if [[ `awk -F "=" '/energygrps/ {print $2}' ${runname}.mdp | awk '{print NF}'` > 1 && `awk -F "=" '/energygrps/ {print $2}' ${runname}.mdp ` != "System" && $GPU == 1 ]] ; then
	message E  "ERROR: Use of GPU and multiple energy groups is not supported, simulations would be downgraded to CPU only. Please either use only 1 energy groups or remove GPU support."
    fi

    ############################
    ## Check on .mdp settings ##
    ############################

    if (( ${mdpoverride:-0} == 0 )) ; # if override is turned off check values
    then
	# Read scalar MDP values robustly. This accepts common GROMACS formatting:
	#   key = value
	#   key    =    value ; comment
	# and compares numerically rather than as strings, so 1 and 1.0 are equal.
	mdp_get_number () {
	    local key="$1"
	    local file="$2"

	    awk -F '=' -v key="$key" '
		{
		    lhs=$1
		    gsub(/^[[:space:]]+|[[:space:]]+$/, "", lhs)
		}
		lhs == key {
		    rhs=$2
		    sub(/;.*/, "", rhs)
		    gsub(/^[[:space:]]+|[[:space:]]+$/, "", rhs)
		    print rhs + 0
		    found=1
		    exit
		}
		END {
		    if (!found) exit 1
		}
	    ' "$file"
	}

	mdp_number_eq () {
	    local observed="$1"
	    local expected="$2"

	    awk -v observed="$observed" -v expected="$expected" '
		BEGIN {
		    diff = observed - expected
		    if (diff < 0) diff = -diff
		    exit (diff <= 1e-6 ? 0 : 1)
		}
	    '
	}

	rvdw_value="$(mdp_get_number "rvdw" "${runname}.mdp" || true)"
	rcoulomb_value="$(mdp_get_number "rcoulomb" "${runname}.mdp" || true)"

	if [[ -z "${rvdw_value}" ]]; then
	    message E "Error: rvdw was not found in ${runname}.mdp. Please define rvdw or set mdpoverride=1 if this is intentional."
	fi

	if [[ -z "${rcoulomb_value}" ]]; then
	    message E "Error: rcoulomb was not found in ${runname}.mdp. Please define rcoulomb or set mdpoverride=1 if this is intentional."
	fi

	## check for force-switch settings with GROMOS
	if grep -Eq '^[[:space:]]*vdw-modifier[[:space:]]*=' "${runname}.mdp" && [[ "${ffID}" == "G54a7pH" ]] ; then
	    message E "Error: vdw-modifier was given in the MDP for the GROMOS force field. Confirm this option is intended. If it is intended, add to your settings file the flag: export mdpoverride=1 "
	fi

	## Check if ff gromos vdw and coul are compatible with the expected defaults
	if [[ "${ffID}" == "G54a7pH" ]] ; then
	    if { ! mdp_number_eq "${rvdw_value}" "1.0" && ! mdp_number_eq "${rvdw_value}" "1.4"; } || \
	       ! mdp_number_eq "${rcoulomb_value}" "${rvdw_value}" ; then
		message E "Error: GROMOS force field selected but rvdw/rcoulomb are not compatible with the default 1.4 protein or 1.0 membrane settings. Detected rvdw=${rvdw_value}, rcoulomb=${rcoulomb_value}. If this is intended, add to your settings file the flag: export mdpoverride=1 "
	    fi
	fi

	## Check if ff charmm vdw and coul = 1.2
	if [[ "${ffID}" == "CHARMM36pH" ]] ; then
	    if ! mdp_number_eq "${rvdw_value}" "1.2" || ! mdp_number_eq "${rcoulomb_value}" "1.2" ; then
		message E "Error: CHARMM force field selected but rvdw and rcoulomb are not the default 1.2. Detected rvdw=${rvdw_value}, rcoulomb=${rcoulomb_value}. If this is intended, add to your settings file the flag: export mdpoverride=1 "
	    fi
	fi

	## Check if ff Amber vdw and coul = 1.0
	if [[ "${ffID}" == "Amber14SBpH" ]] ; then
	    if ! mdp_number_eq "${rvdw_value}" "1.0" || ! mdp_number_eq "${rcoulomb_value}" "1.0" ; then
		message E "Error: AMBER force field selected but rvdw and rcoulomb are not the default 1.0. Detected rvdw=${rvdw_value}, rcoulomb=${rcoulomb_value}. If this is intended, add to your settings file the flag: export mdpoverride=1 "
	    fi
	fi
    fi

    #################################
    ## Check RT values if RT is on ##
    #################################
    if [[ $ReduceTitration == 1 ]] ; then
	if [ -z $RTInterval ] ; then
	    message E "Error: Reduced Titration has been turned on but no interval was given."
	else
	    message W "Reduced titration has been selected to run at intervals of $RTInterval cycles with a threshold of $RTThreshold "
	fi
    fi
    
    ######################################################
    ## Validate optional fixgro and PLUMED inputs in use ##
    ######################################################
    if cphmd_fixgro_required; then
        cphmd_validate_fixgro_file "${runname}.fixgro"
    fi

    if [[ ${plumed:-0} != "0" ]]; then
        cphmd_validate_plumed_file "${runname}_plumed.dat"
    fi

    #############################################################
    ## Optional ASIC/generalized behaviour controls            ##
    #############################################################
    # MembraneCenteringProtocol="tails" keeps the historical membrane
    # centering workflow and requires the Onetail/Monotail/Bitail index
    # groups. MembraneCenteringProtocol="central_atoms" uses the ASIC
    # workflow and requires CentralAtom/CentralAtoms groups.
    : "${MembraneCenteringProtocol:=tails}"
    if [[ $MembraneCenteringProtocol != "tails" && $MembraneCenteringProtocol != "central_atoms" ]]; then
	message E "MembraneCenteringProtocol=$MembraneCenteringProtocol is not supported. Use tails or central_atoms."
    fi

    # UseInputOffset=0 preserves the historical automatic offset from the
    # last protein residue number. UseInputOffset=1 keeps the user-provided
    # offset value from the settings file.
    : "${UseInputOffset:=0}"
    if [[ $UseInputOffset != 0 && $UseInputOffset != 1 ]]; then
	message E "UseInputOffset=$UseInputOffset is not supported. Use 0 or 1."
    fi

    : "${PBMCDebugPDB:=0}"
    if [[ $PBMCDebugPDB != 0 && $PBMCDebugPDB != 1 ]]; then
	message E "PBMCDebugPDB=$PBMCDebugPDB is not supported. Use 0 or 1."
    fi
	
}


# make_auxiliary_files
#   Build the initial temporary CpHMD working set (TMP_CpHMD.*): processed
#   coordinates, topology, index, TPR, empty accumulation files and auxiliary
#   coordinate representations required by PB/MC. Existing stale temporary
#   outputs are removed before generation.
make_auxiliary_files ()
{
    # Rename original files
    cp -f $TOPin TMP_CpHMD.top
    cp -f $NDXin TMP_CpHMD.ndx
    cp -f $GROin TMP_effective.gro
    
  
    # Check if there are different ionic strength values in MD and PB
    if [[ $ionicstrMD == "" ]]; then ionicstrMD=$ionicstr; fi
    # Convert Ionic Strength from Molar to molecule/nm^3.
    ionicstrMolecule=$(awk -v i=${ionicstrMD} 'BEGIN{print i*0.6022}')
    #
    # Correct Ionic Strength and Temperature in the .mdp file according 
    # to your parameters
    sed "s/\(ionicstrength *= \).*\$/\1 $ionicstrMolecule/" \
        ${runname}.mdp > TMP_aux1.mdp
    message W "Ionic strength in ${runname}.mdp file was changed to $ionicstrMolecule molecules/nm^3 (which corresponds to $ionicstr Molar)."
    #
    tcgroups=`awk '$1=="tc-grps"{for (i=3;i<=NF;i++) if ($i !~ /^;/) {n++; if ($i ~ /;/) exit} else exit}END{print n}' ${runname}.mdp` #SC 25-11-2011
    awk -v e=$tcgroups -v t=$temp '!/^ref_t *=/{print $0}; 
        /^ref_t *=/{printf "ref_t               =  "; 
        for (i=1;i<=e;i++) printf "%-9s", t ;print ""}' \
            TMP_aux1.mdp > ${runname}.mdp #SC 25-11-2011
        message W "Temperature in ${runname}.mdp file was changed to $temp."

    # Make relaxation .mdp gromacs file
    # Remove constraints & change NPT to NVT
    sed "s/\(nsteps *= \).*\$/\1 $RelaxSteps/ 
         s/\(Pcoupl *= \).*\$/\1 No/" \
        ${runname}.mdp > TMP_relax.mdp
    #
    #### Adition on 13/01/2021 to run in either group and verlet cutoff based runs ####
    ### Add Different protocol depending on the cut-off scheme ###
    # Check which cut-off is being used #
    scheme=`awk '/cutoff-scheme/ {print $3}' ${runname}.mdp `
    if [[ $scheme == "verlet" ]]
    then
        echo -e "\nfreezegrps          =  Protein\nfreezedim           =\
        Y Y Y\n" >> TMP_relax.mdp
    elif [[ $scheme == "group" ]]
    then
        echo -e "\nfreezegrps          =  Protein\nfreezedim           =\
        Y Y Y\n\nenergygrp_excl      =  Protein Protein" >> TMP_relax.mdp
    else
        message W "Scheme not group nor verlet... Performing method as group based cut off"
        echo -e "\nfreezegrps          =  Protein\nfreezedim           =\
        Y Y Y\n\nenergygrp_excl      =  Protein Protein" >> TMP_relax.mdp

    fi
    #
    # Make effective .mdp gromacs file
    sed  "s/\(nsteps *= \).*\$/\1 $EffectiveSteps/" ${runname}.mdp > TMP_effective.mdp
    #
    # Modify Effective MDP file in case of using Position Restraints
    if [[ -f $PosRe && $PosRe != *(" ") ]]; then 
        sed -i "s/\(define *= \).*\$/\1 -DPOSRES/" TMP_effective.mdp
    else
        if [ ! -f $PosRe ]; then message W \
        "File $PosRe is missing. Position Restraints will not be used"; fi
    fi
    #
    ### 24/04 changed -p $TOPin to -p TMP_CpHMD.top for the grompp
    ##to run using the folder's FF
    if ! "$GroDIR" grompp \
        -f TMP_effective.mdp \
        -po TMP_effective_out.mdp \
        -c "$GROin" \
        -p TMP_CpHMD.top \
        -pp TMP_processed.top \
        -n TMP_CpHMD.ndx \
        -o TMP_CpHMD.tpr \
        -maxwarn 1000 \
        -quiet
    then
        message E "Initial grompp failed before PB/MC."
    fi

    cphmd_require_nonempty_files \
        "initial grompp did not create TMP_CpHMD.tpr" \
        TMP_CpHMD.tpr || return $?

    rm -f TMP_aux*.mdp

    #################################################################
    ### Getting the correct pbp file for PBMC                     ###
    #################################################################
    if [[ $PBdim == 0 ]]
    then
	sed "s/\(export ionicstr*=\).*\$/\1$ionicstr/
             s/\(export epsin*=\).*\$/\1$epsin/" \
		 "$DelphiDir"/DELPHI_Prot.pbp > DELPHI.pbp
    elif [[ $PBdim == 2 ]]
    then
	 sed "s/\(export ionicstr*=\).*\$/\1$ionicstr/
             s/\(export epsin*=\).*\$/\1$epsin/" \
		 "$DelphiDir"/DELPHI_Memb.pbp > DELPHI.pbp
    fi
    
}

# make_sites SETTINGS_FILE
#   Identify titratable residues and generate the .sites definitions used by
#   Delphi/PETIT. It also prepares reduced-titration site files when that mode
#   is enabled. These site definitions persist across scheduler resubmissions.
make_sites ()
{
    echo "Protein" | "$GroDIR" editconf -f $GROin \
				   -o TMP_protein.pdb -n TMP_CpHMD.ndx -quiet
    
    
    "$CpHDIR"/scripts/make_sites $1 TMP_protein.pdb 
    #
    # Correct TMP_CpHMD.top with NTR and CTR
    #
    #### Previously ######
    ## Would only correct it if NT was in sites, now it will try to do it
    ## Always if it has the H3 ATOM

    if [[ $ffID ==  G54a7pH ]];then 
        awk '
            BEGIN{
		top=ARGV[1]; 
            	while (getline < top)
            	if ($5=="H3") ntr[++j]=$3;
            	close(top)}; 
        	END{
            	while (getline < top) {
            	write=1; 
              	for (i in ntr) {
                    if ($3==ntr[i] && $5=="N") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="H1") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="H2") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="H3") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="CA") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="C") 
                       {print substr($0,1,27),"NTR", substr($0,33); write=0};
                    if ($3==ntr[i] && $5=="O") 
                      {print substr($0,1,27),"NTR", substr($0,33); write=0};
              	};
              	if (write) print $0}
        }' TMP_CpHMD.top > TMP_aux.top

    elif [[ $ffID == Amber14SBpH ]] ;then 
	awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="H3") ntr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ntr) {
                if ($3==ntr[i] && $5=="N") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H1") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H2") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H3") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="CA") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
		if ($3==ntr[i] && $5=="HA") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="C") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="O") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_CpHMD.top > TMP_aux.top

    elif [[ $ffID == CHARMM36pH ]] ;then 
	awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="H3") ntr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ntr) {
                if ($3==ntr[i] && $5=="N") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H1") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H2") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="H3") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
                if ($3==ntr[i] && $5=="CA") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
		if ($3==ntr[i] && $5=="HA") 
                   {print substr($0,1,27),"NTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_CpHMD.top > TMP_aux.top
	
    fi
	
    #if egrep CT ${runname}.sites; then
    if [[ $ffID ==  G54a7pH ]];then 
        awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="O1") ctr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ctr) {
                if ($3==ctr[i] && $5=="C") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="O1") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="O2") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO11") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO12") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO21") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HO22") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_aux.top > TMP_CpHMD_charge.top
	
    elif [[ $ffID == Amber14SBpH ]] ;then #added AMBER case with CTR names
	awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="OC1") ctr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ctr) {
	      	if ($3==ctr[i] && $5=="N") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="H") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="CA") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HA") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="C") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="OC1") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="OC2") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC11") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC12") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC21") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HC22") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_aux.top > TMP_CpHMD_charge.top
    elif [[ $ffID == CHARMM36pH ]] ;then #added AMBER case with CTR names
	awk '
        BEGIN{
            top=ARGV[1]; 
            while (getline < top)
            if ($5=="OT1") ctr[++j]=$3;
            close(top)}; 
        END{
            while (getline < top) {
            write=1; 
              for (i in ctr) {
                if ($3==ctr[i] && $5=="C") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="OT1") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="O2") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT11") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT12") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT21") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
                if ($3==ctr[i] && $5=="HT22") 
                   {print substr($0,1,27),"CTR", substr($0,33); write=0};
              };
              if (write) print $0}
        }' TMP_aux.top > TMP_CpHMD_charge.top
    else
	### May become bugged when unknown ffs are given but
	### Solution added for testing the XOL3 ff from tomás
	### It is always expecting a TMP_CpHMD_charge.top
	### hence it was needed for creation
	cat TMP_CpHMD.top > TMP_CpHMD_charge.top
    fi

    ##################### Adition for offset check ############################
    ## Define Offset and verify if it is able the system is able to compute  ##
    ###########################################################################
    ## UseInputOffset=0 keeps the historical automatic offset.               ##
    ## UseInputOffset=1 keeps the value supplied in the settings file.       ##
    ###########################################################################
    
    if [[ ${UseInputOffset:-0} -eq 1 ]]; then
	export offset=${offset:-2000}
	message W "Using user-defined PB/MC offset from settings: $offset"
    else
	last_res=`awk '/ATOM/ {print substr($0,23,4)}' TMP_protein.pdb | tail -n 1`
	export offset=`echo $last_res | awk '{print $1 +10}'`
	message W "Using automatically determined PB/MC offset: $offset"
    fi

    last_tit=`awk 'NF{value=$1} END{print value+0}' ${runname}.sites`

    #############################
    if [ $(( (offset * 2) +500 + last_tit )) -ge 10000 ] ; then
	message E "System simulated will wield a PB/MC renumbering larger than 10 000 which will crash the PBMC cycle. Your system is currently unfeasible to calculate with CpHMD."
    fi

    # make_sites is a legacy external helper and writes ${runname}.sites. From
    # this point onward CpHMD uses explicit names: -all.sites for the immutable
    # definition and -RT-active.sites for the current reduced subset.
    mv -f -- ${runname}.sites "${CPHMD_ALL_SITES_FILE}"
    
}

# make_delphi_DB
#   Select and prepare the Delphi charge/radius databases for the chosen
#   force field and current titration state. The generated database inputs are
#   consumed by run_PBMC during each constant-pH cycle.
make_delphi_DB ()
{
    echo "atom__res_radius_" > ./DataBaseT.siz
    echo "atom__resnumbc_charge_" > ./DataBaseT.crg
    for type in crg siz ; do
	#####################################################################
	## Resolve and validate the Delphi database used for this run.     ##
	## This keeps the historical DatabaseDIR interface, but makes       ##
	## container/bind-mount mistakes explicit in the error messages.    ##
	#####################################################################
	db_file="$DatabaseDIR/DataBaseT_${ffID}.${type}"

	if [[ -z "$DatabaseDIR" ]]; then
            message E "DatabaseDIR is not defined. Set export DatabaseDIR to the directory containing DataBaseT_${ffID}.${type}."
	fi

	if [[ ! -f "$db_file" ]]; then
            message E "Delphi database file not found: $db_file"
	fi

	message W "Using Delphi database file: $db_file"

	## Insert the terminals in the database
	for ter in CT NT; do
	    case $type in
		crg)   awk -v t=$ter '$2~t {printf"%-6s%-9s%6.3f\n", $1,$2,$3}' "$db_file" >> ./DataBaseT.${type} ;;
		siz)   awk -v t=$ter '$2~t {printf"%-6s%-6s%-6.3f\n", $1,$2,$3}' "$db_file" >> ./DataBaseT.${type} ;;
	    esac
	done
	## Insert other residues
	for db_res in `awk '/ATOM/ {print substr($0,18,4)}' TMP_protein.pdb | sort | uniq` ; do	    
	    if [[ ! -z `awk -v d=$db_res '$2==d {print}' "$db_file"` ]] ; then
		## If the third character of the residue name is a number (hence a CpHMD residue most likely) ##
		if [[ `echo $db_res | awk '{if (substr($1,3,1) ~ /^[0-9]/) {print 1}}' ` == 1 ]] ;then
		    ## Check for the residue in question to not be in the database already ##
		    if [[ -z `awk  -v d=${db_res} '$2==d {print}' ./DataBaseT.${type}` ]] ; then
			case $type in
			    crg)   awk -v d=${db_res} '$2~substr(d,0,2) && substr($2,3,1)~/[0-9]/ {printf"%-6s%-9s%6.3f\n", $1,$2,$3}' "$db_file" >> ./DataBaseT.${type} ;;
			    siz)   awk -v d=${db_res} '$2~substr(d,0,2) && substr($2,3,1)~/[0-9]/ {printf"%-6s%-6s%-6.3f\n", $1,$2,$3}' "$db_file" >> ./DataBaseT.${type} ;;
			esac
		    fi
		else
		    case $type in
			crg)   awk -v d=${db_res} '$2~d {printf"%-6s%-9s%6.3f\n", $1,$2,$3}' "$db_file" >> ./DataBaseT.${type} ;;
			siz)   awk -v d=${db_res} '$2~d {printf"%-6s%-6s%-6.3f\n", $1,$2,$3}' "$db_file" >> ./DataBaseT.${type} ;;
		    esac
		fi
	    else
		message E  "Residue identifier $db_res is not present in Delphi database file $db_file. Program will crash"
	    fi
	done
	## Other option to avoid duplicate lines in database is using an associative array in awk##
	# awk '!a[$0]++' ./DataBaseT.${type} >> ./DataBaseT.${type} but the naming for the databases on the for cycles
	# need to be renamed
	
	## check if size of database isn't too large ##
	if [[ `cat ./DataBaseT.${type} | wc -l` -ge 1000 ]] ; then
	    message E  "DataBaseT.${type} exceeds 1000 lines. Fortran will crash!!!"
	fi
    done

    #### Routine to deal with 4 letter residues ####
    #### Identify if there are residue names with 4 letters ####
    # Create an array with the alphabet to create the replace name
    alphabet=({A..Z})
    ## make function to run over all molecules for PB and get only 4 char residue names. ##
    #declare -A res4
    n=0
    echo `awk '/ATOM/ {print substr($0,18,4)}' TMP_protein.pdb | sort | uniq | awk '$1 ~ /^[a-zA-Z0-9]{4}$/ {print}'`

    for res in `awk '/ATOM/ {print substr($0,18,4)}' TMP_protein.pdb | sort | uniq | awk '$1 ~ /^[a-zA-Z0-9]{4}$/ {print}' ` ; do	 
	## Useless since arrays can't be exported. 
	#res4[$res]=$newname
	newname=`echo ${alphabet[$n]}${alphabet[$n]}${alphabet[$n]}`
	
	echo $res $newname >> ./4letterkey.dat

	sed -i "s/${res}/${newname} /" DataBaseT.crg
	sed -i "s/${res}/${newname} /" DataBaseT.siz
	
	#message W "Residue $res will be recognized by the PB/MC cycle as:  $newname "
	((++n))
	echo $n
    done

    ## make CRG_FILE ##
    awk 'NF==3 {printf"%-6s%-9s%6.3f\n", $1,$2,0.000}' ./DataBaseT.crg > CRG_FILE
    
    # Generate charges database for Delphi
    "$DelphiDir"/gen_charge.awk "${CPHMD_ALL_SITES_FILE}" TMP_CpHMD_charge.top TMP_delphi.crg

}

# run_PBMC SITES_FILE PHASE_LABEL
#   Prepare PB coordinates, run Delphi through a fixed short internal basename,
#   convert the electrostatic result and run PETIT. Every command and expected
#   file is validated before protonation-state or RT logic can continue.
#
#   The short basename is essential: legacy Delphi/Fortran filename handling can
#   fail silently for long user-facing SysName values. User output names remain
#   unchanged; only disposable PB/MC files use CPHMD_PBMC_BASENAME (PBMC).
run_PBMC ()
{
    local sites_file="$1"
    local phase_label="${2:-standard}"
    local expected_sites
    local pbmc_name="${CPHMD_PBMC_BASENAME:-PBMC}"
    local rc

    expected_sites="$(cphmd_count_nonblank "${sites_file}")" || {
        cphmd_pbmc_error input "PB/MC site definition ${sites_file} is missing."
        return $?
    }
    (( expected_sites > 0 )) || {
        cphmd_pbmc_error input "PB/MC was called with an empty site definition ${sites_file}."
        return $?
    }

    if [[ ! "${pbmc_name}" =~ ^[A-Za-z0-9_]{1,8}$ ]]; then
        cphmd_pbmc_error input "Internal PB/MC basename '${pbmc_name}' must contain 1-8 letters, digits or underscores."
        return $?
    fi

    rm -f -- TMP_delphi.out TMP_delphi.err TMP_convert.err TMP_MCarlo.out TMP_MCarlo.err \
        "${pbmc_name}.sites" "${pbmc_name}.pkcrg" "${pbmc_name}.g" \
        "${pbmc_name}.dat" "TMP_${pbmc_name}.gro"

    if [ "$memb" == 1 ]; then
        if [[ ${MembraneCenteringProtocol:-tails} == "central_atoms" ]]; then
            echo -e "CentralAtom\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux1.gro -n TMP_CpHMD.ndx \
                -center -pbc atom || {
                    cphmd_pbmc_error coordinates "First central-atom centering step failed."
                    return $?
                }
            echo -e "CentralAtoms\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective_aux1.gro -s TMP_CpHMD.tpr \
                -o TMP_${runname}.gro -n TMP_CpHMD.ndx \
                -center -pbc atom || {
                    cphmd_pbmc_error coordinates "Second central-atom centering step failed."
                    return $?
                }
        else
            echo -e "Onetail\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux1.gro -n TMP_CpHMD.ndx \
                -center -pbc atom || {
                    cphmd_pbmc_error coordinates "Historical membrane-centering step failed."
                    return $?
                }
            echo -e "Monotail\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective_aux1.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux2.gro -n TMP_CpHMD.ndx \
                -center -pbc atom || {
                    cphmd_pbmc_error coordinates "Historical membrane-centering step failed."
                    return $?
                }
            echo -e "Bitail\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective_aux2.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux3.gro -n TMP_CpHMD.ndx \
                -center -pbc atom || {
                    cphmd_pbmc_error coordinates "Historical membrane-centering step failed."
                    return $?
                }
            echo -e "Protein\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective_aux3.gro -s TMP_CpHMD.tpr \
                -o TMP_${runname}.gro -n TMP_CpHMD.ndx \
                -center -pbc atom || {
                    cphmd_pbmc_error coordinates "Historical membrane-centering step failed."
                    return $?
                }
        fi

        if [[ $PBdim -eq 0 ]]; then
            echo -e "Protein" | "$GroDIR" trjconv \
                -f TMP_${runname}.gro -s TMP_CpHMD.tpr \
                -n TMP_CpHMD.ndx -o TMP_aux1.gro -pbc mol -quiet || {
                    cphmd_pbmc_error coordinates "Molecule-whole PB/MC conversion failed."
                    return $?
                }

            if [[ ${MembraneCenteringProtocol:-tails} == "central_atoms" ]]; then
                cp -f TMP_aux1.gro TMP_${runname}.gro || {
                    cphmd_pbmc_error coordinates "Could not finalize the central-atom PB/MC coordinate file."
                    return $?
                }
            else
                rm -f TMP_auxcenter.dat
                while read -r site_number _; do
                    awk -v i="${site_number}" 'substr($0,1,5)+0==i {print substr($0,21,24)}' \
                        TMP_aux1.gro >> TMP_auxcenter.dat
                done < "${sites_file}"

                cphmd_require_nonempty_files "PB/MC centering coordinates are empty" TMP_auxcenter.dat || return $?
                protx="$(awk '{x+=$1;n++}END{if(n) print x/n}' TMP_auxcenter.dat)"
                proty="$(awk '{x+=$2;n++}END{if(n) print x/n}' TMP_auxcenter.dat)"
                halfsizex="$(tail -n 1 TMP_aux1.gro | awk '{print $1/2}')"
                halfsizey="$(tail -n 1 TMP_aux1.gro | awk '{print $2/2}')"
                XCoor="$(awk -v p="${protx}" -v h="${halfsizex}" 'BEGIN{print h-p}')"
                YCoor="$(awk -v p="${proty}" -v h="${halfsizey}" 'BEGIN{print h-p}')"

                "$GroDIR" editconf -f TMP_aux1.gro -o TMP_aux2.gro \
                    -translate "${XCoor}" "${YCoor}" 0 -quiet || {
                        cphmd_pbmc_error coordinates "PB/MC xy translation failed."
                        return $?
                    }
                echo -e "Protein" | "$GroDIR" trjconv \
                    -f TMP_aux2.gro -s TMP_CpHMD.tpr -n TMP_CpHMD.ndx \
                    -o TMP_${runname}.gro -pbc atom -quiet || {
                        cphmd_pbmc_error coordinates "Final PB/MC atom-PBC conversion failed."
                        return $?
                    }
            fi
        fi
    else
        if [[ $multiple_solutes -eq 1 || $include_itp -eq 1 ]]; then
            echo -e "Protein\nProtein" | "$GroDIR" trjconv \
                -f TMP_effective.gro -s TMP_CpHMD.tpr \
                -o TMP_effective_aux1.gro -n TMP_CpHMD.ndx \
                -pbc mol -center || {
                    cphmd_pbmc_error coordinates "Non-membrane molecule-whole centering failed."
                    return $?
                }
            cphmd_validate_fixgro_file "${runname}.fixgro"
            "${CpHDIR}/scripts/fixgro" TMP_effective_aux1.gro \
                "${runname}.fixgro" > "TMP_${runname}.gro" || {
                    cphmd_pbmc_error coordinates "fixgro failed while preparing the PB/MC structure."
                    return $?
                }
        else
            echo "Protein" | "$GroDIR" trjconv -f TMP_effective.gro \
                -o TMP_${runname}.gro -s TMP_CpHMD.tpr \
                -n TMP_CpHMD.ndx -pbc mol -quiet || {
                    cphmd_pbmc_error coordinates "Non-membrane molecule-whole conversion failed."
                    return $?
                }
        fi

        if [[ ${PBMCDebugPDB:-0} -eq 1 ]]; then
            "$GroDIR" editconf -f TMP_${runname}.gro -o TMP_aux_debug.pdb || {
                cphmd_pbmc_error coordinates "PB/MC debug-PDB conversion failed."
                return $?
            }
            cat TMP_aux_debug.pdb >> ${blockname}_PQR.pdb
        fi
    fi

    cphmd_require_nonempty_files "PB/MC coordinate preparation failed" "TMP_${runname}.gro" || return $?

    # Use a fixed short internal basename for every legacy Delphi/PETIT file.
    cp -f -- "${sites_file}" "${pbmc_name}.sites" || {
        cphmd_pbmc_error input "Could not stage the short-name PB/MC site file."
        return $?
    }
    cp -f -- "TMP_${runname}.gro" "TMP_${pbmc_name}.gro" || {
        cphmd_pbmc_error input "Could not stage the short-name PB/MC coordinate file."
        return $?
    }

    "$DelphiDir"/delphiT "$nCPU" "${pbmc_name}" "$PBdim" \
        > TMP_delphi.out 2> TMP_delphi.err
    rc=$?
    if (( rc != 0 )); then
        cphmd_pbmc_error delphi "delphiT exited with status ${rc} during ${phase_label} PB/MC."
        return $?
    fi
    if grep -qiE 'FORTRAN STOP|segmentation fault|fatal' TMP_delphi.err; then
        cphmd_pbmc_error delphi "delphiT reported a fatal diagnostic during ${phase_label} PB/MC."
        return $?
    fi
    cphmd_require_nonempty_files "delphiT did not produce required electrostatic files" \
        "${pbmc_name}.pkcrg" "${pbmc_name}.g" || {
            cphmd_pbmc_error delphi \
                "delphiT returned success but did not produce ${pbmc_name}.pkcrg and ${pbmc_name}.g."
            return $?
        }

    "$DelphiDir"/convert "${pbmc_name}.pkcrg" "${pbmc_name}.g" "$temp" \
        > "${pbmc_name}.dat" 2> TMP_convert.err
    rc=$?
    if (( rc != 0 )); then
        cphmd_pbmc_error convert "convert exited with status ${rc}."
        return $?
    fi
    cphmd_require_nonempty_files "convert did not produce a PETIT input" "${pbmc_name}.dat" || {
        cphmd_pbmc_error convert "convert returned success but produced no PETIT input."
        return $?
    }

    echo "$pH,$pH,1" -E "$pot,$pot,1" -T "$temp" -c 2 -r "$seed" -q 1000 100000
    "$PetitDIR"/petit -H "$pH,$pH,1" -E "$pot,$pot,1" -T "$temp" \
        -c 2 -r "$seed" -q 1000 100000 < "${pbmc_name}.dat" \
        > TMP_MCarlo.out 2> TMP_MCarlo.err
    rc=$?
    if (( rc != 0 )); then
        cphmd_pbmc_error petit "PETIT exited with status ${rc}."
        return $?
    fi

    cphmd_validate_petit_output TMP_MCarlo.out "${expected_sites}" || return $?

    if [[ $write_states == "y" && $ReduceTitration == 1 ]]; then
        gawk -v c="$Cycle" '
            /^\./ && $4!~"tot" {s=$4; gsub($1" +"$2" +"$3" +"$4,"");p[s]=$0}
            /^>/ {n[$3]=$2}
            END {print "# Cycle "c; for(s in n) printf "%-13s %s\n",n[s],p[s]; print "#"}
        ' TMP_MCarlo.out >> ${runname}.pocc_RT
        awk '/^f/ {print}' TMP_MCarlo.out >> ${runname}.pocc_RT
    fi

    rm -f -- "${pbmc_name}.sites" "${pbmc_name}.pkcrg" "${pbmc_name}.g" \
        "${pbmc_name}.dat" "TMP_${pbmc_name}.gro" \
        "${pbmc_name}.summ" "${pbmc_name}.pkint" "${pbmc_name}.out" \
        "${pbmc_name}.pqr" "${pbmc_name}.pqr2" *.potat TMP_aux*

    message W "Validated ${phase_label} PB/MC result for ${expected_sites} site(s)."
}

# update_topology
#   Apply the full state vector prepared for this cycle to every all-site entry.
#   In RT mode the vector already merges fixed and active sites, so topology
#   updates are complete and consistent even when the active subset is empty.
update_topology ()
{
    local all_count
    local state_count

    cphmd_require_nonempty_files "topology update inputs are incomplete" \
        TMP_CpHMD.top "${CPHMD_ALL_SITES_FILE}" "${CPHMD_RT_CURRENT_OCC_FILE}" || return $?

    all_count="$(cphmd_count_nonblank "${CPHMD_ALL_SITES_FILE}")" || return 42
    state_count="$(awk '{print NF; exit}' "${CPHMD_RT_CURRENT_OCC_FILE}")"
    [[ "${all_count}" == "${state_count}" ]] || {
        cphmd_output_error \
            "Topology update has ${state_count} states for ${all_count} all-site definitions."
        return $?
    }

    mv TMP_CpHMD.top TMP_CpHMD-pre.top || return 42
    awk '{print $1, substr($2,1,3)}' "${CPHMD_ALL_SITES_FILE}" > aux_sites.sites
    awk '{for (i=1; i<=NF; i++) print $i}' "${CPHMD_RT_CURRENT_OCC_FILE}" > tmp_states.out
    paste tmp_states.out aux_sites.sites > states_file.dat

    cphmd_require_nonempty_files "topology state mapping is empty" states_file.dat || return $?

    "$CpHDIR"/scripts/update-top states_file.dat \
        ${ffID}.ff/protstates.dic TMP_CpHMD-pre.top > TMP_CpHMD.top || {
            cphmd_output_error "update-top failed for cycle ${Cycle}."
            return $?
        }
    cphmd_require_nonempty_files "update-top produced no topology" TMP_CpHMD.top || return $?
}

# run_dynamics PREFIX
#   Run one effective MD phase using PREFIX.mdp/.top/.gro. It removes stale
#   outputs, runs grompp and the prebuilt CPHMD_MDRUN_CMD, and requires fresh
#   non-empty GRO/EDR/XTC/LOG/TPR outputs before returning success.
run_dynamics ()
{
    local stage="$1"
    local input_stage="$2"
    local rc
    local plumed_file=""
    local -a md_args

    # A stale output can make a failed mdrun look successful. Remove every file
    # that this invocation is expected to create before grompp/mdrun starts.
    rm -f -- \
        "TMP_${stage}.tpr" "TMP_${stage}_out.mdp" \
        "TMP_${stage}.xtc" "TMP_${stage}.gro" "TMP_${stage}.edr" \
        "TMP_${stage}.log" "TMP_${stage}.trr"

    if ! "$GroDIR" grompp -f "TMP_${stage}.mdp" -po "TMP_${stage}_out.mdp" \
        -c "TMP_${input_stage}.gro" -p TMP_CpHMD.top -pp TMP_processed.top \
        -n TMP_CpHMD.ndx -o "TMP_${stage}.tpr" -maxwarn 1000 -quiet
    then
        cphmd_output_error "grompp failed for the ${stage} stage in cycle ${Cycle}."
        return $?
    fi

    cphmd_require_nonempty_files "grompp did not create a usable TPR" \
        "TMP_${stage}.tpr" || return $?

    md_args=(
        -s "TMP_${stage}.tpr"
        -x "TMP_${stage}.xtc"
        -c "TMP_${stage}.gro"
        -e "TMP_${stage}.edr"
        -g "TMP_${stage}.log"
        -o "TMP_${stage}.trr"
        -nice 19
    )

    if [[ ${stage} == "effective" && ${plumed:-0} -eq 1 ]]; then
        case ${plumedtype} in
            grid)
                if (( Cycle == 1 )); then
                    sed '/RESTART /d' "${runname}_plumed.dat" > "${runname}_first.dat"
                    sed -i "s/GRID_RFILE=${grid_name}//" "${runname}_first.dat"
                    plumed_file="${runname}_first.dat"
                else
                    plumed_file="${runname}_plumed.dat"
                fi
                ;;
            hill)
                if (( Cycle == 1 )); then
                    sed '/RESTART /d' "${runname}_plumed.dat" > "${runname}_first.dat"
                    sed -i "s/GRID_RFILE=${grid_name}//" "${runname}_first.dat"
                    plumed_file="${runname}_first.dat"
                else
                    plumed_file="${runname}_plumed.dat"
                fi
                ;;
            static)
                plumed_file="${runname}_plumed.dat"
                ;;
            *)
                cphmd_output_error "Unsupported plumedtype '${plumedtype}'."
                return $?
                ;;
        esac
        md_args+=( -plumed "${plumed_file}" )
    fi

    # CPHMD_MDRUN_CMD is assembled once in CpHMD.sh as a Bash array. This avoids
    # re-parsing resource options and still permits CPU/GPU settings overrides.
    "${CPHMD_MDRUN_CMD[@]}" "${md_args[@]}"
    rc=$?
    if (( rc != 0 )); then
        cphmd_output_error "mdrun failed for the ${stage} stage in cycle ${Cycle} (mdrun exit code ${rc})."
        return $?
    fi

    # Effective dynamics must always provide the files consumed by data_append.
    # Relaxation does not necessarily write an XTC because RelaxSteps can be
    # shorter than nstxout-compressed, so only its restart-relevant outputs are
    # mandatory.
    if [[ ${stage} == "effective" ]]; then
        cphmd_require_nonempty_files "effective MD output contract failed in cycle ${Cycle}" \
            "TMP_${stage}.tpr" "TMP_${stage}.gro" "TMP_${stage}.edr" \
            "TMP_${stage}.xtc" "TMP_${stage}.log" || return $?
    else
        cphmd_require_nonempty_files "relaxation MD output contract failed in cycle ${Cycle}" \
            "TMP_${stage}.tpr" "TMP_${stage}.gro" "TMP_${stage}.edr" \
            "TMP_${stage}.log" || return $?
    fi

    if [[ ${stage} == "effective" && ${plumed:-0} -eq 1 ]]; then
        if grep -q "PLUMED error" "TMP_${stage}.log"; then
            cphmd_output_error "PLUMED reported an error in cycle ${Cycle}."
            return $?
        fi

        if [[ ${plumedtype} == "grid" ]]; then
            if [[ ! -s ${hills} ]]; then
                cphmd_output_error "PLUMED grid run did not create ${hills} in cycle ${Cycle}."
                return $?
            fi
            mv -f -- "${hills}" "${hills}_${Cycle}"
        fi
    fi

    rm -f \#*
}

# run_relaxation
#   Perform the short solvent/structure relaxation that follows a protonation
#   change, then prepare the relaxed coordinates for the effective MD phase.
run_relaxation ()
{
    # Solvent relaxation.
    run_dynamics relax effective || return $?

    # Prepare input GRO for effective dynamics. The solvent marker remains
    # historical behaviour and is not altered by the output-safety patch.
    local SOL1st="SOL"
    awk -v s="${SOL1st}" '$1 ~ s {exit};{print $0}' TMP_effective.gro > TMP_aux.gro
    awk -v s="${SOL1st}" '$1 ~ s {a=1};a' TMP_relax.gro >> TMP_aux.gro
    mv -f TMP_aux.gro TMP_relax.gro

    cphmd_require_nonempty_files "failed to prepare relaxed coordinates" TMP_relax.gro
}

# data_append
#   Merge the current cycle into the scheduler-segment accumulators. EDR and XTC
#   concatenation is transactional and retains the required final extensions so
#   GROMACS writes the expected filenames. It also appends logs, occupation and
#   optional pulling/PLUMED outputs, then verifies the complete cycle contract.
data_append ()
{
    # Build the accumulated EDR/XTC transactionally. The previous accumulated
    # file is never used as the output path, so a failed GROMACS tool cannot
    # destroy the last valid segment data.
    local segment_start="${CPHMD_SEGMENT_START_CYCLE:-1}"
    local current_start
    # GROMACS determines output type from the final filename extension.
    local transaction_id="${BASHPID:-$$}"
    local edr_tmp="TMP_CpHMD.edr.next.${transaction_id}.edr"
    local xtc_tmp="TMP_CpHMD.xtc.next.${transaction_id}.xtc"
    local -a edr_inputs xtc_inputs set_times

    cphmd_require_nonempty_files "data_append received incomplete effective MD output" \
        TMP_effective.gro TMP_effective.tpr TMP_effective.edr \
        TMP_effective.xtc TMP_effective.log || return $?

    if (( Cycle == segment_start )); then
        InitTime=`echo $sim_time-$WriteTime | bc -l`
        edr_inputs=(TMP_effective.edr)
        xtc_inputs=(TMP_effective.xtc)
        set_times=("${InitTime}")
    else
        cphmd_require_nonempty_files "previous accumulated segment output is missing" \
            TMP_CpHMD.edr TMP_CpHMD.xtc || return $?
        current_start=`echo $sim_time-$WriteTime | bc -l`
        edr_inputs=(TMP_CpHMD.edr TMP_effective.edr)
        xtc_inputs=(TMP_CpHMD.xtc TMP_effective.xtc)
        set_times=("${InitTime}" "${current_start}")
    fi

    rm -f -- "${edr_tmp}" "${xtc_tmp}"

    if ! printf '%s\n' "${set_times[@]}" | \
        "$GroDIR" eneconv -o "${edr_tmp}" -f "${edr_inputs[@]}" -settime -quiet
    then
        rm -f -- "${edr_tmp}" "${xtc_tmp}"
        cphmd_output_error "eneconv failed while appending cycle ${Cycle}."
        return $?
    fi
    cphmd_require_nonempty_files "eneconv produced no accumulated energy file" "${edr_tmp}" || {
        rm -f -- "${edr_tmp}" "${xtc_tmp}"
        return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
    }

    if ! printf '%s\n' "${set_times[@]}" | \
        "$GroDIR" trjcat -f "${xtc_inputs[@]}" -o "${xtc_tmp}" -settime
    then
        rm -f -- "${edr_tmp}" "${xtc_tmp}"
        cphmd_output_error "trjcat failed while appending cycle ${Cycle}."
        return $?
    fi
    cphmd_require_nonempty_files "trjcat produced no accumulated trajectory" "${xtc_tmp}" || {
        rm -f -- "${edr_tmp}" "${xtc_tmp}"
        return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"
    }

    mv -f -- "${edr_tmp}" TMP_CpHMD.edr
    mv -f -- "${xtc_tmp}" TMP_CpHMD.xtc

    # Append only the log created by this effective run. Stale backup logs were
    # removed before mdrun, so duplicate GROMACS records cannot enter the segment.
    cat TMP_effective.log >> TMP_CpHMD.log || {
        cphmd_output_error "Could not append the effective MD log in cycle ${Cycle}."
        return $?
    }

    cp -f TMP_effective.gro TMP_CpHMD.gro || return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"

    # Keep the TPR corresponding to the latest effective dynamics/topology. The
    # historical code retained the initial TMP_CpHMD.tpr instead, which was not a
    # faithful endpoint record after protonation-state changes.
    cp -f TMP_effective.tpr TMP_CpHMD.tpr || return "${CPHMD_EXIT_OUTPUT_FAILURE:-42}"

    cphmd_require_nonempty_files "accumulated CpHMD output contract failed after cycle ${Cycle}" \
        TMP_CpHMD.gro TMP_CpHMD.top TMP_CpHMD.tpr TMP_CpHMD.edr \
        TMP_CpHMD.xtc TMP_CpHMD.log TMP_CpHMD.occ TMP_CpHMD.mocc || return $?

    # Append HILLS for GRID setup.
    if [[ -f ${hills}_${Cycle} && ${plumedtype:-} == "grid" ]]; then
        cat "${hills}_${Cycle}" >> "${hills}_curr_seg"
    fi

    # Append pulling outputs when present.
    if [[ -f pullx.xvg ]]; then
        local e
        for e in x f; do
            awk -v t=`echo "$EffectiveSteps*$TimeStep*($Cycle-1)" | bc -l` \
                '!/^\#/ && !/^\@/ && $1!=0.0000{print $1+t,$2}' "pull${e}.xvg" \
                >> "TMP_CpHMD_pull${e}.xvg"
        done
    fi

    rm -f TMP_aux* \#* TMP_effective*.log pull?.xvg
}

# cphmd_write_checkpoint FINISHED_CYCLE
#   Snapshot the latest complete cycle into segments/, restart/ and log-files/.
#   In RT mode the active subset and fixed-state templates are part of the
#   authoritative restart package and are copied atomically on every cycle.
cphmd_write_checkpoint ()
{
    local finished_cycle="$1"
    local segment_start="${CPHMD_SEGMENT_START_CYCLE:-1}"
    local segment_name="${blockname}_cycles${segment_start}-${finished_cycle}"
    local segment_dir="../${CPHMD_SEGMENT_DIR:-segments}"
    local restart_dir="../${CPHMD_RESTART_DIR:-restart}"
    local log_dir="../${CPHMD_LOG_DIR:-log-files}"
    local segment_base="${segment_dir}/${segment_name}"
    local log_base="${log_dir}/${segment_name}"
    local e

    mkdir -p -- "${segment_dir}" "${restart_dir}" "${log_dir}" || {
        cphmd_output_error "Could not create organized CpHMD output directories."
        return $?
    }

    cphmd_require_nonempty_files "refusing to checkpoint an incomplete CpHMD cycle ${finished_cycle}" \
        TMP_CpHMD.gro TMP_CpHMD.top TMP_CpHMD.tpr TMP_CpHMD.edr \
        TMP_CpHMD.xtc TMP_CpHMD.log TMP_CpHMD.occ TMP_CpHMD.mocc \
        "${blockname}.info" || return $?

    find "${segment_dir}" -maxdepth 1 -type f \
        -name "${blockname}_cycles${segment_start}-*" -delete 2>/dev/null || true
    find "${log_dir}" -maxdepth 1 -type f \
        -name "${blockname}_cycles${segment_start}-*" -delete 2>/dev/null || true

    cphmd_copy_atomic TMP_CpHMD.gro "${restart_dir}/${blockname}.gro" || return $?
    cphmd_copy_atomic TMP_CpHMD.top "${restart_dir}/${blockname}.top" || return $?
    cphmd_copy_atomic TMP_CpHMD.tpr "${restart_dir}/${blockname}.tpr" || return $?

    if (( ${ReduceTitration:-0} == 1 )); then
        [[ -e "${CPHMD_RT_ACTIVE_SITES_FILE}" ]] || {
            cphmd_output_error "Missing active RT site set at checkpoint cycle ${finished_cycle}."
            return $?
        }
        cphmd_require_nonempty_files "missing RT templates at checkpoint cycle ${finished_cycle}" \
            "${CPHMD_RT_TEMPLATE_OCC_FILE}" "${CPHMD_RT_TEMPLATE_MOCC_FILE}" || return $?
        cphmd_copy_atomic "${CPHMD_RT_ACTIVE_SITES_FILE}" \
            "${restart_dir}/${blockname}.RT-active.sites" || return $?
        cphmd_copy_atomic "${CPHMD_RT_TEMPLATE_OCC_FILE}" \
            "${restart_dir}/${blockname}.RT-template.occ" || return $?
        cphmd_copy_atomic "${CPHMD_RT_TEMPLATE_MOCC_FILE}" \
            "${restart_dir}/${blockname}.RT-template.mocc" || return $?
    else
        rm -f -- "${restart_dir}/${blockname}.RT-active.sites" \
            "${restart_dir}/${blockname}.RT-template.occ" \
            "${restart_dir}/${blockname}.RT-template.mocc"
    fi

    for e in gro top tpr edr xtc occ mocc; do
        cphmd_copy_atomic "TMP_CpHMD.${e}" "${segment_base}.${e}" || return $?
    done

    cphmd_copy_atomic TMP_CpHMD.log "${log_base}.gromacs.log" || return $?
    cphmd_copy_atomic "${blockname}.info" "${log_base}.info" || return $?
    if [[ -s ${blockname}_MD.log ]]; then
        cphmd_copy_atomic "${blockname}_MD.log" "${log_base}.driver.log" || return $?
    fi

    [[ -e ${CPHMD_ALL_SITES_FILE} ]] && \
        cphmd_copy_atomic "${CPHMD_ALL_SITES_FILE}" "../${CPHMD_ALL_SITES_FILE}" || true

    if [[ -f Eb_calculation.dat ]]; then
        cphmd_copy_atomic Eb_calculation.dat "${segment_base}.ene" || return $?
    fi

    if [[ -f TMP_CpHMD_pullx.xvg ]]; then
        for e in x f; do
            [[ -f TMP_CpHMD_pull${e}.xvg ]] && \
                cphmd_copy_atomic "TMP_CpHMD_pull${e}.xvg" "${segment_base}_pull${e}.xvg" || true
        done
    fi

    if [[ -f ${CPHMD_RT_HISTORY_FILE} ]]; then
        cphmd_copy_atomic "${CPHMD_RT_HISTORY_FILE}" "../${SysName}_RT-sites.dat" || return $?
    fi

    if (( ${CPHMD_KEEP_RT_DEBUG:-0} == 1 )) && [[ -f ${runname}.pocc_RT ]]; then
        cphmd_copy_atomic "${runname}.pocc_RT" "${log_dir}/${SysName}_RT-debug.pocc_RT" || return $?
    fi

    if [[ ${plumed:-0} == "1" ]]; then
        [[ -f ${colvar_name} ]] && cphmd_copy_atomic "${colvar_name}" "../${colvar_name}" || true
        if [[ -f ${hills} && ${plumedtype} != "grid" ]]; then
            cphmd_copy_atomic "${hills}" "../${hills}" || return $?
        fi
        if [[ -f ${hills}_curr_seg && ${plumedtype} == "grid" ]]; then
            cphmd_copy_atomic "${hills}_curr_seg" "../${hills}" || return $?
        fi
        [[ -f ${grid_name} && ${plumedtype} == "grid" ]] && \
            cphmd_copy_atomic "${grid_name}" "../${grid_name}" || true
    fi

    cphmd_write_state_atomic "${finished_cycle}" 0 validated_cycle_checkpoint
}

# cphmd_finalize_segment FINISHED_CYCLE
#   Commit a scheduler segment: verify required files, run `gmx check` on XTC and
#   EDR, compress configured immutable outputs, remove optional logs when asked,
#   and finally write checkpoint.state with FINALIZED=1. A later job may resume
#   only after this function succeeds.
cphmd_finalize_segment ()
{
    local finished_cycle="$1"
    local segment_start="${CPHMD_SEGMENT_START_CYCLE:-1}"
    local segment_name="${blockname}_cycles${segment_start}-${finished_cycle}"
    local segment_dir="../${CPHMD_SEGMENT_DIR:-segments}"
    local restart_dir="../${CPHMD_RESTART_DIR:-restart}"
    local log_dir="../${CPHMD_LOG_DIR:-log-files}"
    local segment_base="${segment_dir}/${segment_name}"
    local log_base="${log_dir}/${segment_name}"
    local validation_log="${log_base}.validation.log"

    mkdir -p -- "${segment_dir}" "${restart_dir}" "${log_dir}" || {
        cphmd_output_error "Could not create organized CpHMD output directories."
        return $?
    }

    cphmd_require_nonempty_files "cannot finalize incomplete scheduler segment ${segment_start}-${finished_cycle}" \
        "${segment_base}.gro" "${segment_base}.top" "${segment_base}.tpr" \
        "${segment_base}.edr" "${segment_base}.xtc" \
        "${segment_base}.occ" "${segment_base}.mocc" \
        "${log_base}.gromacs.log" "${log_base}.info" \
        "${restart_dir}/${blockname}.gro" "${restart_dir}/${blockname}.top" \
        "${restart_dir}/${blockname}.tpr" || return $?

    if (( ${ReduceTitration:-0} == 1 )); then
        [[ -e "${restart_dir}/${blockname}.RT-active.sites" ]] || {
            cphmd_output_error "Finalized RT segment is missing ${blockname}.RT-active.sites."
            return $?
        }
        cphmd_require_nonempty_files "finalized RT restart templates are incomplete" \
            "${restart_dir}/${blockname}.RT-template.occ" \
            "${restart_dir}/${blockname}.RT-template.mocc" || return $?
    fi

    # Full binary validation occurs once per scheduler segment.
    if (( ${CPHMD_VALIDATE_SEGMENT_BINARIES:-1} == 1 )); then
        (
            echo "=== XTC validation: $(date --iso-8601=seconds) ==="
            "$GroDIR" check -f "${segment_base}.xtc"
            xtc_rc=$?
            echo
            echo "=== EDR validation: $(date --iso-8601=seconds) ==="
            "$GroDIR" check -e "${segment_base}.edr"
            edr_rc=$?
            exit $((xtc_rc != 0 ? xtc_rc : edr_rc))
        ) > "${validation_log}" 2>&1

        local check_rc=$?
        if (( check_rc != 0 )); then
            cphmd_output_error "GROMACS binary validation failed for segment ${segment_start}-${finished_cycle}; see ${validation_log}."
            return $?
        fi
    fi

    # XTC is already compressed and EDR remains directly readable. Immutable
    # structures, topology and reproducibility TPR can be compressed safely;
    # restart files remain uncompressed for immediate continuation.
    if (( ${CPHMD_COMPRESS_SEGMENT_LOGS:-1} == 1 )); then
        cphmd_gzip_atomic "${log_base}.gromacs.log" || return $?
        if [[ -s ${log_base}.driver.log ]]; then
            cphmd_gzip_atomic "${log_base}.driver.log" || return $?
        fi
    fi

    if (( ${CPHMD_COMPRESS_SEGMENT_TPR:-1} == 1 )); then
        cphmd_gzip_atomic "${segment_base}.tpr" || return $?
    fi

    if (( ${CPHMD_COMPRESS_SEGMENT_STRUCTURES:-1} == 1 )); then
        cphmd_gzip_atomic "${segment_base}.gro" || return $?
        cphmd_gzip_atomic "${segment_base}.top" || return $?
    fi

    if (( ${CPHMD_KEEP_DRIVER_LOG:-1} == 0 )); then
        rm -f -- "${log_base}.driver.log" "${log_base}.driver.log.gz"
    fi

    if (( ${CPHMD_KEEP_VALIDATION_LOG:-1} == 0 )); then
        rm -f -- "${validation_log}"
    fi

    # Mark the segment committed only after validation and compression succeed.
    cphmd_write_state_atomic "${finished_cycle}" 1 validated_finalized_segment || return $?

    message W "Finalized validated CpHMD segment ${segment_start}-${finished_cycle}."
}

# clean_up
#   Remove temporary PB/MC, GROMACS and force-field working files from the inner
#   CpHMD-run_<pid> directory. Persistent segment/restart/log outputs live one
#   directory above and are intentionally not matched by this cleanup.
clean_up ()
{    
    rm -rf ${runname}.{summ,g,pkint,pkcrg,dat,pqr,out} state*\
        *.st aminoacids.dat FF.dat specbond.dat ff* traj.trr \
        ${runname}_cpu* PBMC* TMP_PBMC* TMP_delphi* TMP_convert* \
        TMP_RT_current.* TMP_{statepqr.in,${runname}.pqr,effective,relax,CpHMD,aux,mead,MCarlo,template,posre,processed,pqr}*

   
}
