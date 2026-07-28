#!/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
repo_root="$(cd "${script_dir}/.." 2>/dev/null && pwd)"

main_script="${repo_root}/CpHMD/scripts/CpHMD.sh"
functions="${repo_root}/CpHMD/scripts/functions.sh"
defaults="${repo_root}/CpHMD/CpHMD-default.settings"
advanced="${repo_root}/template/CpHMD-advanced.settings"

failures=0

section() {
    echo
    echo "==> $1"
}

pass() {
    echo "    PASS: $1"
}

fail() {
    echo "ERROR: $1" >&2
    failures=$((failures + 1))
}

require_file() {
    file="$1"
    description="$2"

    if [ -f "${file}" ]; then
        pass "${description}"
    else
        fail "Missing ${description}: ${file}"
    fi
}

require_pattern() {
    pattern="$1"
    file="$2"
    description="$3"

    if grep -Eq -- "${pattern}" "${file}"; then
        pass "${description}"
    else
        fail "${description}. Pattern: ${pattern}"
    fi
}

echo "==> Running reduced titration logic checks from: ${repo_root}"

section "Checking required files"

require_file "${main_script}" "CpHMD controller"
require_file "${functions}" "CpHMD function library"
require_file "${defaults}" "default settings"
require_file "${advanced}" "advanced settings template"

if [ "${failures}" -ne 0 ]; then
    echo
    echo "==> Reduced titration checks failed before source inspection."
    echo "Failures: ${failures}"
    exit 1
fi

section "Checking reduced titration settings"

for variable in ReduceTitration RTInterval RTThreshold; do
    require_pattern \
        "^[[:space:]]*export[[:space:]]+${variable}=" \
        "${defaults}" \
        "default ${variable} setting"

    require_pattern \
        "^[[:space:]]*export[[:space:]]+${variable}=" \
        "${advanced}" \
        "advanced-template ${variable} setting"
done

section "Checking reduced titration helper definitions"

require_pattern \
    '^cphmd_rt_full_refresh_cycle[[:space:]]*\(\)' \
    "${functions}" \
    "full-refresh scheduling helper"

require_pattern \
    '\(cycle[[:space:]]*-[[:space:]]*1\)[[:space:]]*%[[:space:]]*RTInterval[[:space:]]*==[[:space:]]*0' \
    "${functions}" \
    "cycle-1 modulo RTInterval scheduling rule"

require_pattern \
    '^cphmd_initialize_or_restore_rt_state[[:space:]]*\(\)' \
    "${functions}" \
    "persistent reduced-titration initialization and restore helper"

require_pattern \
    '^cphmd_validate_petit_output[[:space:]]*\(\)' \
    "${functions}" \
    "PETIT output validation helper"

require_pattern \
    '^cphmd_refresh_reduced_titration_state[[:space:]]*\(\)' \
    "${functions}" \
    "full-refresh active-set selection helper"

require_pattern \
    '^cphmd_merge_reduced_titration_cycle[[:space:]]*\(\)' \
    "${functions}" \
    "reduced active-set merge helper"

require_pattern \
    '^cphmd_materialize_fixed_rt_cycle[[:space:]]*\(\)' \
    "${functions}" \
    "zero-active-site fixed-state materialization helper"

require_pattern \
    '^cphmd_materialize_nonrt_cycle[[:space:]]*\(\)' \
    "${functions}" \
    "non-reduced all-site materialization helper"

require_pattern \
    '^cphmd_record_current_fractions[[:space:]]*\(\)' \
    "${functions}" \
    "current-cycle fraction recording helper"

section "Checking PB/MC reduced-titration control flow"

require_pattern \
    '^[[:space:]]*cphmd_initialize_or_restore_rt_state' \
    "${main_script}" \
    "reduced-titration state initialization before cycling"

require_pattern \
    '^[[:space:]]*if[[:space:]]+cphmd_rt_full_refresh_cycle[[:space:]]+"\$\{Cycle\}"' \
    "${main_script}" \
    "full-refresh scheduling through cphmd_rt_full_refresh_cycle"

require_pattern \
    '^[[:space:]]*run_PBMC[[:space:]]+.*CPHMD_ALL_SITES_FILE.*all-site refresh' \
    "${main_script}" \
    "explicit all-site PB/MC refresh"

require_pattern \
    '^[[:space:]]*cphmd_refresh_reduced_titration_state' \
    "${main_script}" \
    "active-site selection after full refresh"

require_pattern \
    '^[[:space:]]*run_PBMC[[:space:]]+.*CPHMD_RT_ACTIVE_SITES_FILE.*reduced active-set' \
    "${main_script}" \
    "reduced active-set PB/MC execution"

require_pattern \
    '^[[:space:]]*cphmd_merge_reduced_titration_cycle' \
    "${main_script}" \
    "reduced PETIT result merge"

require_pattern \
    '^[[:space:]]*cphmd_materialize_fixed_rt_cycle' \
    "${main_script}" \
    "fixed-state reuse when the active set is empty"

require_pattern \
    '^[[:space:]]*cphmd_materialize_nonrt_cycle' \
    "${main_script}" \
    "ordinary all-site fraction materialization"

require_pattern \
    '^[[:space:]]*cphmd_record_current_fractions' \
    "${main_script}" \
    "full current-cycle OCC/MOCC recording"

section "Checking persistent reduced-titration state"

require_pattern \
    'CPHMD_RT_LAST_FULL_CYCLE' \
    "${main_script}" \
    "last-full-cycle checkpoint variable in the controller"

require_pattern \
    'CPHMD_RT_ACTIVE_COUNT' \
    "${main_script}" \
    "active-site-count checkpoint variable in the controller"

require_pattern \
    'CPHMD_RT_LAST_FULL_CYCLE' \
    "${functions}" \
    "last-full-cycle state handling in the function library"

require_pattern \
    'CPHMD_RT_ACTIVE_COUNT' \
    "${functions}" \
    "active-site-count state handling in the function library"

require_pattern \
    'CPHMD_RT_HISTORY_FILE' \
    "${functions}" \
    "persistent reduced-titration history handling"

require_pattern \
    'CPHMD_RT_ACTIVE_SITES_FILE' \
    "${functions}" \
    "persistent active-site file handling"

section "Checking threshold and termini decisions"

require_pattern \
    'threshold="\$\{RTThreshold\}"' \
    "${functions}" \
    "RTThreshold supplied to active-site selection"

require_pattern \
    'NTR.*CTR|CTR.*NTR' \
    "${functions}" \
    "terminal-site retention logic"

section "Testing full-refresh scheduling behavior"

helper_definition="$(
    awk '
        /^cphmd_rt_full_refresh_cycle[[:space:]]*\(\)/ {
            capture=1
        }

        capture {
            print
        }

        capture && /^}/ {
            exit
        }
    ' "${functions}"
)"

if [ -z "${helper_definition}" ]; then
    fail "Could not extract cphmd_rt_full_refresh_cycle for behavioral testing."
else
    bash -c '
        '"${helper_definition}"'

        RTInterval=2

        for cycle in 1 3 5; do
            cphmd_rt_full_refresh_cycle "${cycle}" || exit 1
        done

        for cycle in 2 4 6; do
            if cphmd_rt_full_refresh_cycle "${cycle}"; then
                exit 1
            fi
        done

        RTInterval=1

        for cycle in 1 2 3 4; do
            cphmd_rt_full_refresh_cycle "${cycle}" || exit 1
        done

        exit 0
    '
    helper_rc=$?

    if [ "${helper_rc}" -eq 0 ]; then
        pass "full-refresh helper selects the expected cycles"
    else
        fail "full-refresh helper behavioral test failed"
    fi
fi

echo

if [ "${failures}" -eq 0 ]; then
    echo "==> Reduced titration logic checks completed successfully."
    exit 0
fi

echo "==> Reduced titration logic checks failed."
echo "Failures: ${failures}"
exit 1
