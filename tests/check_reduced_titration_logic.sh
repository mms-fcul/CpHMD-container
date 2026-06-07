#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Reduced titration static/control-flow checks.
#
# This script checks that the reduced titration implementation is present and
# internally wired through the expected files.
#
# It does not run a CpHMD simulation. It checks that:
#   - default settings expose ReduceTitration, RTInterval, RTThreshold
#   - CpHMD.sh contains the expected all-site/reduced-site scheduling logic
#   - functions.sh contains the expected reduced titration helper logic
#   - RT debug outputs are still copied at the end of the workflow
# ---------------------------------------------------------------------------

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "==> Running reduced titration logic checks from: ${repo_root}"
echo

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local file="$1"
    [[ -f "${file}" ]] || fail "Required file not found: ${file}"
}

require_pattern() {
    local pattern="$1"
    local file="$2"
    local description="$3"

    if ! grep -Eq "${pattern}" "${file}"; then
        fail "Missing ${description} in ${file}. Pattern: ${pattern}"
    fi
}

settings="CpHMD/CpHMD-default.settings"
template="template/CpHMD-advanced.settings"
main_script="CpHMD/scripts/CpHMD.sh"
functions="CpHMD/scripts/functions.sh"

require_file "${settings}"
require_file "${template}"
require_file "${main_script}"
require_file "${functions}"

echo "==> Checking reduced titration settings"

require_pattern 'ReduceTitration' "${settings}" "ReduceTitration default setting"
require_pattern 'RTInterval' "${settings}" "RTInterval default setting"
require_pattern 'RTThreshold' "${settings}" "RTThreshold default setting"

require_pattern 'ReduceTitration' "${template}" "ReduceTitration template documentation"
require_pattern 'RTInterval' "${template}" "RTInterval template documentation"
require_pattern 'RTThreshold' "${template}" "RTThreshold template documentation"

echo "    Reduced titration settings are present."
echo

echo "==> Checking reduced titration validation in functions.sh"

require_pattern 'ReduceTitration[[:space:]]*==[[:space:]]*1' "${functions}" \
    "ReduceTitration enabled check"

require_pattern 'RTInterval' "${functions}" \
    "RTInterval use in functions.sh"

require_pattern 'RTThreshold' "${functions}" \
    "RTThreshold use in functions.sh"

require_pattern 'Reduced titration has been selected' "${functions}" \
    "reduced titration user-facing warning"

echo "    Validation hooks are present."
echo

echo "==> Checking PB/MC reduced-titration control flow"

require_pattern 'run_PBMC[[:space:]]+red' "${main_script}" \
    "all-site PB/MC call that establishes reduced titration"

require_pattern 'Cycle[[:space:]]*%[[:space:]]*RTInterval' "${main_script}" \
    "cycle interval scheduling"

require_pattern '\$\{runname\}-all\.sites' "${main_script}" \
    "all-sites file handling in CpHMD.sh"

require_pattern '\$\{runname\}\.sites' "${main_script}" \
    "active-sites file handling in CpHMD.sh"

require_pattern 'write_fractions_all_sites' "${main_script}" \
    "all-site fraction writing in CpHMD.sh"

require_pattern 'write_fractions' "${main_script}" \
    "reduced-site fraction writing in CpHMD.sh"

echo "    Main control flow looks OK."
echo

echo "==> Checking reduced-titration helper functions"

require_pattern 'write_fractions_all_sites[[:space:]]*\(\)' "${functions}" \
    "write_fractions_all_sites function"

require_pattern 'write_fractions[[:space:]]*\(\)' "${functions}" \
    "write_fractions function"

require_pattern 'TMP_template_occ' "${functions}" \
    "template occupation file handling"

require_pattern 'TMP_template_mocc' "${functions}" \
    "template mean occupation file handling"

require_pattern 'TMP_MCarlo_mod\.out' "${functions}" \
    "modified Monte Carlo output for fixed-site topology updates"

require_pattern 'reducedtitration\.sites' "${functions}" \
    "reduced titration sites log"

require_pattern 'pocc_RT' "${functions}" \
    "reduced titration debug occupation output"

echo "    Helper functions look OK."
echo

echo "==> Checking topology update behavior"

require_pattern 'update_topology[[:space:]]*\(\)' "${functions}" \
    "update_topology function"

require_pattern '\$\{runname\}-all\.sites' "${functions}" \
    "all-sites topology update path"

require_pattern '\$\{runname\}\.sites' "${functions}" \
    "reduced-sites topology update path"

require_pattern 'states_file\.dat' "${functions}" \
    "states file generation"

require_pattern 'update-top' "${functions}" \
    "update-top helper call"

echo "    Topology update logic is present."
echo

echo "==> Checking final RT output copy hooks"

require_pattern 'cphmd_write_checkpoint' "${main_script}" \
    "checkpoint/finalization call in CpHMD.sh"

require_pattern 'cphmd_write_checkpoint[[:space:]]*\(\)' "${functions}" \
    "cphmd_write_checkpoint function definition"

require_pattern '_RT-sites\.dat' "${functions}" \
    "RT sites output copy in checkpoint/finalization function"

require_pattern '_RT-debug\.pocc_RT' "${functions}" \
    "RT debug occupation output copy in checkpoint/finalization function"

echo "    Final RT output hooks are present."
echo

echo "==> Reduced titration logic checks completed successfully."
