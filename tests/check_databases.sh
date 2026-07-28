#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Delphi database sanity checks.
#
# This script verifies that the Amber14SBpH Delphi charge/radius databases
# contain entries required by the current CpHMD workflows, including:
#   - standard database presence
#   - CYX support
#   - POPC / CHARMM-GUI split lipid-fragment support
#   - basic format sanity
#
# Important database-format note:
#   These database files may start with comment lines beginning with "!".
#   Therefore, the required database header is not necessarily on line 1.
#   The checks below look for the exact header anywhere in the file and skip
#   comments/header lines during data-line validation.
#
# This script does not validate the physical correctness of each parameter.
# It only catches accidental deletion, malformed files, or missing required
# entries needed by the current workflows.
# ---------------------------------------------------------------------------

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "==> Running database checks from: ${repo_root}"
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

check_minimum_lines() {
    local file="$1"
    local min_lines="$2"
    local n_lines

    n_lines="$(wc -l < "${file}")"

    if (( n_lines < min_lines )); then
        fail "${file} has only ${n_lines} lines; expected at least ${min_lines}."
    fi
}

check_exact_header_once() {
    local file="$1"
    local header="$2"
    local description="$3"
    local count

    count="$(grep -Ec "^${header}$" "${file}" || true)"

    if (( count == 0 )); then
        fail "${file} does not contain the expected ${description} header: ${header}"
    fi

    if (( count > 1 )); then
        fail "${file} contains ${count} ${description} headers; expected exactly 1."
    fi
}

check_three_column_data_lines() {
    local file="$1"

    awk '
        # Skip comments, blank lines, and database header lines.
        /^!/ { next }
        /^#/ { next }
        /^[[:space:]]*$/ { next }
        /^atom__/ { next }

        NF != 3 {
            printf "Malformed line in %s at line %d: expected 3 columns, found %d\n", FILENAME, NR, NF > "/dev/stderr"
            print $0 > "/dev/stderr"
            exit 1
        }

        $3 !~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ {
            printf "Malformed numeric value in %s at line %d: third column is not numeric\n", FILENAME, NR > "/dev/stderr"
            print $0 > "/dev/stderr"
            exit 1
        }
    ' "${file}"
}

crg="Databases/DataBaseT_Amber14SBpH.crg"
siz="Databases/DataBaseT_Amber14SBpH.siz"
crg_header="atom__resnumbc_charge_"
siz_header="atom__res_radius_"

echo "==> Checking database files exist"

require_file "${crg}"
require_file "${siz}"

echo "    Database files exist."
echo

echo "==> Checking database sizes"

# These thresholds are intentionally conservative. They are meant to catch
# empty/truncated files, not enforce an exact database size.
check_minimum_lines "${crg}" 900
check_minimum_lines "${siz}" 900

echo "    Database sizes look reasonable."
echo

echo "==> Checking database headers"

# The header may not be on line 1 because the databases can begin with
# comment lines such as "!crg file modified...". Search for the exact
# expected header anywhere in the file and require it to appear once.
check_exact_header_once "${crg}" "${crg_header}" "charge-database"
check_exact_header_once "${siz}" "${siz_header}" "radius-database"

echo "    Headers look OK."
echo

echo "==> Checking basic database formatting"

check_three_column_data_lines "${crg}"
check_three_column_data_lines "${siz}"

echo "    Basic formatting looks OK."
echo

echo "==> Checking required residue support"

# CYX was required for the ASIC workflow and should exist in both databases.
require_pattern '(^|[[:space:]])CYX([[:space:]]|$)' "${crg}" "CYX charge entries"
require_pattern '(^|[[:space:]])CYX([[:space:]]|$)' "${siz}" "CYX radius entries"

# POPC support may be represented as split CHARMM-GUI residue fragments.
# The exact naming depends on the database convention, so we accept any of
# these expected lipid-related residue identifiers.
require_pattern '(^|[[:space:]])(POPC|PPA|PPC|POL|PA|PC|OL)([[:space:]]|$)' "${crg}" \
    "POPC or split-lipid charge entries"

require_pattern '(^|[[:space:]])(POPC|PPA|PPC|POL|PA|PC|OL)([[:space:]]|$)' "${siz}" \
    "POPC or split-lipid radius entries"

echo "    Required residue support is present."
echo

echo "==> Checking for obvious duplicate data lines"

# Duplicate full data lines are not always fatal in legacy database files, but
# they are usually accidental. Report them as warnings for now so this check
# does not break older valid databases. This can be made fatal later if desired.
for file in "${crg}" "${siz}"; do
    duplicates="$({ grep -Ev '^!|^#|^atom__|^[[:space:]]*$' "${file}" || true; } | sort | uniq -d | wc -l)"
    if (( duplicates > 0 )); then
        echo "    WARNING: ${file} contains ${duplicates} duplicate data line(s)."
    else
        echo "    ${file}: no duplicate data lines detected."
    fi
done

echo

echo "==> Database checks completed successfully."
