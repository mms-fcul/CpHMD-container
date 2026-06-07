#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# ASIC-support static checks.
#
# This script checks that the optional ASIC/generalized support added on top
# of the reduced-titration container is present and backwards-compatible.
#
# It verifies:
#   - default settings exist for new ASIC-related options
#   - settings template documents those options
#   - functions.sh uses these options in the expected locations
#   - historical/default behavior remains available
# ---------------------------------------------------------------------------

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "==> Running ASIC-support checks from: ${repo_root}"
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
functions="CpHMD/scripts/functions.sh"

require_file "${settings}"
require_file "${template}"
require_file "${functions}"

echo "==> Checking ASIC-related default settings"

require_pattern 'UseInputOffset' "${settings}" \
    "UseInputOffset default setting"

require_pattern 'MembraneCenteringProtocol' "${settings}" \
    "MembraneCenteringProtocol default setting"

require_pattern 'PBMCDebugPDB' "${settings}" \
    "PBMCDebugPDB default setting"

echo "    ASIC-related defaults are present."
echo

echo "==> Checking ASIC-related template documentation"

require_pattern 'UseInputOffset' "${template}" \
    "UseInputOffset template documentation"

require_pattern 'MembraneCenteringProtocol' "${template}" \
    "MembraneCenteringProtocol template documentation"

require_pattern 'PBMCDebugPDB' "${template}" \
    "PBMCDebugPDB template documentation"

require_pattern 'central_atoms' "${template}" \
    "central_atoms documentation"

require_pattern 'tails' "${template}" \
    "tails documentation"

echo "    ASIC-related template documentation is present."
echo

echo "==> Checking optional offset logic"

require_pattern 'UseInputOffset' "${functions}" \
    "UseInputOffset logic in functions.sh"

require_pattern 'Using user-defined PB/MC offset from settings' "${functions}" \
    "user-defined offset reporting message"

require_pattern 'Using automatically determined PB/MC offset' "${functions}" \
    "automatic offset reporting message"

require_pattern 'last_res' "${functions}" \
    "automatic offset fallback"

require_pattern 'export offset' "${functions}" \
		"offset export logic"

echo "    Offset logic is present."
echo

echo "==> Checking optional membrane centering logic"

require_pattern 'MembraneCenteringProtocol' "${functions}" \
    "MembraneCenteringProtocol logic in functions.sh"

require_pattern 'central_atoms' "${functions}" \
    "central_atoms centering branch"

require_pattern 'CentralAtom' "${functions}" \
    "CentralAtom index group use"

require_pattern 'CentralAtoms' "${functions}" \
    "CentralAtoms index group use"

require_pattern 'Onetail' "${functions}" \
    "historical Onetail centering branch"

require_pattern 'Monotail' "${functions}" \
    "historical Monotail centering branch"

require_pattern 'Bitail' "${functions}" \
    "historical Bitail centering branch"

echo "    Membrane-centering logic is present."
echo

echo "==> Checking optional PB/MC debug PDB logic"

require_pattern 'PBMCDebugPDB' "${functions}" \
    "PBMCDebugPDB logic in functions.sh"

require_pattern 'TMP_aux_debug\.pdb' "${functions}" \
    "debug PDB output generation"

require_pattern '_PQR\.pdb' "${functions}" \
    "debug PDB append output"

echo "    PB/MC debug PDB logic is present."
echo

echo "==> Checking improved Delphi database diagnostics"

require_pattern 'db_file=' "${functions}" \
    "db_file variable in make_delphi_DB"

require_pattern 'DatabaseDIR is not defined' "${functions}" \
    "missing DatabaseDIR diagnostic"

require_pattern 'Delphi database file not found' "${functions}" \
    "missing database-file diagnostic"

require_pattern 'Using Delphi database file' "${functions}" \
    "database file reporting"

require_pattern 'Residue identifier .*database file' "${functions}" \
    "specific missing-residue diagnostic"

echo "    Delphi database diagnostics are present."
echo

echo "==> Checking backwards-compatible defaults"

# These checks are intentionally simple. They verify that the new settings
# are present and appear to default to non-ASIC behavior unless explicitly
# requested by the user.
require_pattern 'UseInputOffset=.*0|UseInputOffset="?0"?' "${settings}" \
    "UseInputOffset defaulting to disabled"

require_pattern 'MembraneCenteringProtocol=.*tails|MembraneCenteringProtocol="?tails"?' "${settings}" \
    "MembraneCenteringProtocol defaulting to tails"

require_pattern 'PBMCDebugPDB=.*0|PBMCDebugPDB="?0"?' "${settings}" \
    "PBMCDebugPDB defaulting to disabled"

echo "    Backwards-compatible defaults look OK."
echo

echo "==> ASIC-support checks completed successfully."
