#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Static shell checks for the CpHMD container repository.
#
# This script performs fast checks that should run on every commit/PR:
#   - repository layout sanity checks
#   - bash syntax checks for Bash scripts and sourced Bash files
#   - optional shellcheck checks, when shellcheck is installed
#   - executable-bit checks for helper programs/scripts
#
# It intentionally does not run GROMACS, Delphi, PETIT, or the container.
# ---------------------------------------------------------------------------

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "==> Running static checks from: ${repo_root}"
echo

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local file="$1"
    [[ -f "${file}" ]] || fail "Required file not found: ${file}"
}

require_executable() {
    local file="$1"
    [[ -x "${file}" ]] || fail "Required executable file is missing executable permission: ${file}"
}

is_shell_script() {
    local file="$1"

    # True for files with a shell shebang.
    head -n 1 "${file}" 2>/dev/null | grep -Eq '^#!.*\b(bash|sh|env bash|env sh)\b'
}

bash_syntax_check_if_shell() {
    local file="$1"

    if is_shell_script "${file}"; then
        echo "    bash -n ${file}"
        bash -n "${file}"
    else
        echo "    skipping bash -n for non-shell executable: ${file}"
    fi
}

shellcheck_if_shell() {
    local file="$1"

    if is_shell_script "${file}"; then
        shellcheck -x "${file}" || true
    fi
}

echo "==> Checking repository layout"

require_file "CpHMD/scripts/CpHMD.sh"
require_file "CpHMD/scripts/functions.sh"
require_file "CpHMD/CpHMD-default.settings"
require_file "template/CpHMD-advanced.settings"

require_file "CpHMD/scripts/fixgro"
require_file "CpHMD/scripts/make_sites"
require_file "CpHMD/scripts/update-top"

require_file "Databases/DataBaseT_Amber14SBpH.crg"
require_file "Databases/DataBaseT_Amber14SBpH.siz"

echo "    Repository layout looks OK."
echo

echo "==> Checking executable permissions"

require_executable "CpHMD/scripts/CpHMD.sh"
require_executable "CpHMD/scripts/fixgro"
require_executable "CpHMD/scripts/make_sites"
require_executable "CpHMD/scripts/update-top"

echo "    Executable permissions look OK."
echo

echo "==> Running bash syntax checks"

# Main Bash entrypoint.
echo "    bash -n CpHMD/scripts/CpHMD.sh"
bash -n "CpHMD/scripts/CpHMD.sh"

# functions.sh is sourced by CpHMD.sh. It may not have a shebang,
# but it is Bash syntax and should always be checked.
echo "    bash -n CpHMD/scripts/functions.sh"
bash -n "CpHMD/scripts/functions.sh"

# Check additional *.sh scripts, excluding the two already checked above.
for script in CpHMD/scripts/*.sh; do
    [[ -e "${script}" ]] || continue

    case "${script}" in
        CpHMD/scripts/CpHMD.sh|CpHMD/scripts/functions.sh)
            continue
            ;;
    esac

    bash_syntax_check_if_shell "${script}"
done

# Known helper executables. Some are compiled binaries, so only syntax-check
# them if they have a shell shebang.
bash_syntax_check_if_shell "CpHMD/scripts/fixgro"
bash_syntax_check_if_shell "CpHMD/scripts/make_sites"
bash_syntax_check_if_shell "CpHMD/scripts/update-top"

echo "    Bash syntax checks passed."
echo

echo "==> Running shellcheck, if available"

if command -v shellcheck >/dev/null 2>&1; then
    # For now shellcheck is advisory. The existing code base may trigger
    # style warnings that are not fatal. Once the code is cleaned further,
    # remove '|| true' to make shellcheck mandatory.
    shellcheck_if_shell "CpHMD/scripts/CpHMD.sh"

    # functions.sh has no shebang but is Bash code sourced by CpHMD.sh.
    shellcheck -x "CpHMD/scripts/functions.sh" || true

    for script in CpHMD/scripts/*.sh; do
        [[ -e "${script}" ]] || continue

        case "${script}" in
            CpHMD/scripts/CpHMD.sh|CpHMD/scripts/functions.sh)
                continue
                ;;
        esac

        shellcheck_if_shell "${script}"
    done

    shellcheck_if_shell "CpHMD/scripts/fixgro"
    shellcheck_if_shell "CpHMD/scripts/make_sites"
    shellcheck_if_shell "CpHMD/scripts/update-top"

    echo "    Shellcheck completed in advisory mode."
else
    echo "    shellcheck not installed; skipping advisory shellcheck step."
fi

echo
echo "==> Static checks completed successfully."
