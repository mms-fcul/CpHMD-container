#!/usr/bin/env python3
"""
Audit CpHMD reduced-titration decisions and restart-state continuity.

The audit follows one interpretation of the available output:

1. Only cycles represented in *_RT-debug.pocc_RT are audited at cycle level.
2. A cycle also recorded in *_RT-sites.dat is treated as a full refresh.
3. Every other debug cycle is treated as a reduced cycle and must use exactly
   the active subset selected by the most recent preceding full refresh.
4. Historical full-refresh snapshots absent from the debug file are context,
   not missing audit data and not failures.
5. The persistent restart active-site file and checkpoint metadata must match
   the latest full-refresh snapshot in *_RT-sites.dat.

This intentionally provides no compatibility mode requiring historical debug
sections that are not present in the available debug file.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


class AuditInputError(ValueError):
    """Raised when an audit input exists but is malformed or inconsistent."""


@dataclass(frozen=True)
class DebugSite:
    residue_label: str
    residue_key: str
    mocc: float
    state_probabilities: tuple[float, ...]

    @property
    def max_state_probability(self) -> float:
        return max(self.state_probabilities)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit CpHMD reduced-titration full-refresh decisions, reduced-cycle "
            "active-set continuity, and persistent restart state."
        )
    )

    parser.add_argument(
        "--run-dir",
        default=".",
        help="Directory containing CpHMD output files. Default: current directory.",
    )
    parser.add_argument(
        "--sys-name",
        required=True,
        help="CpHMD SysName used in output filenames.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help="Reduced-titration threshold matching RTThreshold. Default: 0.01.",
    )
    parser.add_argument(
        "--from-cycle",
        type=int,
        default=None,
        help=(
            "Audit only debug cycles at or after this cycle. Earlier full-refresh "
            "snapshots remain available as context for reduced cycles."
        ),
    )
    parser.add_argument(
        "--show-fixed",
        action="store_true",
        help=(
            "Print a short summary of residues predicted to be fixed for each "
            "audited full-refresh cycle."
        ),
    )

    args = parser.parse_args()

    if not math.isfinite(args.threshold) or not 0.0 <= args.threshold <= 1.0:
        parser.error("--threshold must be a finite number between 0 and 1 inclusive")

    if args.from_cycle is not None and args.from_cycle < 1:
        parser.error("--from-cycle must be a positive integer")

    return args


def normalize_residue_key(token: str) -> str:
    """Normalize labels such as GLU-25 or 25 to the .sites residue key."""
    token = token.strip()
    match = re.search(r"(\d+)$", token)
    if match:
        return match.group(1)
    return token


def residue_key_from_sites_line(line: str) -> str:
    fields = line.split()
    if not fields:
        raise AuditInputError("encountered an empty .sites record")
    return normalize_residue_key(fields[0])


def natural_key(value: str) -> tuple[int, int | str]:
    if value.isdigit():
        return (0, int(value))
    return (1, value)


def is_terminal_residue(residue_label: str) -> bool:
    """Mirror the runtime rule that NTR and CTR sites are always retained."""
    upper = residue_label.upper()
    return upper.startswith("NTR") or upper.startswith("CTR")


def parse_debug_pocc(path: Path) -> dict[int, dict[str, DebugSite]]:
    """Parse *_RT-debug.pocc_RT into cycle -> residue-key -> DebugSite."""
    cycles: dict[int, dict[str, DebugSite]] = {}
    current_cycle: int | None = None
    cycle_re = re.compile(r"^#\s*Cycle\s+(\d+)\s*$")

    with path.open(encoding="utf-8") as fh:
        for line_number, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue

            cycle_match = cycle_re.match(line)
            if cycle_match:
                current_cycle = int(cycle_match.group(1))
                if current_cycle in cycles:
                    raise AuditInputError(
                        f"{path}: duplicate debug section for cycle {current_cycle} "
                        f"at line {line_number}"
                    )
                cycles[current_cycle] = {}
                continue

            if current_cycle is None:
                continue

            if line == "#" or line.startswith("f "):
                continue

            fields = line.split()
            if len(fields) < 3:
                continue

            residue_label = fields[0]
            residue_key = normalize_residue_key(residue_label)

            try:
                numeric_values = tuple(float(value) for value in fields[1:])
            except ValueError:
                continue

            if len(numeric_values) < 2:
                continue
            if not all(math.isfinite(value) for value in numeric_values):
                raise AuditInputError(
                    f"{path}: non-finite probability data for {residue_label} "
                    f"in cycle {current_cycle} at line {line_number}"
                )
            if any(value < 0.0 or value > 1.0 for value in numeric_values[1:]):
                raise AuditInputError(
                    f"{path}: state probability outside [0, 1] for {residue_label} "
                    f"in cycle {current_cycle} at line {line_number}"
                )
            if residue_key in cycles[current_cycle]:
                raise AuditInputError(
                    f"{path}: duplicate residue key {residue_key!r} in cycle "
                    f"{current_cycle} at line {line_number}"
                )

            cycles[current_cycle][residue_key] = DebugSite(
                residue_label=residue_label,
                residue_key=residue_key,
                mocc=numeric_values[0],
                state_probabilities=numeric_values[1:],
            )

    return cycles


def parse_rt_sites(path: Path) -> dict[int, set[str]]:
    """Parse full-refresh active-subset snapshots from *_RT-sites.dat."""
    cycles: dict[int, set[str]] = {}
    current_cycle: int | None = None
    cycle_re = re.compile(
        r"^Active reduced-titration sites after full update at Cycle =\s*(\d+)\s*$"
    )

    with path.open(encoding="utf-8") as fh:
        for line_number, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue

            cycle_match = cycle_re.match(line)
            if cycle_match:
                current_cycle = int(cycle_match.group(1))
                if current_cycle in cycles:
                    raise AuditInputError(
                        f"{path}: duplicate full-refresh snapshot for cycle "
                        f"{current_cycle} at line {line_number}"
                    )
                cycles[current_cycle] = set()
                continue

            if current_cycle is None or line.startswith("Active sites:"):
                continue

            residue_key = residue_key_from_sites_line(line)
            if residue_key in cycles[current_cycle]:
                raise AuditInputError(
                    f"{path}: duplicate residue key {residue_key!r} in the "
                    f"cycle {current_cycle} snapshot at line {line_number}"
                )
            cycles[current_cycle].add(residue_key)

    return cycles


def parse_sites_keys(path: Path) -> set[str]:
    keys: set[str] = set()

    with path.open(encoding="utf-8") as fh:
        for line_number, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            residue_key = residue_key_from_sites_line(line)
            if residue_key in keys:
                raise AuditInputError(
                    f"{path}: duplicate residue key {residue_key!r} at line {line_number}"
                )
            keys.add(residue_key)

    return keys


def parse_checkpoint(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}

    with path.open(encoding="utf-8") as fh:
        for line_number, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            if "=" not in line:
                raise AuditInputError(
                    f"{path}: malformed checkpoint line {line_number}: {line!r}"
                )
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            if not key:
                raise AuditInputError(
                    f"{path}: empty checkpoint key at line {line_number}"
                )
            if key in values:
                raise AuditInputError(
                    f"{path}: duplicate checkpoint key {key!r} at line {line_number}"
                )
            values[key] = value

    return values


def parse_checkpoint_integer(checkpoint: dict[str, str], key: str) -> int:
    if key not in checkpoint:
        raise AuditInputError(f"checkpoint is missing required key {key}")

    try:
        value = int(checkpoint[key])
    except ValueError as exc:
        raise AuditInputError(
            f"checkpoint key {key} is not an integer: {checkpoint[key]!r}"
        ) from exc

    if value < 0:
        raise AuditInputError(f"checkpoint key {key} must not be negative: {value}")

    return value


def format_site(site: DebugSite | None) -> str:
    if site is None:
        return "not found in debug section"
    return (
        f"{site.residue_label:12s} "
        f"mocc={site.mocc:.8f} "
        f"max_state_prob={site.max_state_probability:.8f}"
    )


def print_first_keys(
    title: str,
    keys: Iterable[str],
    debug_sites: dict[str, DebugSite] | None = None,
) -> None:
    ordered = sorted(keys, key=natural_key)
    if not ordered:
        return

    print(title)
    for key in ordered[:20]:
        if debug_sites is None:
            print(f"    {key}")
        else:
            print(f"    {format_site(debug_sites.get(key))}")


def latest_snapshot_before(
    snapshots: dict[int, set[str]], cycle: int
) -> tuple[int, set[str]] | None:
    eligible = [snapshot_cycle for snapshot_cycle in snapshots if snapshot_cycle < cycle]
    if not eligible:
        return None
    snapshot_cycle = max(eligible)
    return snapshot_cycle, snapshots[snapshot_cycle]


def require_files(paths: Iterable[Path]) -> bool:
    missing = [path for path in paths if not path.exists()]
    if not missing:
        return True

    for path in missing:
        print(f"ERROR: missing required file: {path}", file=sys.stderr)
    return False


def main() -> int:
    args = parse_args()

    run_dir = Path(args.run_dir).resolve()
    sys_name = args.sys_name
    threshold = args.threshold
    fixed_cutoff = 1.0 - threshold

    debug_file = run_dir / "log-files" / f"{sys_name}_RT-debug.pocc_RT"
    rt_sites_file = run_dir / f"{sys_name}_RT-sites.dat"
    all_sites_file = run_dir / f"{sys_name}-all.sites"
    checkpoint_block = f"{sys_name}_CpHrun"
    active_sites_file = run_dir / "restart" / f"{checkpoint_block}.RT-active.sites"
    checkpoint_file = run_dir / "resub_maint" / "checkpoint.state"

    required_files = (
        debug_file,
        rt_sites_file,
        all_sites_file,
        active_sites_file,
        checkpoint_file,
    )
    if not require_files(required_files):
        return 2

    try:
        debug_cycles_all = parse_debug_pocc(debug_file)
        site_cycles = parse_rt_sites(rt_sites_file)
        all_site_keys = parse_sites_keys(all_sites_file)
        persistent_active_keys = parse_sites_keys(active_sites_file)
        checkpoint = parse_checkpoint(checkpoint_file)
        checkpoint_active_count = parse_checkpoint_integer(
            checkpoint, "CPHMD_RT_ACTIVE_COUNT"
        )
        checkpoint_last_full = parse_checkpoint_integer(
            checkpoint, "CPHMD_RT_LAST_FULL_CYCLE"
        )
    except (OSError, AuditInputError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.from_cycle is None:
        debug_cycles = debug_cycles_all
    else:
        debug_cycles = {
            cycle: sites
            for cycle, sites in debug_cycles_all.items()
            if cycle >= args.from_cycle
        }

    print("Reduced-titration audit")
    print("=======================")
    print(f"Run directory          : {run_dir}")
    print(f"System name            : {sys_name}")
    print(f"RT debug file          : {debug_file.relative_to(run_dir)}")
    print(f"RT history file        : {rt_sites_file.relative_to(run_dir)}")
    print(f"All-site file          : {all_sites_file.relative_to(run_dir)}")
    print(f"Persistent active file : {active_sites_file.relative_to(run_dir)}")
    print(f"Checkpoint file        : {checkpoint_file.relative_to(run_dir)}")
    print(f"RTThreshold            : {threshold}")
    print(f"Fixed-state cutoff     : {fixed_cutoff}")
    print(
        "Audit scope            : "
        + (
            "all available debug cycles"
            if args.from_cycle is None
            else f"available debug cycles >= {args.from_cycle}"
        )
    )
    print("Decision rule          : retain if max_state_probability <= 1 - RTThreshold")
    print("                         fix if max_state_probability  > 1 - RTThreshold")
    print("Termini rule           : NTR/CTR sites are always retained")
    print("History interpretation : snapshots absent from debug are context, not failures")
    print("Note                   : mocc is reported but not used for the decision")
    print()

    print("Available state")
    print("---------------")
    print(f"All-site count              : {len(all_site_keys)}")
    print(f"Full-refresh snapshots      : {len(site_cycles)}")
    print(f"Debug cycles in file        : {len(debug_cycles_all)}")
    print(f"Debug cycles in audit scope : {len(debug_cycles)}")
    print(f"Persistent active count     : {len(persistent_active_keys)}")
    print(f"Checkpoint active count     : {checkpoint_active_count}")
    print(f"Checkpoint last full cycle  : {checkpoint_last_full}")

    overall_ok = True

    if checkpoint.get("CPHMD_REDUCED_TITRATION") != "1":
        print(
            "ERROR: checkpoint does not report CPHMD_REDUCED_TITRATION=1 "
            f"(found {checkpoint.get('CPHMD_REDUCED_TITRATION', 'missing')!r})."
        )
        overall_ok = False

    if not site_cycles:
        print("ERROR: RT history contains no full-refresh snapshots.")
        overall_ok = False

    if not debug_cycles:
        print("ERROR: no debug cycles are present in the selected audit scope.")
        overall_ok = False

    print()
    print("RT-history structural consistency")
    print("---------------------------------")
    if not site_cycles:
        print("No full-refresh snapshots are available.")
    else:
        for cycle in sorted(site_cycles):
            snapshot = site_cycles[cycle]
            unexpected = snapshot - all_site_keys
            print(
                f"Cycle {cycle}: {len(snapshot)} active sites; "
                f"outside all-site definition: {len(unexpected)}"
            )
            if unexpected:
                overall_ok = False
                print_first_keys("  First unexpected history keys:", unexpected)

    print()
    print("Persistent active-state consistency")
    print("-----------------------------------")
    if site_cycles:
        latest_full_cycle = max(site_cycles)
        latest_snapshot = site_cycles[latest_full_cycle]
        missing_persistent = latest_snapshot - persistent_active_keys
        unexpected_persistent = persistent_active_keys - latest_snapshot

        print(f"Latest recorded full refresh : cycle {latest_full_cycle}")
        print(f"Recorded active sites        : {len(latest_snapshot)}")
        print(f"Persistent active sites      : {len(persistent_active_keys)}")
        print(f"Missing from persistent file : {len(missing_persistent)}")
        print(f"Unexpected persistent sites  : {len(unexpected_persistent)}")

        if missing_persistent or unexpected_persistent:
            overall_ok = False
            print_first_keys(
                "  First sites missing from persistent file:", missing_persistent
            )
            print_first_keys(
                "  First unexpected persistent sites:", unexpected_persistent
            )

        if checkpoint_last_full != latest_full_cycle:
            print(
                "ERROR: checkpoint last-full cycle does not match RT history "
                f"({checkpoint_last_full} != {latest_full_cycle})."
            )
            overall_ok = False

    if checkpoint_active_count != len(persistent_active_keys):
        print(
            "ERROR: checkpoint active count does not match the persistent "
            f"active-site file ({checkpoint_active_count} != "
            f"{len(persistent_active_keys)})."
        )
        overall_ok = False

    print()
    print("Debug section sizes")
    print("-------------------")
    if debug_cycles:
        for cycle in sorted(debug_cycles):
            cycle_type = "full refresh" if cycle in site_cycles else "reduced"
            print(
                f"Cycle {cycle:>4}: {len(debug_cycles[cycle]):>4} residues "
                f"({cycle_type})"
            )
    else:
        print("No debug cycles are in the selected audit scope.")

    full_refresh_cycles = sorted(set(debug_cycles) & set(site_cycles))

    print()
    print("Full-refresh decision audit")
    print("---------------------------")
    if not full_refresh_cycles:
        print("No full-refresh cycle is represented in the selected debug scope.")

    for cycle in full_refresh_cycles:
        debug_sites = debug_cycles[cycle]
        debug_keys = set(debug_sites)
        missing_from_debug = all_site_keys - debug_keys
        unexpected_in_debug = debug_keys - all_site_keys

        predicted_kept: set[str] = set()
        predicted_fixed: set[str] = set()

        for residue_key, site in debug_sites.items():
            if (
                is_terminal_residue(site.residue_label)
                or site.max_state_probability <= fixed_cutoff
            ):
                predicted_kept.add(residue_key)
            else:
                predicted_fixed.add(residue_key)

        observed_kept = site_cycles[cycle]
        missing_from_history = predicted_kept - observed_kept
        unexpected_in_history = observed_kept - predicted_kept

        print(f"Cycle {cycle}")
        print(f"  All-site definitions     : {len(all_site_keys)}")
        print(f"  Debug residues analysed  : {len(debug_sites)}")
        print(f"  Missing from debug       : {len(missing_from_debug)}")
        print(f"  Unexpected in debug      : {len(unexpected_in_debug)}")
        print(f"  Predicted retained       : {len(predicted_kept)}")
        print(f"  Predicted fixed          : {len(predicted_fixed)}")
        print(f"  Recorded retained        : {len(observed_kept)}")
        print(f"  Missing from RT history  : {len(missing_from_history)}")
        print(f"  Unexpected in RT history : {len(unexpected_in_history)}")

        if args.show_fixed:
            print_first_keys(
                "  First predicted fixed residues:",
                predicted_fixed,
                debug_sites,
            )

        if missing_from_debug or unexpected_in_debug:
            overall_ok = False
            print_first_keys(
                "  First all-site residues missing from debug:",
                missing_from_debug,
                debug_sites,
            )
            print_first_keys(
                "  First unexpected debug residues:",
                unexpected_in_debug,
                debug_sites,
            )

        if missing_from_history or unexpected_in_history:
            overall_ok = False
            print_first_keys(
                "  First predicted residues missing from RT history:",
                missing_from_history,
                debug_sites,
            )
            print_first_keys(
                "  First unexpected RT-history residues:",
                unexpected_in_history,
                debug_sites,
            )

        print()

    reduced_cycles = sorted(set(debug_cycles) - set(site_cycles))

    print("Reduced-cycle active-set continuity")
    print("-----------------------------------")
    if not reduced_cycles:
        print("No reduced cycles are present in the selected debug scope.")
    else:
        for cycle in reduced_cycles:
            context = latest_snapshot_before(site_cycles, cycle)
            if context is None:
                print(
                    f"Cycle {cycle}: ERROR: no preceding full-refresh snapshot "
                    "defines the expected active subset."
                )
                overall_ok = False
                continue

            full_cycle, expected_active = context
            debug_sites = debug_cycles[cycle]
            observed_active = set(debug_sites)
            missing_from_debug = expected_active - observed_active
            unexpected_in_debug = observed_active - expected_active

            print(
                f"Cycle {cycle}: expected {len(expected_active)} sites from full "
                f"refresh cycle {full_cycle}; observed {len(observed_active)}"
            )
            print(f"  Missing active sites    : {len(missing_from_debug)}")
            print(f"  Unexpected active sites : {len(unexpected_in_debug)}")

            if missing_from_debug or unexpected_in_debug:
                overall_ok = False
                print_first_keys(
                    "  First expected sites missing from debug:",
                    missing_from_debug,
                    debug_sites,
                )
                print_first_keys(
                    "  First unexpected debug sites:",
                    unexpected_in_debug,
                    debug_sites,
                )

    print()
    if overall_ok:
        print(
            "PASS: available full-refresh decisions, reduced-cycle active sets, "
            "and persistent RT restart state are mutually consistent."
        )
        return 0

    print(
        "FAIL: reduced-titration decisions, active-set continuity, or persistent "
        "restart state are inconsistent."
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
