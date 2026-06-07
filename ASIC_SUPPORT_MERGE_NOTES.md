# asic-support merge notes

This implementation starts from the git-reduced CpHMD container code and adds the ASIC/large-system support that was developed locally, while keeping historical defaults for existing users.

## Main changes

1. `CpHMD/scripts/CpHMD.sh`
   - Added cycle-level checkpoint/restart support controlled by environment variables:
     - `CPHMD_CHECKPOINT_FILE` default `.cphmd_checkpoint_state`
     - `CPHMD_MAX_SECONDS` default `0` (disabled)
     - `CPHMD_WALL_START_EPOCH`
   - Restarts from `${SysName}_CpHrun.restart.gro` and `.top` when the checkpoint file says cycles are complete.
   - Starts from `FirstCycle = CPHMD_DONE_CYCLES + 1` instead of always starting at cycle 1.
   - Writes checkpoint/segment outputs after each completed cycle through `cphmd_write_checkpoint`.
   - Exits with code `10` at a clean cycle boundary when `CPHMD_MAX_SECONDS` is exceeded.
   - Preserves reduced titration and PLUMED logic from the git-reduced version.

2. `CpHMD/scripts/functions.sh`
   - Added optional membrane centering modes:
     - `MembraneCenteringProtocol="tails"` preserves the historical Onetail/Monotail/Bitail workflow.
     - `MembraneCenteringProtocol="central_atoms"` uses the ASIC-compatible CentralAtom/CentralAtoms workflow.
   - Added optional offset control:
     - `UseInputOffset=0` preserves automatic offset from the last protein residue.
     - `UseInputOffset=1` uses the `offset` value from the settings file.
   - Added `PBMCDebugPDB=0/1` to control PB/MC debug PDB output.
   - Improved `make_delphi_DB` so it validates and prints the exact Delphi database file used.
   - Added segment-aware `data_append` compatible with checkpoint restarts.
   - Added `cphmd_write_checkpoint`, including reduced-titration and PLUMED output preservation.

3. Settings files
   - Added documented defaults for `UseInputOffset`, `MembraneCenteringProtocol`, and `PBMCDebugPDB`.
   - Historical defaults preserve existing behaviour.

## ASIC settings to enable

For the ASIC membrane system, set these in the simulation settings file:

```bash
export memb=1
export UseInputOffset=1
export offset=2000
export MembraneCenteringProtocol="central_atoms"
export PBMCDebugPDB=0
```

The index must contain `CentralAtom` and `CentralAtoms` groups when `MembraneCenteringProtocol="central_atoms"` is used.
