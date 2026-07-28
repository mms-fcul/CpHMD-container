# CpHMD settings and scheduler reference

This document contains the detailed explanations intentionally omitted from the
concise `CpHMD-default.settings` template.

## Optional embedded blocks

The settings template may always contain `#fixgro#`, `#plumed#`, and `#mdp#`
lines. They are comments to the shell and are extracted by `CpHMD.sh` only when
needed.

### `#fixgro#`

The fixgro block is active only when all of the following are true:

- `memb=0`; and
- `multiple_solutes=1` or `include_itp=1`.

Membrane workflows (`memb=1`) use `MembraneCenteringProtocol` instead. Keeping a
fixgro block in a membrane settings file is harmless because it is ignored.
When active, the block must contain system-specific molecule/group definitions;
blank blocks and placeholders such as `TOTALATOM` are rejected before PB/MC.

### `#plumed#`

The PLUMED block is active only when `plumed=1`. With `plumed=0`, it is ignored
and no PLUMED input file is created. When active, `plumedtype` must be `grid`,
`hill`, or `static`. Template tokens are replaced from `colvar_stride`,
`colvar_name`, `hills`, and (for grid mode) `grid_name`. Obvious placeholder
paths or atom selections are rejected before PB/MC; final PLUMED/GROMACS syntax
is exercised by the zero-step `mdrun` preflight.

## GROMACS safety flags

- `CPHMD_STRICT_GROMACS_MODE=1`: CPU mode requires the CPU build and GPU mode
  requires the CUDA build.
- `CPHMD_GPU_PREFLIGHT=1`: GPU mode requires a visible NVIDIA device and valid
  logical GPU ID before PB/MC.
- `CPHMD_MDRUN_PREFLIGHT=1`: initialize the configured TPR with `-nsteps 0`
  before PB/MC.

## Organized output contract

- `segments/`: immutable scientific segment files.
- `restart/`: latest uncompressed restart GRO/TOP/TPR.
- `log-files/`: GROMACS, driver, validation and failure logs.
- `resub_maint/`: checkpoint state and scheduler status markers.
- `slurm_files/`: generated SLURM scripts and scheduler output.

## Scratch policy (scheduler environment)

`CPHMD_SCRATCH_ROOT` selects the host scratch root before the container starts.
It defaults to `$TMPDIR`, then `/tmp`.

`CPHMD_SCRATCH_POLICY` accepts exactly:

- `require-local`: reject detected network filesystems and unknown filesystem
  type. Suitable for clusters that require node-local scratch.
- `prefer-local`: default. Warn for network-backed or unknown scratch, but
  continue. This keeps cloud and shared-filesystem deployments usable.
- `allow-network`: deliberately accept network-backed scratch without warning.

Examples:

```bash
# Local cluster production
export CPHMD_SCRATCH_ROOT=/tmp
export CPHMD_SCRATCH_POLICY=require-local

# Portable/default behavior
export CPHMD_SCRATCH_POLICY=prefer-local

# AWS EFS or FSx for Lustre
export CPHMD_SCRATCH_ROOT=/fsx/cphmd
export CPHMD_SCRATCH_POLICY=allow-network
```

AWS instance store and EBS formatted as `ext4`/`xfs` are accepted as block
filesystems. EFS (`nfs4`) and FSx for Lustre (`lustre`) warn under
`prefer-local`, fail under `require-local`, and run under `allow-network`.
