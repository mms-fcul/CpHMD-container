# CpHMD-container 0.7.0

## Main changes

- x86-64 CPU and NVIDIA CUDA execution modes.
- GROMACS 2026.3 in isolated CPU and CUDA installations.
- CUDA runtime 12.6.3.
- PLUMED support in both GROMACS builds.
- Explicit CPU/CUDA launcher selection.
- Reduced-titration checkpoint and restart support.
- Organized, finalized cycle-segment output layout.
- Compact persistent restart state.
- Transactional output validation and synchronization.
- Clean scheduler-boundary exit code 10.
- Output-contract failure code 42.
- PB/MC failure code 43.
- Deterministic cross-job continuation.
- Automatic scheduler self-resubmission.
- Generic, adaptable SLURM scheduler templates.

## Architecture

Version 0.7.0 supports x86-64. ARM/aarch64 and A64FX/SVE are planned for a
subsequent release.

## Validation

The release candidate was validated in CPU and CUDA modes using a private
membrane-protein integration system. Validation covered reduced-titration
threshold decisions, persistent restart state, deterministic continuation,
automatic self-resubmission, and binary readability of finalized GROMACS
trajectory and energy outputs.
