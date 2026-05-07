# AArch64 porting notes

## Goal
Port the CpHMD Singularity container from x86_64 to aarch64 for Deucalion ARM nodes.

## Confirmed facts
- Deucalion ARM nodes report `aarch64`.
- Singularity is available on the ARM partition.
- Current public repo contains container-related build directories.

## To audit
- Definition file location
- GROMACS build source and flags
- Delphi build source and flags
- Petit build source and flags
- Python dependency architecture issues
- Hardcoded x86_64 assumptions
- Embedded apps / runscript behavior

## Tests
- Current x86_64 SIF on ARM
- Minimal ARM smoke build with fakeroot
- Minimal ARM container with CpHMD scripts only
- Add GROMACS
- Add Delphi
- Add Petit
