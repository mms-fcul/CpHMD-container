# AArch64 build plan

## Phase 1
Build a minimal ARM-compatible container with:
- CpHMD scripts
- force fields
- state/tautomer files
- databases
- analysis scripts
- Delphi payload copied as-is for inspection only

Excluded in phase 1:
- GROMACS
- GROMACS GPU
- Petit
- x86_64 CUDA libraries
- copied x86 shared libraries

## Goal
Confirm that the container definition is portable and can be built natively on ARM once fakeroot access is enabled.
