# Architecture support

## Version 0.7.0

CpHMD-container 0.7.0 supports x86-64 systems.

The distributed x86-64 image contains:

- a CPU-only GROMACS 2026.3 installation;
- a CUDA-enabled GROMACS 2026.3 installation for NVIDIA GPUs;
- PLUMED support in both GROMACS builds.

CPU jobs use `/gromacs/bin/gmx` and do not require Apptainer `--nv`.

NVIDIA GPU jobs use `/gromacs/bin/gmx-gpu` and require Apptainer `--nv`.

ARM/aarch64 and A64FX/SVE images are not included or supported in version
0.7.0. ARM integration, synchronization with the 0.7 runtime, and native ARM
validation are planned for version 0.8.0.
