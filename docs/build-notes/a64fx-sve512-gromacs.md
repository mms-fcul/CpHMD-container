# A64FX SVE-512 GROMACS build

This branch enables ARM SVE SIMD for the ARM/A64FX GROMACS 2024.3 build.

Relevant CMake options:
-DGMX_SIMD=ARM_SVE
-DGMX_SIMD_ARM_SVE_LENGTH=512

Deucalion ARM node check:
lscpu
cat /proc/sys/abi/sve_default_vector_length

Expected runtime check:
singularity exec CpHMD-aarch64-phase6-sve512.sif /gromacs/bin/gmx --version

Expected output includes:
SIMD instructions: ARM_SVE
