# Phase 6 ARM milestone

Confirmed on Deucalion ARM:
- Petit rebuilt natively for ARM
- GROMACS rebuilt natively for ARM
- geom_center rebuilt natively for ARM
- calc_pkg_low rebuilt natively for ARM
- delphi_pgf90 replaced with ARM DelPhi solver build
- NVHPC runtime libraries bundled into container
- Delphi runtime binaries execute successfully inside ARM container

Current status:
- ARM container port completed at runtime level
- next step is first minimal end-to-end CpHMD smoke test
