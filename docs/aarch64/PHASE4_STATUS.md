# Phase 4 ARM status

Current ARM container status:
- portable CpHMD assets work
- Petit rebuilt natively for ARM
- GROMACS rebuilt natively for ARM
- Delphi bundled x86 libc removed

Remaining unresolved ARM blockers:
- /Delphi/geom_center
- /Delphi/calc_pkg_low
- /Delphi/delphi_pgf90

Repo audit indicates these binaries are required by delphiT, but their source/build instructions are not present in the public repository.
