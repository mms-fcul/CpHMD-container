# Phase 3 ARM build success

Successfully built and tested on Deucalion ARM:
- Petit native aarch64 binary at /programs/petit1.6.1/petit
- GROMACS 2024.3 native aarch64 binary at /gromacs/bin/gmx

Confirmed:
- `gmx --version` runs successfully
- `petit` runs successfully and prints usage
- CpHMD default settings still point to:
  - /gromacs/bin/gmx
  - /programs/petit1.6.1
