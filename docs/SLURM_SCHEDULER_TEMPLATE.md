# Generic SLURM scheduler template

CpHMD-container provides two user-facing scheduler-template files:

- `template/slurm/CpHMD-slurm.settings`
- `template/slurm/generate_CpHMD_slurm.sh`

The scheduler settings file contains cluster-specific paths, partitions,
allocation settings, bind mounts, and scratch policy. The generator reads
`GPU`, `nCPU`, and `SysName` from the simulation's `01_CpHMD.settings` and
creates a complete SLURM script.

The generator does not submit by default:

```bash
cp template/slurm/CpHMD-slurm.settings ./my-run-slurm.settings
# Edit my-run-slurm.settings.
template/slurm/generate_CpHMD_slurm.sh ./my-run-slurm.settings
Review the generated script, then submit it with the printed sbatch command.

The generated job:

selects CPU or CUDA mode from GPU=0 or GPU=1;

adds --nv only for CUDA execution;

uses the runtime embedded in the released SIF by default;

runs active CpHMD files on scheduler scratch mounted as /out;

preserves compact restart and checkpoint state;

preserves reduced-titration history and optional diagnostic history;

preserves HILLS and COLVAR continuation files;

synchronizes finalized segments transactionally;

automatically resubmits only after clean-boundary exit code 10;

does not resubmit output-contract or PB/MC failures.

Input paths referenced in 01_CpHMD.settings must be made visible inside the
container through CPHMD_BIND_PATHS.
