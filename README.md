# CpHMD-container

CpHMD-container provides the stochastic-titration constant-pH molecular
dynamics workflow as an Apptainer SIF image.

The container packages the CpHMD runtime, molecular-dynamics software,
Poisson–Boltzmann and Monte Carlo tools, force fields, tautomer definitions,
databases, preparation utilities, and analysis support needed by the
workflow.

## Current release

The current stable release is **v0.7.0**.

Version 0.7.0 includes:

- Ubuntu 24.04 as the container base;
- GROMACS 2026.3 CPU-only and NVIDIA CUDA installations;
- CUDA runtime 12.6.3;
- PLUMED support in both GROMACS installations;
- the CpHMD runtime and preparation utilities;
- GROMOS54a7, CHARMM36, and AMBER14SB CpHMD force-field support;
- Delphi-based Poisson–Boltzmann calculations;
- PETIT 1.6.1 Monte Carlo calculations;
- persistent reduced-titration state;
- deterministic restart and scheduler continuation;
- organized, finalized output segments;
- portable SLURM templates.

### Architecture support

Version 0.7.0 supports:

- x86-64 CPU execution;
- x86-64 execution with NVIDIA CUDA GPUs.

ARM/aarch64 and A64FX/SVE are not supported by version 0.7.0.

See [Architecture support](docs/ARCHITECTURE_SUPPORT.md) for details.

## Download version 0.7.0

### Release page

- [CpHMD-container v0.7.0 release](https://github.com/mms-fcul/CpHMD-container/releases/tag/v0.7.0)
- [All releases](https://github.com/mms-fcul/CpHMD-container/releases)

### Direct downloads

- [Download the v0.7.0 SIF image](https://github.com/mms-fcul/CpHMD-container/releases/download/v0.7.0/CpHMD-container-v0.7.0-x86_64-gromacs-2026.3-cuda-12.6.3.sif)
- [Download the SHA-256 checksum](https://github.com/mms-fcul/CpHMD-container/releases/download/v0.7.0/CpHMD-container-v0.7.0-x86_64-gromacs-2026.3-cuda-12.6.3.sif.sha256)

### Command-line download

```bash
version="v0.7.0"
image="CpHMD-container-v0.7.0-x86_64-gromacs-2026.3-cuda-12.6.3.sif"
release_base="https://github.com/mms-fcul/CpHMD-container/releases/download/${version}"

container_dir="${HOME}/containers/CpHMD"
mkdir -p "${container_dir}"
cd "${container_dir}"

curl \
    --fail \
    --location \
    --retry 3 \
    --output "${image}" \
    "${release_base}/${image}"

curl \
    --fail \
    --location \
    --retry 3 \
    --output "${image}.sha256" \
    "${release_base}/${image}.sha256"

sha256sum \
    --check \
    "${image}.sha256"
```

A successful checksum verification prints:

```text
CpHMD-container-v0.7.0-x86_64-gromacs-2026.3-cuda-12.6.3.sif: OK
```

Users with GitHub CLI may alternatively run:

```bash
gh release download \
    v0.7.0 \
    --repo mms-fcul/CpHMD-container \
    --pattern 'CpHMD-container-v0.7.0-*.sif' \
    --pattern 'CpHMD-container-v0.7.0-*.sif.sha256'
```

## Verify the container

Define the absolute image path:

```bash
image="${HOME}/containers/CpHMD/CpHMD-container-v0.7.0-x86_64-gromacs-2026.3-cuda-12.6.3.sif"
```

Inspect and test the image:

```bash
apptainer inspect "${image}"

apptainer test \
    --cleanenv \
    "${image}"
```

Check the CPU GROMACS installation:

```bash
apptainer exec \
    --cleanenv \
    "${image}" \
    /gromacs/bin/gmx --version
```

On a host with an available NVIDIA GPU, check the CUDA installation:

```bash
nvidia-smi

apptainer exec \
    --cleanenv \
    --nv \
    "${image}" \
    /gromacs/bin/gmx-gpu --version
```

The container intentionally provides separate launchers:

- CPU: `/gromacs/bin/gmx`
- NVIDIA CUDA: `/gromacs/bin/gmx-gpu`

## Prepare a CpHMD run

A production directory normally contains the prepared molecular system,
CpHMD settings, molecular-dynamics inputs, and any optional PLUMED or
system-specific files.

The settings templates are available at:

- [`CpHMD/CpHMD-default.settings`](CpHMD/CpHMD-default.settings)
- [`template/CpHMD-advanced.settings`](template/CpHMD-advanced.settings)

Detailed settings documentation is available in the
[CpHMD settings reference](docs/CPHMD_SETTINGS_REFERENCE.md).

The preparation wiki remains available for workflow-specific guidance:

- [System preparation](https://github.com/mms-fcul/CpHMD-container/wiki/System-Preparation)
- [CpHMD-container wiki](https://github.com/mms-fcul/CpHMD-container/wiki)

## Run CpHMD directly

Use an absolute path for the SIF and bind the production directory to a
clear path inside the container.

Do not bind a host directory over `/CpHMD`, because version 0.7.0 uses the
validated runtime embedded in the released image.

### CPU execution

Set `GPU=0` in the CpHMD settings file.

```bash
image="${HOME}/containers/CpHMD/CpHMD-container-v0.7.0-x86_64-gromacs-2026.3-cuda-12.6.3.sif"
run_dir="/absolute/path/to/production-directory"
settings="01_CpHMD.settings"

cd "${run_dir}"

apptainer exec \
    --cleanenv \
    --bind "${run_dir}:/work" \
    --pwd /work \
    "${image}" \
    /CpHMD/scripts/CpHMD.sh "./${settings}"
```

### NVIDIA GPU execution

Set `GPU=1` in the CpHMD settings file and add Apptainer's `--nv` option.

```bash
image="${HOME}/containers/CpHMD/CpHMD-container-v0.7.0-x86_64-gromacs-2026.3-cuda-12.6.3.sif"
run_dir="/absolute/path/to/production-directory"
settings="01_CpHMD.settings"

cd "${run_dir}"

apptainer exec \
    --cleanenv \
    --nv \
    --bind "${run_dir}:/work" \
    --pwd /work \
    "${image}" \
    /CpHMD/scripts/CpHMD.sh "./${settings}"
```

The simulation settings and container invocation must agree:

- `GPU=0`: CPU execution without `--nv`;
- `GPU=1`: NVIDIA CUDA execution with `--nv`.

## Run through SLURM

Version 0.7.0 includes a portable scheduler template:

- [`template/slurm/CpHMD-slurm.settings`](template/slurm/CpHMD-slurm.settings)
- [`template/slurm/generate_CpHMD_slurm.sh`](template/slurm/generate_CpHMD_slurm.sh)

Create a local scheduler-settings file:

```bash
cp \
    template/slurm/CpHMD-slurm.settings \
    ./my-run-slurm.settings
```

Edit `my-run-slurm.settings` for the cluster, image location, partitions,
allocation, bind paths, and scratch policy.

Generate the job script:

```bash
template/slurm/generate_CpHMD_slurm.sh \
    ./my-run-slurm.settings
```

The generator does not submit by default. Review the generated script and
then use the `sbatch` command printed by the generator.

See the
[SLURM scheduler-template documentation](docs/SLURM_SCHEDULER_TEMPLATE.md)
for the full workflow.

## Container applications

List the applications embedded in the image:

```bash
apptainer inspect \
    --list-apps \
    "${image}"
```

Display help for the PDB-renaming preparation application:

```bash
apptainer run \
    --app pdb2cphmd \
    "${image}" \
    -h
```

Other embedded applications include utilities for extracting force fields
and tautomer definitions. Use `apptainer inspect --list-apps` to obtain the
authoritative application list for the downloaded image.

## Output and restart behavior

Version 0.7.0 organizes run data into dedicated locations including:

- `segments/` for finalized scientific segments;
- `restart/` for the latest compact restart state;
- `log-files/` for execution and validation logs;
- `resub_maint/` for checkpoint and scheduler state;
- `slurm_files/` for generated scheduler files and outputs.

See the
[organized output-layout documentation](docs/ORGANIZED_CPHMD_OUTPUT_LAYOUT.md)
for details.

## Troubleshooting

### The image cannot be found

Use an absolute SIF path:

```bash
readlink -f "${image}"
```

### Checksum verification fails

Remove both the incomplete SIF and checksum file, download them again, and
repeat `sha256sum --check`.

### The NVIDIA GPU is unavailable

Confirm all three conditions:

```bash
nvidia-smi
grep -E '^[[:space:]]*GPU=' 01_CpHMD.settings
apptainer exec --nv "${image}" /gromacs/bin/gmx-gpu --version
```

GPU execution requires `GPU=1`, a working host NVIDIA driver, and
Apptainer's `--nv` option.

### CPU mode starts the CUDA launcher, or GPU mode starts the CPU launcher

Check the `GPU` setting and do not override the embedded `/gromacs/bin`
launchers with host bind mounts.

## Documentation

- [Architecture support](docs/ARCHITECTURE_SUPPORT.md)
- [CpHMD settings reference](docs/CPHMD_SETTINGS_REFERENCE.md)
- [SLURM scheduler template](docs/SLURM_SCHEDULER_TEMPLATE.md)
- [Organized output layout](docs/ORGANIZED_CPHMD_OUTPUT_LAYOUT.md)
- [PB/MC and reduced-titration reference](docs/CPHMD_PBMC_RT_REFERENCE.md)
- [Version 0.7.0 release notes](docs/RELEASE_NOTES_v0.7.0.md)
- [Project wiki](https://github.com/mms-fcul/CpHMD-container/wiki)

## License

CpHMD-container is distributed under the terms in [LICENSE](LICENSE).
