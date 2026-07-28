# CpHMD PB/MC and reduced-titration reference — v6

## User-facing site files

- `<SysName>-all.sites`: immutable definition of every titratable site.
- `<SysName>_RT-sites.dat`: history of the active subset after each full-site
  refresh.
- `restart/<block>.RT-active.sites`: current active subset applicable after the
  last full-site refresh. The file may legitimately be empty.
- `restart/<block>.RT-template.occ`: full-length fixed-state template. A dash
  marks a currently active site.
- `restart/<block>.RT-template.mocc`: matching fixed-population template.

The old ambiguous `<SysName>.sites` persistent output is no longer created.
The external `make_sites` helper may create that name temporarily inside the
private working directory, after which CpHMD immediately renames it.

## Full refresh schedule

For `RTInterval=N`, full-site PB/MC runs at cycles:

```text
1, 1+N, 1+2N, ...
```

The implementation uses `(cycle - 1) % RTInterval == 0`, including the special
case `RTInterval=1`.

## PB/MC fail-fast contract

CpHMD uses the fixed internal basename `PBMC` for disposable Delphi/PETIT files.
This avoids legacy Fortran filename failures when `SysName` is long.

Before topology or MD can continue, CpHMD requires:

1. successful PB-coordinate preparation;
2. successful `delphiT` exit and nonempty `PBMC.pkcrg`/`PBMC.g`;
3. successful `convert` exit and nonempty `PBMC.dat`;
4. successful PETIT exit;
5. exactly one PETIT `f` state line with the expected number of sites;
6. one site header and one population record for every expected site.

A failure returns application exit code 43, writes a flat PB/MC diagnostic under
`log-files/`, preserves node-local scratch and is never resubmitted.

## Applied states in RT mode

At a full refresh:

- sites with maximum state population greater than `1-RTThreshold` are fixed in
  their most populated state;
- the remaining sites stay active;
- NTR and CTR sites always stay active;
- the current cycle uses sampled states for active sites and maximum-population
  states for fixed sites.

At intermediate cycles, PETIT states for the active subset are merged into the
full fixed templates before topology updating and OCC/MOCC recording. If the
active subset is legitimately empty, PB/MC is skipped and the validated fixed
state is reused.

## Restart contract

A finalized RT checkpoint records:

```text
CPHMD_REDUCED_TITRATION=1
CPHMD_RT_LAST_FULL_CYCLE=<cycle>
CPHMD_RT_ACTIVE_COUNT=<count>
```

The scheduler requires the three RT files in `restart/`. On continuation, CpHMD
regenerates and compares the all-site definition, restores the active set and
templates, verifies their sizes and resumes the interval without silently
resetting reduced titration.
