# Organized CpHMD output layout — PB/MC/RT v6

```text
run-directory/
├── 01_CpHMD.settings
├── <SysName>-all.sites
├── <SysName>_RT-sites.dat              # only when RT is enabled
├── segments/
├── restart/
│   ├── <block>.gro
│   ├── <block>.top
│   ├── <block>.tpr
│   ├── <block>.RT-active.sites         # may be empty
│   ├── <block>.RT-template.occ
│   └── <block>.RT-template.mocc
├── log-files/
├── resub_maint/
└── slurm_files/
```

`log-files/` remains flat. PB/MC failures create files such as:

```text
<block>_PBMC-failure_cycle6_stage-delphi_pid1234.log
failure_job-5125.log
```

The old persistent `<SysName>.sites` name is removed because it did not explain
that the file represented the active subset selected at the last RT refresh.
