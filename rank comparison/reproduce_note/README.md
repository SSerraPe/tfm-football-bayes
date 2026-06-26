# Reproducing the note

This folder is optional. It reproduces the simulation tables in
`../methodology_report.pdf`. Run R from this directory.

For a quick Stan smoke run:

```powershell
Rscript run_main_simulation.R
Rscript run_supplementary_comparisons.R
$env:PLAYER_MODEL_SUMMARY_PROFILE = "quick"
Rscript summarize_note_tables.R
```

To reproduce the longer runs reported in the note:

```powershell
$env:PLAYER_MODEL_SIM_PROFILE = "thorough"
Rscript run_main_simulation.R
$env:PLAYER_MODEL_TASK_PROFILE = "thorough"
Rscript run_supplementary_comparisons.R
Rscript summarize_note_tables.R
```

The scripts write their files to `../results/reproduction/`. They use the
same two Stan files as the real-data workflow.
