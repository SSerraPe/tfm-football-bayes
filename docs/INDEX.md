# docs/ — documentation index

All written documentation for the project now lives under `docs/` (consolidated 2026-06-26).
`CLAUDE.md` and `README.md` remain at the project root by convention.

| Path | What it is |
|---|---|
| `PROJECT_RUNDOWN_AND_RUNTIME.txt` | Plain-language process walkthrough + runtime deep-dive + per-issue fix notes + Claude-Project upload list. |
| `journal/pipeline_journal.qmd` | Canonical narrative log (Entries 1–4). Read the latest entry when resuming. |
| `professor/professor_summary.qmd` (+ `.pdf`) | Polished results summary for the professor. **PDF is a tracked, shareable deliverable.** |
| `thesis/thesis_technical_document.qmd` (+ `.pdf`, `references.bib`) | Formal thesis write-up. **PDF is a tracked, shareable deliverable.** |
| `model/additive_model.{md,tex}`, `model/modelling_steps_03_to_08.md` | Technical model documentation. |
| `ideas/ideas.Rmd`, `ideas/TODO.Rmd` | Professor's ideas / TODO sources. |
| `analysis_notes/ANALYSIS_INDEX.md`, `analysis_notes/note_*.md`, validation writeups | Hand-authored interpretation notes. |

## Conventions

- **Rendered PDFs are gitignored** except `professor_summary.pdf` and
  `thesis_technical_document.pdf` (the two shareable deliverables). Regenerate any other PDF on
  demand, e.g. `quarto render docs/journal/pipeline_journal.qmd`.
- The three render-source `.qmd` files set their **knit root two levels up** to the project root,
  so their `outputs/...` figure/table includes resolve after the move into `docs/`.
- **Generated** stage notes (written by pipeline scripts) stay in `outputs/notes/` — they are
  script outputs, not hand-authored docs.
- Legacy/historical material is under `archive/` (not active).
