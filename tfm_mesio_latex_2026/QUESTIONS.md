# QUESTIONS.md — Advisor revision pass

Four items surfaced during Phase 1 that needed flagging rather than guessing. None of these
block starting Phase 2 on the other tasks — they're scoped narrowly to the items below.

---

## 1. Task 8 cites "§8"/"§8.2"; the actual current section is §9

Your prompt says "In §8 (Rank selection)..." and "In §8.2, delete the passage...". The section
titled `\section{Rank selection}` in `methodology.tex` is currently **§9** (line 331), with the
target subsection `\subsection{Selected rank}` (`sec:rank_selection`) at line 357 — not §8/§8.2.

I matched the location by the exact quoted text (found verbatim, unambiguous), not by the section
number, so Tasks 7 and 8 in `PLAN.md` target the right place regardless. Flagging only in case
the stale number means you actually meant a *different* section that happens to not exist under
that number — let me know if so, otherwise I'll proceed against the content match.

## 2. Deleting Remark 4.19 stray-orphans a second paragraph you didn't name

Your prompt says: "also delete the earlier forward-reference to it around the exploratory-factor
discussion" — I found that reference at `methodology.tex:113` (matches your quoted phrase "need
not, and does not, match the raw factor count" exactly).

There is a **second** reference to Remark 4.19, at `methodology.tex:343` — the opening paragraph
of `\subsection{A qualitative argument for the rank}`, which exists *solely* to say "the full
reasoning is in Remark~\ref{rem:k_not_comparable} below." Once the remark is deleted, this
paragraph points at nothing. I'm treating its deletion as a required, unavoidable consequence
(PLAN.md task 9c) rather than something to leave dangling — `pdflatex` would otherwise report an
undefined reference. Flagging so you know it's happening and can veto it if you'd rather I keep a
trimmed version of that paragraph without the dead reference instead of deleting it outright.

## 3. Moving Future work breaks the Introduction's own roadmap sentence

Not one of your 14 numbered tasks, but a direct consequence of Task 14: `main.tex:215` currently
tells the reader that Chapter 5 (Results) "closes with the model's limitations and directions for
future work." Once Future work moves to Chapter 6, that sentence is simply wrong.

I've drafted a one-clause fix in `PLAN.md` (task 14e) that splits the sentence so Chapter 6 gets
credited with the future-work discussion instead. Flagging in case you'd rather leave the
roadmap paragraph untouched and accept the (small) inaccuracy, or word the fix differently
yourself.

## 4. One reference to `def:pca_postproc` is inside a table caption, not body prose

`results.tex:28` (the caption of Table 1, "Posterior convergence summary...") references
`Definition~\ref{def:pca_postproc}` by label. Since Task 11 keeps that label string unchanged
(only its section location and rendered number change), this caption needs **no edit** — I'm
flagging it only so it isn't mistaken for something that broke during the §10→§6 move when you
review the diff.

---

## Choices I made without asking (documented for visibility, not seeking a veto)

- **Task 2** (delete the identifiability-constraint sentence): the prompt allows "at most a
  half-clause" survivor. I dropped it entirely rather than keeping a half-clause, since "Three
  factors are retained on a qualitative argument" (Task 3's replacement) reads cleanly as the
  next sentence without it. If you'd rather keep a short clause (e.g. "...we report identified
  summaries of the factor structure.") tell me and I'll fold it in.
- **Task 11**'s new subsection title: used your suggested *"Identified summaries of the factor
  structure"* verbatim.
- **Task 14a**: folded the four Future-work items into Chapter 6 with no sub-heading (matching
  "no revision history" spirit of a short chapter) rather than keeping a `\section{Future work}`
  heading inside the renamed chapter — since the chapter title itself now says "and future work",
  a repeated section heading seemed redundant. Easy to add back if you'd prefer the explicit
  heading kept.
