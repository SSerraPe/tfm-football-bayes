# PLAN.md — Advisor revision pass

One row (block) per task, in the order given in the prompt. Each block: file, line range, the
exact intended edit (old text quoted verbatim, new text as currently drafted), and a risk flag.
**No edits have been made yet.** See `QUESTIONS.md` for the handful of items that needed
flagging rather than guessing.

---

## Task 1 — Abstract: residual sentence
**File:** `main.tex` **Lines:** 152–156 (part of a longer sentence spanning 152–156)
**Risk:** mechanical

Old (the clause to replace runs from "residuals follow" through "specification."):
> "Each player-season observation is decomposed additively into a stable player effect, a season
> effect, and a residual term. The player effect follows a low-rank factor structure, and
> residuals follow a Student-$t$ distribution whose scale is further adjusted for minutes played,
> correcting a systematic under-fit of low-minute seasons observed under a constant-variance
> specification."

New:
> "Each player-season observation is decomposed additively into a stable player effect, a season
> effect, and a residual term. The player effect follows a low-rank factor structure. Residuals
> are given a Student-$t$ rather than a Gaussian distribution, since thin tails let a handful of
> extreme player-seasons distort the player and season effects; the residual scale is weighted by
> minutes played, reflecting the larger sampling error of per-90 rates computed over few minutes."

---

## Task 2 — Abstract: identifiability-constraint sentence
**File:** `main.tex` **Lines:** 156–158
**Risk:** judgment (one sentence, but wording needs care so the surrounding paragraph still flows)

Old:
> "A standard identifiability constraint on the factor loadings is resolved after fitting through
> a rotation-invariant transformation, recovering interpretable canonical directions of style."

New (deleted; at most a half-clause survives, folded into the next sentence rather than standing
alone — see Task 3 below, which now opens the sentence this one used to precede):
> *(deleted — no replacement sentence; the following sentence, rewritten in Task 3, absorbs
> "we report identified summaries of the factor structure" if needed for flow)*

---

## Task 3 — Abstract: posterior-diagnostics + rank sentences
**File:** `main.tex` **Lines:** 158–161
**Risk:** mechanical

Old:
> "Posterior diagnostics show excellent convergence, and the fitted covariance reproduces the
> observed cross-feature correlation structure closely ($r = 0.975$). A three-factor solution is
> retained on a qualitative variance-share and loading-reliability argument, after comparison
> with alternative ranks."

New:
> "Three factors are retained on a qualitative argument."

*(Combined with Task 2's deletion, the paragraph's second half becomes: "...The player effect
follows a low-rank factor structure. Residuals are given a Student-$t$ rather than a Gaussian
distribution, since thin tails let a handful of extreme player-seasons distort the player and
season effects; the residual scale is weighted by minutes played, reflecting the larger sampling
error of per-90 rates computed over few minutes. Three factors are retained on a qualitative
argument." — I dropped the identifiability-constraint sentence entirely rather than keeping a
half-clause, since "Three factors are retained on a qualitative argument" reads cleanly on its
own and doesn't need it for grammar. Flagging this choice in QUESTIONS.md.)*

---

## Task 4 — Abstract: similarity-score sentence
**File:** `main.tex` **Lines:** 168–169
**Risk:** mechanical

Old:
> "Building on the same shared covariance structure, the thesis introduces a bounded similarity
> score between any two players, calibrated against the population's own typical style distance."

New:
> "Building on the same shared covariance structure, the thesis introduces a bounded similarity
> score between any two players."

---

## Task 5 — Introduction: AI disclosure
**File:** `main.tex` **Line:** 217
**Risk:** mechanical (zero cross-references to this paragraph, confirmed repo-wide)

Old:
> "Part of the R and Stan code supporting this thesis --- model implementation, debugging, and
> pipeline refactoring --- was developed with the assistance of an AI coding assistant (Claude,
> Anthropic), used under the author's direction and review throughout. The statistical modelling
> decisions, the analysis, and the interpretation of results are the author's own."

New:
> "\textbf{AI disclosure:} An AI assistant (Claude, Anthropic) was used in preparing this thesis:
> for R and Stan implementation, debugging and pipeline refactoring, and for drafting and editing
> prose. All output was reviewed and revised by the author. The statistical modelling decisions,
> the analysis, and the interpretation of results are the author's own."

---

## Task 6 — Chapter 4: "cannot quantify uncertainty" motivation
**File:** `methodology.tex` **Line:** 8 (one paragraph)
**Risk:** judgment (self-contained rewrite, same length)

Old:
> "What the exploratory analyses cannot do, however, is quantify uncertainty in the latent
> structure. PCA scores and factor loadings are point estimates computed from the sample; they
> carry no natural notion of sampling variability, posterior distribution, or credible interval.
> This matters especially in a longitudinal setting, where the same player is observed in
> multiple seasons. The exploratory analyses treat all player--season rows as exchangeable,
> whereas a proper account of the data-generating process recognises that two observations
> belonging to the same player are not independent: they share the player's stable latent
> characteristics. Ignoring this structure risks confusing player persistence with seasonal
> noise, and overstating the precision of latent-dimension estimates."

New (draft — argues the model's value is stating a believed structure, not uncertainty
quantification per se):
> "What the exploratory analyses do not do, however, is commit to a structure for how a player's
> statistics arise. PCA scores and factor loadings summarise the sample as given, treating every
> player--season row as exchangeable and saying nothing about why two rows from the same player
> should resemble each other more than two rows from different players. The model developed in
> this chapter states that structure directly: an additive decomposition that separates a stable
> player effect from a season effect, with a residual specification chosen to match how the data
> actually behave rather than assumed by default. This matters especially in a longitudinal
> setting, where the same player is observed in multiple seasons and two such observations are
> not independent: they share the player's stable latent characteristics, which the exploratory
> analyses' exchangeability assumption does not recognise."

---

## Task 7 — Chapter 4 §9: "not a continuous model parameter"
**File:** `methodology.tex` **Line:** 334 (opening sentence of §9 "Rank selection")
**Risk:** mechanical

Old:
> "The factor rank $K$ governs the dimension of the latent space used to model the player
> covariance $\boldsymbol{\Sigma}_a$. Unlike the parameters $\boldsymbol{\Lambda}$,
> $\boldsymbol{\Psi}_a$, and $\nu$, the rank $K$ is not a continuous model parameter that can be
> assigned a prior and sampled within a single model: it determines the parametric family itself.
> The selection of $K$ is therefore a model comparison problem."

New:
> "The factor rank $K$ governs the dimension of the latent space used to model the player
> covariance $\boldsymbol{\Sigma}_a$. Unlike the parameters $\boldsymbol{\Lambda}$,
> $\boldsymbol{\Psi}_a$, and $\nu$, $K$ could in principle be given a prior and estimated jointly
> with the rest of the model; here it is instead fixed by a qualitative, post-hoc argument, since
> $K$ determines the parametric family itself and a within-model treatment was not pursued."

*(Note: this changes methodology.tex:9's stale "§8" numbering context slightly — see
QUESTIONS.md item 1.)*

---

## Task 8 — Chapter 4: delete the K=2/K=3-family passage
**File:** `methodology.tex` **Lines:** 360–364 (full paragraph)
**Risk:** mechanical — **verified exact boundaries, clean deletion**

Old (delete in full):
> "$K=2$ and $K=3$ below are evaluated on the production, minutes-scaled, fixed-$\phi=1/2$
> family; $K=4$ is carried over from the constant-scale family (Decision 2 of the revision that
> corrected this chapter: refitting $K=4$ under minutes-scaling was not undertaken, since the
> existing constant-scale $K=4$ fit already required substantially relaxed sampler settings and
> multi-day runtime for reasons unrelated to residual scaling). The three numbers move only
> marginally between families."

New: *(deleted, no replacement — the paragraph immediately following, line 365, already opens
with "Applying Criteria A and B across $K=2, 3, 4$ gives the following picture," exactly matching
the restart point specified in the prompt. K=4's provenance is stated elsewhere — see
`results.tex:145–169` — so per the prompt's own condition, no footnote is added.)*

---

## Task 9 — Chapter 4: delete Remark 4.19 and its forward-references
**File:** `methodology.tex` **Risk:** mechanical, but touches three locations

**9a. Delete the remark itself — lines 369–372:**
> "\begin{remark}[Comparability with the exploratory factor count of Chapter~\ref{chap:background}]
> \label{rem:k_not_comparable}
> Chapter~\ref{chap:background} reports an eleven-factor solution from classical,
> maximum-likelihood factor analysis of the raw player--season feature matrix. ... This reasoning
> is offered as an explanation for the direction of the discrepancy; it has not been separately
> validated ..., and should be read as such.
> \end{remark}"

New: *(deleted in full)*

**9b. Trim the forward-reference clause — line 113 (part of a longer sentence):**

Old (full sentence):
> "The rank $K$ actually adopted for $\boldsymbol{\Sigma}_a$ is selected separately in
> Section~\ref{sec:rank}, and Remark~\ref{rem:k_not_comparable} explains why it need not, and
> does not, match the raw factor count above."

New:
> "The rank $K$ actually adopted for $\boldsymbol{\Sigma}_a$ is selected separately in
> Section~\ref{sec:rank}."

**9c. Delete the lead-in paragraph — line 343** (flagged in QUESTIONS.md as a necessary
consequence beyond the one location the prompt names):
> "One caveat is worth flagging briefly before making that argument: the rank $K^*=3$ selected
> here is not expected to match the eleven-factor count found in Chapter~\ref{chap:background}'s
> exploratory analysis. The two describe different things, for reasons given in full in
> Remark~\ref{rem:k_not_comparable} below."

New: *(deleted in full)*

---

## Task 10 — Chapter 4: delete §10.1
**File:** `methodology.tex` **Lines:** 410–412
**Risk:** mechanical (confirmed near-duplicate of Remark 4.14, which is kept — externally
referenced from `results.tex:40,315`)

Old:
> "\subsection{Why the identification constraint does not yield interpretable directions}
>
> The LT-PD constraint of Definition~\ref{def:ltpd} resolves the rotational non-identifiability
> by fixing a particular coordinate system in the $K$-dimensional factor space, but it does so in
> a way that is mathematically convenient rather than interpretable. ... The columns of the fitted
> LT-PD loading matrix therefore cannot be reported or interpreted as football archetypes."

New: *(deleted in full — Remark 4.14, `methodology.tex:268–271`, already makes this point and
survives unchanged)*

---

## Task 11 — Chapter 4: dissolve §10 into §6
**File:** `methodology.tex` **Risk:** judgment — flagged per the prompt's own request to confirm
before running

**Delete** the section wrapper — lines 407–408, 414 (headings/labels only, content below moves):
> "\section{Canonical factor directions via spectral decomposition}
> \label{sec:pca_postprocessing}"
> ...
> "\subsection{Spectral decomposition of the common variance matrix}
> \label{sec:spectral}"

**Move**, appended as a new final subsection of §6 "Model identification"
(`methodology.tex:221–271`), titled *"Identified summaries of the factor structure"* (keeping
label `sec:spectral` on this new subsection):
- Trimmed intro (was line 417, one sentence — trim per "two or three sentences" instruction):
  Old: "A rotation-invariant representation of the factor structure is obtained by working
  directly with the common variance matrix $\mathbf{C} = \boldsymbol{\Lambda}\boldsymbol{\Lambda}^\top$,
  which is uniquely determined by the model regardless of which identification constraint was
  used."
  New (draft, ~2-3 sentences): "We report the spectral decomposition of the common variance
  matrix $\mathbf{C} = \boldsymbol{\Lambda}\boldsymbol{\Lambda}^\top$ rather than the raw loadings
  $\boldsymbol{\Lambda}$ themselves, because $\mathbf{C}$ is invariant to the LT-PD identification
  constraint (Proposition~\ref{prop:spectral_inv} below): it is uniquely determined by the model
  regardless of which constraint was used to fit $\boldsymbol{\Lambda}$. What follows are
  therefore identified summaries of the factor structure, not an interpretive rotation."
- Definition 4.22 (lines 419–436) — moved verbatim, label `def:pca_postproc` unchanged.
- Proposition 4.23 + proof (lines 438–444) — moved verbatim, label `prop:spectral_inv` unchanged.
- Line 446 ("Proposition~\ref{prop:spectral_inv} is the key justification...") — folded into the
  trimmed intro above rather than kept as a separate sentence.
- Remark 4.24 (lines 448–451) — moved verbatim, label `rem:sign` unchanged (zero external refs,
  confirmed).

**Repoint** two internal self-references from `\ref{sec:pca_postprocessing}` to
`\ref{sec:spectral}`:
- `methodology.tex:270` (inside Remark 4.14): "...recovered through the post-processing
  procedure described in Section~\ref{sec:pca_postprocessing}." → "...Section~\ref{sec:spectral}."
- `methodology.tex:349` (inside Remark 4.18): "...introduced in
  Section~\ref{sec:pca_postprocessing}." → "...Section~\ref{sec:spectral}."

External references (`results.tex:28,316,318,510` for `def:pca_postproc`/`prop:spectral_inv`) use
labels, not section numbers — no edit needed there; only the rendered theorem *numbers* change,
automatically.

---

## Task 12 — Fix section numbering (Chapters 5/6, all chapters really)
**File:** `main.tex` (preamble) **Risk:** mechanical, one line

**Root cause** (confirmed by reading `amsbook.cls`): the section counter is already
chapter-scoped (`\newcounter{section}[chapter]`), but the *display* macro is
`\thesection = \arabic{section}` — no chapter prefix. No `\numberwithin` needed.

Add to `main.tex`'s preamble (near the other small class-behavior tweaks already there):
```latex
\renewcommand{\thesection}{\thechapter.\arabic{section}}
```
`\thesubsection` already builds on `\thesection` in `amsbook.cls`, so subsections inherit the fix
automatically (e.g. "5.3.1" not "3.1"). No `\flushbottom`/`\raggedbottom` is active anywhere in
the document (confirmed), so there's no vertical-justification interaction.

---

## Task 13 — De-duplicate rank-selection numbers (Ch. 4 keeps argument, Ch. 5 keeps numbers)
**File:** `methodology.tex` **Line:** 365 **Risk:** judgment

Old:
> "Applying Criteria A and B across $K=2, 3, 4$ gives the following picture. The last canonical
> factor's common-variance share declines from $30.7\%$ at $K=2$ to $10.4\%$ at $K=3$ to $7.1\%$
> at $K=4$---a continued but decelerating decline, not a sharp drop to a noise floor at any point
> (Criterion A). Mean residual exceedance under $\tau_{99}$ is $0.77\%$, $0.79\%$, and $0.79\%$ at
> $K=2,3,4$ respectively---all comfortably below the $1\%$ floor and essentially flat, so residual
> adequacy holds at every rank examined and does not discriminate among them (Criterion B)."

New (numbers stripped, qualitative shape kept — table/exact figures live in
`results.tex:130–173` per the resolution):
> "Applying Criteria A and B across $K=2, 3, 4$ gives the following picture. The last canonical
> factor's common-variance share declines steadily from $K=2$ to $K=3$ to $K=4$---a continued but
> decelerating decline, not a sharp drop to a noise floor at any point (Criterion A). Mean
> residual exceedance under $\tau_{99}$ stays below the $1\%$ floor and essentially flat across
> all three ranks, so residual adequacy holds at every rank examined and does not discriminate
> among them (Criterion B)."

*(Line 367, the "Read honestly...honest summary" passage right after, is left untouched by this
task — it's a Phase 3 prose-review item, not a Phase 2 target; see `prose-review.txt`.)*

---

## Task 14 — Move Future work into Chapter 6; rename chapter
**Files:** `results.tex`, `conclusion.tex`, `main.tex` **Risk:** judgment (cross-reference-heavy)

**14a.** Move `results.tex:742–776` (the entire `\section{Future work}`, all four
`\textbf{...}` items, verbatim) to the end of `conclusion.tex` (after its current line 47),
dropping the `\section{Future work}` heading and its label `\label{sec:results-future}` (no
longer needed — see 14d for its replacement target).

**14b.** Rename chapter title — `conclusion.tex:1`:
Old: `\chapter{Conclusion}`
New: `\chapter{Discussion and future work}`
(Label `\label{chap:conclusion}` on line 2 stays unchanged — zero external `\ref` to it today, so
nothing to repoint.)

**14c.** Reconcile the existing 3-of-4 preview sentence — `conclusion.tex:36–39`:

Old:
> "That last result is informative in its own right: a style-only measure should not be expected
> to recover output, and showing exactly where it stops is more useful than a measure that
> quietly conflates the two. It is also the most direct link to the further work identified
> alongside the results themselves, where a companion performance score, uncertainty on the
> similarity scores, and team-level effects are the natural next steps, each building on the
> representation established here rather than replacing it."

New (draft — trims to a short transition rather than a second, competing 3-item summary of the
four items that now follow directly):
> "That last result is informative in its own right: a style-only measure should not be expected
> to recover output, and showing exactly where it stops is more useful than a measure that
> quietly conflates the two. It is also the most direct link to the four directions taken up
> below, each building on the representation established here rather than replacing it."

**14d.** Retarget forward-references in `results.tex`:
- Line 523: `...is future work (\S\ref{sec:results-future}).` → `...is future work
  (Chapter~\ref{chap:conclusion}).`
- Line 694: `...(\S\ref{sec:results-future}) would let this chapter close...` → `...
  (Chapter~\ref{chap:conclusion}) would let this chapter close...`
- Line 730: `...that half of the question is a genuine performance score, taken up as future work
  below, but it is...` → `...that half of the question is a genuine performance score, taken up
  as future work in Chapter~\ref{chap:conclusion}, but it is...`

**14e.** Update the Introduction's roadmap sentence — `main.tex:215` (flagged in
QUESTIONS.md as a necessary consequence, not explicitly named in the prompt):

Old (last clause of the sentence):
> "...Chapter~\ref{chap:results} reports model diagnostics, the resulting canonical factor
> directions and player archetypes, and case studies on player similarity and its scouting
> applications, and closes with the model's limitations and directions for future work,
> connecting this line of work to broader Bayesian modelling in sports
> \cite{baio2010bayesian,ribeiro2025bayesian,swinton2023bayesian}."

New:
> "...Chapter~\ref{chap:results} reports model diagnostics, the resulting canonical factor
> directions and player archetypes, and case studies on player similarity and its scouting
> applications; Chapter~\ref{chap:conclusion} discusses the model's limitations and directions
> for future work, connecting this line of work to broader Bayesian modelling in sports
> \cite{baio2010bayesian,ribeiro2025bayesian,swinton2023bayesian}."

---

## Verification checklist (Phase 2, after edits)

1. Full clean recompile: `rm -f main.aux main.lof main.lot main.toc main.out main.bbl main.blg
   main.pdf && pdflatex -interaction=nonstopmode main.tex && bibtex main && pdflatex
   -interaction=nonstopmode main.tex && pdflatex -interaction=nonstopmode main.tex`.
2. `grep -n "^!" main.log`, `grep -n "undefined" main.log`, `grep -n "??" main.pdf` (via
   `pdftotext`) — zero hits required.
3. Visually confirm (render + inspect): Abstract page; §6's new final subsection (moved §10
   content); Chapter 5's rank-selection section (§ number should read like "5.3"); the renamed
   Chapter 6 with Future work folded in; and that Chapter/section numbers throughout now read
   "X.Y" instead of bare "Y".
