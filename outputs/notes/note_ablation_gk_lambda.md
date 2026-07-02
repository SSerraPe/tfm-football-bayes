# Ablation Note: GK-Only Exclusion Effect on Lambda_a[10,1]

## Setup

Stage-26 fit: Student-t model, K=2, GK rows excluded (N=4586, I=1529), no other
changes. Same Stan model, same sampler settings (500 warmup / 1000 sampling,
delta=0.9, max_depth=10) as stage-18 production fit.

Tested: does removing 358 GK rows (121 GK players) resolve the Lambda_a[10,1]
(per90_passes, factor 1) mixing failure that drove LLt.1.10 ESS=24?

## Results

### Lambda_a[10,1] per-chain means: pre vs. post GK exclusion

| Chain | Pre-GK | Post-GK | Δ |
|---|---|---|---|
| 1 | −0.4476 | −0.6257 | −0.178 |
| 2 | −0.4325 | −0.6161 | −0.184 |
| 3 | −0.4627 | −0.6146 | −0.152 |
| 4 | −0.4657 | −0.6005 | −0.135 |
| **Grand mean** | **−0.4521** | **−0.6142** | **−0.162** |
| **Between-chain SD** | **0.0153** | **0.0104** | **1.5× reduction** |

### ESS / Rhat

| Quantity | Pre-GK | Post-GK |
|---|---|---|
| Lambda_a[10,1] ESS | ~19 | 30 |
| Lambda_a[10,1] Rhat | ~1.18 | 1.061 |
| LLt.1.10 ESS | 24 | 64 |
| LLt.1.10 Rhat | 1.130 | 1.048 |
| LLt.1.2 ESS | 93 | 78 |

## Verdict: PARTIAL GK EFFECT — mixing not resolved

GK exclusion produced measurable but insufficient improvement:
- Between-chain SD dropped 1.5× (0.0153 → 0.0104)
- LLt.1.10 ESS improved 2.7× (24 → 64)
- Rhat improved (LLt.1.10: 1.130 → 1.048)

**But Lambda_a[10,1] ESS=30 and LLt.1.10 ESS=64 remain far below the 200
threshold. The mixing problem is primarily a fundamental posterior geometry
issue, not GK bimodality.**

Observations:
1. The grand mean shifted substantially (−0.452 → −0.614): GKs with high
   passes/low goals diluted the per90_passes loading on the goals factor. The
   post-GK estimate is more meaningful (outfield-player covariance structure).
2. Chain ordering inverted (pre-GK: chains 3,4 most negative; post-GK: chains
   1,2 most negative), indicating poor mixing on both sides — not a stable
   convergence point, just nearby poor-ESS regions.
3. LLt.1.2 (goals×assists) ESS decreased slightly (93 → 78), suggesting GK
   exclusion had no net benefit for this pair.

## Implications for Task K

GK exclusion is the right step regardless (it improves interpretability and
shifts the loading to its correct outfield value). But to fix the Lambda_a[10,1]
mixing, Task K's full bundle should also:
- Apply sqrt transform to per90_passes (stage-23 candidate; skew 1–2 tier),
  which may reshape the distribution's relationship with per90_goals
- Use longer warmup (≥ 1000 iterations, delta=0.95) in the production refit
- If still poor after full bundle: the LT constraint places per90_passes in a
  "pivot" role for factor 1 (column 10 in the LT ordering), which may create
  an intrinsically curved posterior ridge. Document as a known convergence
  limitation.

_Written 2026-07-01, stage 26 (ablation fit, ~4684s wall clock)_
