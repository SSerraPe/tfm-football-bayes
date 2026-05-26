# Comprehensive Validation Analysis

Generated: 2026-05-14

This note summarizes the complete validation state of the additive Bayesian
modelling pipeline after fitting and diagnosing:

- `03_real_diagonal_additive`
- `05_sim_diagonal_additive`
- `06_sim_lowrank_recovery`

It explains what was fitted, what was validated, what the results mean, which
results are reliable, which are not, and what should be tried next. No model
code or fitted results are changed by this note.

## 0. Thesis-Level Theoretical Framing

This analysis is intended to support a master's thesis, so the statistical
claims must be stated with more care than an ordinary exploratory notebook.
The main goal is not simply to fit a model that runs, but to establish which
football-performance quantities can be interpreted as stable player-level
variation, which can be interpreted as season-level variation, and which remain
mostly unexplained row-level variation after accounting for players and seasons.

### 0.0 Notation And Definitions

This section defines the notation, acronyms, diagnostics, and modelling terms
used throughout the note. Later sections use these meanings.

#### Data Units

`Player-season row`

: One observation corresponding to one player in one season, after aggregating
  lower-level Wyscout match-player records.

`Observed player-season`

: A player-season combination that exists in the data after filtering and
  aggregation. The pipeline does not create rows for seasons in which a player
  is absent.

`Panel`

: A dataset where units are observed repeatedly over time. Here the units are
  players and the time periods are seasons.

`Irregular panel`

: A panel where not every unit is observed in every time period. Here, some
  players appear once, some appear repeatedly, some disappear, and some reappear
  after gaps.

`Feature`

: One measured football variable used as a response column, such as
  `per90_shots`, `per90_passes`, or `rate_successful_crosses`.

`per90 feature`

: A count-like football statistic converted to a rate per 90 minutes. For
  example, `per90_shots` is shots per 90 minutes played.

`rate feature`

: A ratio or success-rate feature, usually successful actions divided by
  attempted actions. For example, `rate_successful_crosses` is successful
  crosses divided by total crosses.

`Scaled feature`

: A feature transformed by subtracting its empirical mean and dividing by its
  empirical standard deviation. Scaling makes features comparable on a common
  unitless scale.

`Response`

: The modelled outcome. Here the response is multivariate: each player-season
  row has a vector of 53 scaled football features.

#### Indices And Dimensions

`n`

: Row index. It identifies one observed player-season row.

`p`

: Feature index. It identifies one selected football feature.

`i`

: Player index. It identifies one player.

`j`

: Season index. It identifies one season.

`N`

: Total number of observed player-season rows. In this analysis, `N = 4944`.

`P`

: Total number of selected modelled features. In this analysis, `P = 53`.

`I`

: Total number of players. In this analysis, `I = 1650`.

`S`

: Total number of seasons. In this analysis, `S = 12`.

`i[n]`

: The player index associated with row `n`.

`j[n]`

: The season index associated with row `n`.

`y_n`

: The full response vector for row `n`; it contains all `P` scaled features for
  that player-season.

`y_np`

: The value of feature `p` in row `n`.

#### Additive Model Terms

`mu`

: The population-level mean vector in the scaled feature space.

`a_i`

: The player-effect vector for player `i`, representing that player's stable
  deviation from the population mean within the observed data.

`b_j`

: The season-effect vector for season `j`, representing that season's common
  deviation from the population mean across observed rows.

`epsilon_n`

: The residual vector for row `n`, meaning the remaining player-season
  variation not explained by `mu`, `a_i`, or `b_j`.

`Random effect`

: A group-level deviation estimated with partial pooling. Player and season
  effects are random effects in this pipeline.

`Partial pooling`

: Bayesian sharing of information across units. It prevents every player or
  season effect from being estimated as if it were completely unrelated to all
  others.

`Variance component`

: The amount of model-implied variance attributed to one source, such as player,
  season, or residual variation.

`Player variance`

: Feature-specific variance attributed to stable player-level differences. In
  the diagonal model it is `sigma_a,p^2`.

`Season variance`

: Feature-specific variance attributed to common season-level shifts. In the
  diagonal model it is `sigma_b,p^2`.

`Residual variance`

: Feature-specific variance left after player and season effects. In the
  diagonal model it is `sigma_e,p^2`.

`Variance share`

: A variance component divided by total model-implied variance for the same
  feature.

#### Distribution And Matrix Terms

`Normal distribution`

: A probability distribution described by a mean and variance. The notation
  `Normal(0, sigma^2)` means mean zero and variance `sigma^2`.

`Multivariate`

: Involving multiple variables at the same time. Here the response has 53
  variables.

`Covariance`

: A measure of how two variables vary together.

`Covariance matrix`

: A square matrix with variances on the diagonal and covariances off the
  diagonal.

`Diagonal covariance`

: A covariance matrix with zero off-diagonal covariances. It gives each feature
  its own variance but does not model cross-feature covariance inside that
  component.

`diag(x)`

: A diagonal matrix with vector `x` on its diagonal, or the diagonal entries of
  a matrix, depending on context.

`sigma`

: A standard deviation or scale parameter. Variance is `sigma^2`.

`sigma_a,p`

: Player-effect standard deviation for feature `p`.

`sigma_b,p`

: Season-effect standard deviation for feature `p`.

`sigma_e,p`

: Residual standard deviation for feature `p`.

`Sigma_a`

: Player-effect covariance matrix.

`Sigma_b`

: Season-effect covariance matrix.

`Sigma_epsilon`

: Residual covariance matrix.

`Full rank`

: A matrix property meaning that no selected column is an exact linear
  combination of the others.

`Correlation`

: A standardized linear association between two quantities, ranging from `-1`
  to `1`.

#### Low-Rank And Factor Terms

`Low-rank covariance`

: A covariance structure represented by a small number of latent dimensions
  rather than by an unrestricted full covariance matrix.

`Rank`

: The number of latent dimensions in the low-rank component.

`Latent factor`

: An unobserved dimension used to explain shared variation across observed
  features.

`Loading`

: A coefficient linking a latent factor to an observed feature.

`Loading matrix`

: A matrix collecting all factor loadings. Rows are features and columns are
  latent factors.

`Lambda_a`

: Player-effect loading matrix in the low-rank simulation and recovery model.

`Lambda_b`

: Season-effect loading matrix in the low-rank simulation and recovery model.

`Psi_a`

: Player-effect uniqueness variances, meaning feature-specific player variance
  not explained by player latent factors.

`Psi_b`

: Season-effect uniqueness variances, meaning feature-specific season variance
  not explained by season latent factors.

`Uniqueness`

: Feature-specific variance not explained by common latent factors.

`Anchor variable`

: A feature selected to identify a factor.

`Pure-anchor constraint`

: A factor-identification constraint where each anchor loads positively on its
  assigned factor and has zero loading on other factors.

`Identifiability`

: The ability to learn a parameter uniquely from the model and data. If many
  parameter values imply nearly the same likelihood, the parameter is weakly
  identified.

`Rotation ambiguity`

: A factor-model problem where different rotated loading matrices imply the
  same covariance.

`Sign ambiguity`

: A factor-model problem where multiplying a factor and its loadings by `-1`
  leaves the implied covariance unchanged.

`Loading-uniqueness tradeoff`

: A factor-model problem where the model can explain variance through common
  loadings or through uniqueness terms, creating weakly identified posterior
  regions.

#### Bayesian Computation Terms

`Prior`

: A probability distribution expressing assumptions about a parameter before
  seeing the data.

`Likelihood`

: The probability model for observed data given parameter values.

`Posterior`

: The probability distribution of parameters after combining prior information
  with the observed data.

`Posterior mean`

: The average of a parameter over posterior draws.

`Posterior median`

: The middle posterior value of a parameter.

`Credible interval`

: An interval summarizing posterior uncertainty. In this note, `q5` to `q95`
  is the central 90% posterior interval.

`q5`

: The 5th percentile of posterior draws.

`q95`

: The 95th percentile of posterior draws.

`Draw`

: One sampled value from the posterior distribution.

`Chain`

: One independent Markov chain used to sample from the posterior.

`Warmup`

: Early Stan iterations used to tune the sampler. Warmup draws are not used for
  posterior summaries here.

`Sampling iterations`

: Post-warmup iterations used for posterior inference.

`HMC`

: Hamiltonian Monte Carlo, the Markov chain Monte Carlo algorithm family used
  by Stan.

`NUTS`

: No-U-Turn Sampler, Stan's adaptive HMC algorithm.

`CmdStan`

: The command-line interface to Stan used to compile and sample the models.

`CmdStanR`

: The R interface used to run CmdStan from R scripts.

`CSV`

: Comma-separated values file. CmdStan writes posterior draws and diagnostics
  to CSV files.

`RDS`

: R's serialized single-object file format. This project uses lightweight RDS
  fit pointers rather than huge full fit objects.

#### Diagnostic Terms

`Rhat`

: A convergence diagnostic comparing variation between chains and within
  chains. Values close to 1 are good. Values above 1.01 indicate possible
  convergence issues.

`ESS`

: Effective sample size. It estimates how many independent posterior draws the
  autocorrelated Markov chain draws are worth.

`Bulk ESS`

: Effective sample size for the central part of the posterior distribution.

`Tail ESS`

: Effective sample size for the tails of the posterior distribution.

`Divergence`

: A Stan warning that HMC had difficulty following the posterior geometry.
  Divergences can indicate biased sampling.

`Treedepth`

: A measure of how long NUTS had to simulate a trajectory before proposing a
  draw.

`E-BFMI`

: Energy Bayesian Fraction of Missing Information. It checks whether HMC is
  exploring the energy distribution adequately.

`accept_stat`

: Stan's proposal acceptance statistic.

`adapt_delta`

: Stan tuning parameter controlling target acceptance probability.

`max_treedepth`

: Maximum allowed NUTS tree depth.

#### Recovery And Validation Metrics

`Bias`

: Average signed error: posterior mean minus true value.

`MAE`

: Mean absolute error.

`RMSE`

: Root mean squared error.

`Coverage`

: Proportion of true values falling inside a posterior interval.

`Truth`

: The known parameter value used to generate simulated data. Truth is available
  for simulations only, not for real data.

`Recovery`

: The degree to which posterior summaries match known simulation truth.

#### Predictive Validation Terms

`log_lik`

: Pointwise log-likelihood. It records how much probability the model assigns
  to each observed row under each posterior draw.

`LOO`

: Leave-one-out cross-validation.

`PSIS-LOO`

: Pareto-smoothed importance-sampling leave-one-out cross-validation, an
  approximation to exact LOO that avoids refitting the model once per
  observation.

`elpd_loo`

: Expected log predictive density estimated by LOO. Higher is better when
  comparing models on the same data and when diagnostics are reliable.

`p_loo`

: Effective number of parameters estimated by LOO.

`looic`

: LOO information criterion, equal to `-2 * elpd_loo`. Lower is better when the
  LOO approximation is reliable.

`Pareto k`

: Diagnostic for PSIS approximation stability. Values above 0.7 are
  problematic; values above 1 indicate severe unreliability.

`K-fold cross-validation`

: Predictive validation where the data are split into `K` parts and the model
  is refit `K` times, each time holding out one part.

The theoretical object of interest is a multivariate player-season response:

```text
y_n = (y_n1, ..., y_nP)
```

where each row `n` corresponds to one observed player-season and each column is
one standardized football feature. In this pipeline:

```text
N = 4944 observed player-seasons
P = 53 selected scaled features
I = 1650 players
S = 12 seasons
```

The index functions are:

```text
i[n] = player index for row n
j[n] = season index for row n
```

The central additive decomposition is:

```text
y_n = mu + a_i[n] + b_j[n] + epsilon_n
```

where:

- `mu` is the population-level mean vector in the scaled feature space.
- `a_i` is the player-specific deviation vector for player `i`.
- `b_j` is the season-specific deviation vector for season `j`.
- `epsilon_n` is the remaining player-season residual vector.

Because all features are standardized before modelling, the variance components
are interpretable on a common scale. A variance value of 1 is roughly one
standardized feature variance unit. This is essential for comparing football
features that originally have different units, such as shots per 90, pass
success rates, recoveries, and card counts.

### 0.1 Estimands

For the real-data diagonal additive model, the thesis-level estimands are:

```text
Var_player,p   = sigma_a,p^2
Var_season,p   = sigma_b,p^2
Var_residual,p = sigma_e,p^2
```

for each feature `p`. The proportional variance decomposition is:

```text
Prop_player,p   = Var_player,p / (Var_player,p + Var_season,p + Var_residual,p)
Prop_season,p   = Var_season,p / (Var_player,p + Var_season,p + Var_residual,p)
Prop_residual,p = Var_residual,p / (Var_player,p + Var_season,p + Var_residual,p)
```

These proportions are the most interpretable thesis outputs. They answer:

- How much of a feature's variation is stable across players?
- How much is shared by seasons?
- How much remains row-specific after accounting for player and season?

Importantly, these are model-based variance components, not causal effects.
For example, a large player component for `per90_shots` means that players
systematically differ in shot volume after pooling across observed seasons. It
does not prove that player identity causally produces shot volume independently
of position, team, tactical role, or selection effects.

### 0.2 What The Model Conditions On

The model conditions on the observed panel structure. It does not model:

- why a player appears in a season,
- why a player is absent in another season,
- transfers,
- injuries,
- team identity,
- position,
- coach or tactical system,
- match-level opponent strength,
- minutes selection beyond the minimum-minutes filtering.

Therefore the player effect should be interpreted as a stable player-season
profile effect within the observed data-generating and selection process, not
as a pure intrinsic ability parameter.

The season effect should be interpreted as a common seasonal shift shared
across observed player-season rows. It can capture league-wide tactical,
measurement, event-definition, or temporal changes. It should not be
overinterpreted as a causal calendar-year effect unless additional design
controls are added.

### 0.3 Why Missing Player-Seasons Are Not Imputed

The panel is irregular: players appear, disappear, and reappear. The modelling
unit is an observed player-season row. A missing player-season is not treated
as an unobserved response to impute. This is theoretically important.

If missing player-seasons were imputed without modelling the appearance process,
the analysis would implicitly assume that the absence mechanism is ignorable.
In football data, that is unlikely: absences can be caused by transfers,
injuries, relegation, age, playing time, retirement, or selection. The current
approach avoids making a strong and probably false imputation assumption.

The cost is that conclusions apply to observed eligible player-seasons, not to
all possible player-season combinations.

### 0.4 Diagonal Covariance As A Deliberate Thesis Model

The real-data model assumes diagonal covariance for player, season, and residual
effects:

```text
Cov(a_i)       = diag(sigma_a^2)
Cov(b_j)       = diag(sigma_b^2)
Cov(epsilon_n) = diag(sigma_e^2)
```

This does not mean the football features are truly independent. Instead, the
diagonal model is a deliberate first-order decomposition that estimates one
variance component per feature and per level. Its strength is interpretability:
each feature receives a clear player, season, and residual variance share.

The diagonal model is appropriate for the main thesis claim if the thesis
question is:

```text
Which features are mostly player-stable, season-shifted, or residual/noisy?
```

It is not sufficient if the thesis question is:

```text
What is the latent covariance structure connecting different features?
```

That second question requires a factor or covariance model and stronger
identification work.

### 0.5 Low-Rank Model As A Verification Device

The low-rank simulation and recovery model serve a different theoretical role.
They are not the main real-data interpretive model. They are a validation
device. The simulation generates data from:

```text
Cov(a_i) = Lambda_a Lambda_a' + diag(Psi_a)
Cov(b_j) = Lambda_b Lambda_b' + diag(Psi_b)
```

with known truth. This allows two types of verification:

1. Misspecification robustness:

   Fit the simpler diagonal model to low-rank data and check whether it still
   recovers covariance diagonals. This asks whether the main variance
   decomposition remains useful when true covariance is not diagonal.

2. Recovery:

   Fit a matching low-rank model and check whether known loadings, uniquenesses,
   and covariance diagonals are recovered.

The results support the first point strongly and the second point only
partially.

### 0.6 Identifiability Of Loadings

Factor loading matrices are not automatically identifiable. Without constraints,
many loading matrices can imply the same covariance because of rotations,
reflections, and scaling tradeoffs. The recovery model uses pure-anchor
constraints:

- selected anchor variables have positive own-factor loadings,
- off-factor anchor loadings are fixed to zero,
- non-anchor loadings are free.

This removes some rotational and sign ambiguity, but not all practical
identifiability problems. The posterior can still trade off:

```text
loading magnitude vs latent score magnitude
loading variance vs uniqueness variance
one loading pattern vs another weakly identified pattern
```

The poor Rhat values in the low-rank recovery model show that the current
constraints are not strong enough for stable loading-level interpretation.
This matters for the thesis: the low-rank loading results should be framed as
a failed or incomplete recovery experiment, not as substantive football-factor
findings.

### 0.7 Difference Between Covariance Recovery And Loading Recovery

Recovering a covariance diagonal is easier than recovering a loading matrix.
A covariance diagonal is invariant to many transformations that change loadings.
For example, two different loading configurations can imply similar total
feature variance. Therefore:

- good recovery of `Sigma_a_diag` does not imply good recovery of `Lambda_a`;
- good recovery of residual scales does not imply identifiable latent factors;
- poor loading recovery does not necessarily invalidate diagonal variance
  decomposition.

This distinction is central for the thesis. The current pipeline supports
variance-decomposition claims more strongly than latent-factor interpretation.

### 0.8 Posterior Diagnostics As Part Of The Scientific Argument

The diagnostic quantities are not merely software checks. They determine which
claims are scientifically defensible.

Rhat measures whether chains agree. Values near 1 indicate that chains are
sampling similar posterior regions. Large values indicate non-convergence or
multimodality.

Effective sample size measures how much independent posterior information the
draws contain after autocorrelation. Low ESS means posterior means and intervals
are less stable.

Divergences and treedepth hits diagnose HMC geometry. In this project, the
geometry is clean, so the main concern is not divergent exploration but
identifiability and mixing.

For thesis reporting:

- parameters with clean Rhat and adequate ESS can be interpreted normally;
- parameters with mild Rhat/ESS issues can be interpreted cautiously;
- parameters with severe Rhat/ESS issues should not be used for substantive
  conclusions.

### 0.9 PSIS-LOO As Predictive Validation

PSIS-LOO estimates expected out-of-sample predictive accuracy by reweighting
posterior draws to approximate leave-one-observation-out refits. It is useful
only when the importance weights are stable. Pareto-k diagnostics assess that
stability.

In this project, PSIS-LOO was computed, but Pareto-k diagnostics failed badly.
Therefore the LOO table should be reported as a diagnostic exercise, not as
decisive model-selection evidence.

The theoretical reason is that hierarchical models with player effects can
make individual rows highly influential. Leaving out a row can substantially
change the relevant player effect, especially for one-season players. Standard
row-wise PSIS-LOO is often fragile in this setting.

For a master's thesis, the safest statement is:

```text
PSIS-LOO was attempted and computed, but Pareto-k diagnostics indicate that
the approximation is unreliable. Exact K-fold cross-validation is recommended
for predictive model comparison.
```

### 0.10 Hierarchy Of Evidence In This Analysis

The results should be weighted in the following order:

1. Sampler geometry:

   All three models have clean geometry. This supports the claim that the
   samplers ran without major HMC pathologies.

2. Posterior convergence:

   The real diagonal and simulated diagonal models are acceptable. The low-rank
   recovery model is not acceptable for loading-level interpretation.

3. Simulation recovery:

   Diagonal variance recovery is strong. Low-rank loading recovery is weak.

4. Real-data interpretation:

   Real-data variance decomposition is defensible for player/season/residual
   shares, with mild caveats.

5. PSIS-LOO:

   Computed but unreliable because Pareto-k diagnostics fail.

This hierarchy should guide how the thesis discusses results.

## 1. Data And Modelling Setup

The modelling unit is an observed player-season row. Missing player-seasons are
not imputed. The pipeline uses scaled feature values, player indices, and season
indices built from the same observed Wyscout-derived panel.

The model object summary is:

| quantity | value |
|---|---:|
| observed player-season rows | 4944 |
| selected features | 53 |
| players | 1650 |
| seasons | 12 |
| one-season players | 605 |
| repeated players | 1045 |
| players with gaps | 346 |
| feature matrix rank | 53 |
| max absolute feature correlation | 0.972 |

The feature matrix is full rank. The largest retained correlation is high
(`0.972`), between closely related pass-volume variables. This does not break
the diagonal additive models, but it is a warning that some features carry
nearly redundant information.

The simulation uses the same `N`, `P`, players, seasons, and row mapping:

| quantity | value |
|---|---:|
| simulation seed | 20260913 |
| rows | 4944 |
| features | 53 |
| players | 1650 |
| seasons | 12 |
| player low-rank dimension | 2 |
| season low-rank dimension | 1 |

The simulation was intentionally generated from known low-rank plus diagonal
player and season covariance components. This lets us check two things:

1. Whether a simpler diagonal model can recover diagonal variance components
   even under low-rank misspecification.
2. Whether the low-rank recovery model can recover known loadings and
   uniquenesses.

## 2. Models Fitted

### 2.1 Real Diagonal Additive Model

The real-data model is:

```text
y_n = mu + a_player[n] + b_season[n] + epsilon_n
```

with diagonal covariance in the scaled feature space:

```text
a_i       ~ Normal(0, diag(sigma_a^2))
b_j       ~ Normal(0, diag(sigma_b^2))
epsilon_n ~ Normal(0, diag(sigma_e^2))
```

The player and season effects are non-centered and column-centered so that
`mu` remains identifiable.

This model is intended as a descriptive additive variance decomposition:

- player-level stable differences
- season-level common shifts
- residual row-level variation

### 2.2 Simulated Diagonal Additive Model

This is the same diagonal additive model fitted to simulated data. The
simulation truth is low-rank plus diagonal, so this model is knowingly
misspecified. The goal is not to recover loadings. The goal is to check whether
the pipeline can recover diagonal covariance components reasonably well even
when the true covariance is not diagonal.

### 2.3 Simulated Low-Rank Recovery Model

The low-rank recovery model matches the simulation more closely:

```text
Sigma_a = Lambda_a Lambda_a' + diag(Psi_a)
Sigma_b = Lambda_b Lambda_b' + diag(Psi_b)
Sigma_e = diag(sigma_epsilon^2)
```

The model estimates:

- player loading matrix `Lambda_a`
- season loading matrix `Lambda_b`
- player uniqueness scales
- season uniqueness scales
- residual scales
- covariance diagonals implied by the above

It uses pure-anchor loading constraints: the first `Q` variables are positive
own-factor anchors and have zero off-factor loadings. This is meant to remove
rotation and sign ambiguity enough for direct comparison with known simulation
truth.

## 3. Computational Status

All three models completed and wrote durable fit artifacts:

- lightweight fit pointers in `model/fits/`
- persistent CmdStan CSV files in `model/fits/csv/`
- posterior summaries in `model/outputs/tables/`
- sampler diagnostics in `model/outputs/diagnostics/`

Relevant output files include:

- `outputs/tables/03_real_diagonal_additive_posterior_summary.csv`
- `outputs/tables/05_sim_diagonal_additive_posterior_summary.csv`
- `outputs/tables/06_sim_lowrank_recovery_posterior_summary.csv`
- `outputs/diagnostics/08_posterior_health_summary.csv`
- `outputs/diagnostics/08_sampler_health_summary.csv`
- `outputs/tables/08_loo_summary.csv`
- `outputs/tables/08_sim_diagonal_recovery_metrics.csv`
- `outputs/tables/08_sim_lowrank_loading_recovery_metrics.csv`
- `outputs/tables/08_sim_lowrank_variance_recovery_metrics.csv`

The earlier memory failures were not sampling failures. They happened during
post-processing when R tried to load or summarize very large draw objects. The
current validation path avoids full draw loading except where specifically
requested for LOO.

## 4. Sampler Geometry

All three models have clean HMC geometry:

| model | draws | divergences | max treedepth | treedepth hits | min E-BFMI | mean accept |
|---|---:|---:|---:|---:|---:|---:|
| real diagonal | 4000 | 0 | 8 | 0 | 0.561 | 0.953 |
| simulated diagonal | 4000 | 0 | 8 | 0 | 0.551 | 0.951 |
| simulated low-rank | 4000 | 0 | 8 | 0 | 0.586 | 0.951 |

Interpretation:

- There are no divergent transitions.
- There are no treedepth saturations at the configured limit of 12.
- E-BFMI values are above common warning thresholds.
- Average acceptance rates are close to the target `adapt_delta = 0.95`.

Therefore the main problems are not HMC geometry problems. The main problems
are posterior mixing, identifiability, and predictive validation diagnostics.

## 5. Posterior Health

Posterior summary diagnostics:

| model | summarized parameters | max Rhat | Rhat > 1.01 | Rhat > 1.05 | min bulk ESS | min tail ESS |
|---|---:|---:|---:|---:|---:|---:|
| real diagonal | 531 | 1.035 | 13 | 0 | 159 | 175 |
| simulated diagonal | 531 | 1.010 | 3 | 0 | 779 | 1412 |
| simulated low-rank | 688 | 1.739 | 233 | 208 | 6 | 107 |

### 5.1 Real Diagonal Model

The real diagonal model is mostly usable. It has mild Rhat and ESS warnings,
mainly in a few player variance parameters. The worst parameters are:

| parameter | feature | posterior mean | sd | Rhat | bulk ESS | tail ESS |
|---|---|---:|---:|---:|---:|---:|
| `var_a[11]` | `per90_forward_passes` | 0.827 | 0.032 | 1.035 | 159 | 290 |
| `sigma_a[11]` | `per90_forward_passes` | 0.909 | 0.018 | 1.035 | 159 | 290 |
| `prop_e[11]` | `per90_forward_passes` | 0.150 | 0.006 | 1.028 | 202 | 642 |
| `prop_a[11]` | `per90_forward_passes` | 0.846 | 0.006 | 1.027 | 179 | 376 |
| `var_a[44]` | `per90_sliding_tackles` | 0.652 | 0.029 | 1.016 | 894 | 1643 |
| `var_a[16]` | `per90_progressive_passes` | 0.752 | 0.029 | 1.015 | 509 | 1746 |
| `var_a[24]` | `rate_successful_crosses` | 0.022 | 0.013 | 1.012 | 197 | 177 |

These are not catastrophic. No Rhat exceeds 1.05. But for a thesis, these
parameters should be described as mildly under-mixed rather than perfectly
converged.

### 5.2 Simulated Diagonal Model

The simulated diagonal model is healthy. The maximum Rhat is 1.010 and no ESS
values are below 400. This is strong evidence that the diagonal model and
pipeline can sample well on data with the same structure as the observed panel.

### 5.3 Simulated Low-Rank Recovery Model

The low-rank recovery model is not healthy for loading interpretation. It has:

- max Rhat 1.739
- 233 parameters above Rhat 1.01
- 208 parameters above Rhat 1.05
- minimum bulk ESS around 6

The worst parameters are low-rank loadings and uniqueness terms:

| parameter | posterior mean | sd | Rhat | bulk ESS |
|---|---:|---:|---:|---:|
| `lambda_b_free[37]` / `Lambda_b[38,1]` | -0.010 | 0.271 | 1.739 | 6.1 |
| `lambda_a_free[41]` / `Lambda_a[43,1]` | 0.000 | 0.111 | 1.738 | 6.2 |
| `lambda_a_free[23]` / `Lambda_a[25,1]` | 0.000 | 0.210 | 1.738 | 6.1 |
| `lambda_b_free[24]` / `Lambda_b[25,1]` | 0.018 | 0.328 | 1.738 | 6.1 |
| `psi_b[1]` | 0.332 | 0.246 | 1.737 | 6.2 |
| `psi_a[1]` | 0.418 | 0.127 | 1.737 | 6.2 |

This is a strong identifiability and mixing warning. The model has clean HMC
trajectories but chains are not exploring the same posterior region for many
loading parameters.

## 6. Real-Data Variance Decomposition

The real diagonal additive model decomposes scaled feature variance into player,
season, and residual components. Averaged across features:

| component | mean share |
|---|---:|
| player | 0.634 |
| season | 0.047 |
| residual | 0.319 |

This means the dominant source of variation in most selected features is stable
player-level heterogeneity. Season-level variation exists but is much smaller
on average. Residual variation remains important, especially for rate features.

### 6.1 Features Dominated By Player Effects

Highest player-level variance shares:

| feature | player | season | residual |
|---|---:|---:|---:|
| `per90_recoveries` | 0.884 | 0.016 | 0.100 |
| `per90_touch_in_box` | 0.872 | 0.001 | 0.127 |
| `per90_offensive_duels` | 0.869 | 0.056 | 0.075 |
| `per90_shots` | 0.860 | 0.008 | 0.132 |
| `per90_forward_passes` | 0.846 | 0.004 | 0.150 |
| `per90_crosses` | 0.846 | 0.006 | 0.148 |
| `per90_long_passes` | 0.843 | 0.005 | 0.152 |
| `per90_clearances` | 0.837 | 0.034 | 0.129 |
| `per90_passes` | 0.829 | 0.011 | 0.160 |
| `per90_back_passes` | 0.823 | 0.008 | 0.169 |

Interpretation: volume and style variables are strongly player-specific. This
is plausible because players have stable roles, positions, and tactical usage.

### 6.2 Features With The Largest Season Effects

Highest season-level variance shares:

| feature | player | season | residual |
|---|---:|---:|---:|
| `per90_smart_passes` | 0.440 | 0.376 | 0.183 |
| `per90_pressing_duels` | 0.455 | 0.320 | 0.225 |
| `per90_accelerations` | 0.598 | 0.164 | 0.238 |
| `per90_dribbles_against` | 0.693 | 0.136 | 0.171 |
| `rate_defensive_duels_won` | 0.244 | 0.129 | 0.627 |
| `per90_through_passes` | 0.583 | 0.121 | 0.296 |
| `per90_sliding_tackles` | 0.606 | 0.119 | 0.275 |
| `per90_loose_ball_duels` | 0.718 | 0.113 | 0.168 |
| `per90_dribbles` | 0.800 | 0.101 | 0.100 |
| `rate_successful_progressive_passes` | 0.484 | 0.084 | 0.432 |

Interpretation: the strongest season effects occur in features that plausibly
reflect changing tactical environments, league-wide style, event definitions,
or tactical fashion over seasons. The most prominent examples are smart passes
and pressing duels.

### 6.3 Features Dominated By Residual Effects

Highest residual variance shares:

| feature | player | season | residual |
|---|---:|---:|---:|
| `rate_successful_crosses` | 0.022 | 0.003 | 0.975 |
| `rate_successful_sliding_tackles` | 0.057 | 0.011 | 0.932 |
| `rate_shots_on_target` | 0.167 | 0.001 | 0.832 |
| `rate_goal_conversion` | 0.180 | 0.003 | 0.818 |
| `rate_dribbles_against_won` | 0.184 | 0.008 | 0.808 |
| `rate_successful_dribbles` | 0.231 | 0.009 | 0.760 |
| `rate_defensive_duels_won` | 0.244 | 0.129 | 0.627 |
| `per90_assists` | 0.391 | 0.010 | 0.599 |
| `rate_offensive_duels_won` | 0.412 | 0.016 | 0.573 |
| `per90_yellow_cards` | 0.441 | 0.017 | 0.542 |

Interpretation: rate variables and sparse outcome variables are noisy. They
depend heavily on opportunity counts and row-level context. They are much less
stable than volume/style features.

## 7. Simulation Validation: Diagonal Model On Low-Rank Data

The simulated diagonal model is intentionally misspecified because the data
were generated from low-rank plus diagonal covariance. Despite this, it recovers
diagonal variance components well.

Recovery metrics:

| component | bias | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|---:|
| player | 0.0028 | 0.0078 | 0.0112 | 0.996 | 0.906 |
| residual | -0.0009 | 0.0040 | 0.0055 | 0.998 | 0.925 |
| season | 0.0250 | 0.0300 | 0.0468 | 0.885 | 0.868 |

This validates the diagonal additive machinery for variance decomposition:

- player variance recovery is excellent
- residual variance recovery is excellent
- season variance recovery is weaker but still useful

The season component is hardest because there are only 12 seasons. This is a
small number of season-level units for estimating feature-by-feature season
variance.

Worst season variance errors in the simulated diagonal model:

| feature | posterior mean | truth | error |
|---|---:|---:|---:|
| `rate_successful_sliding_tackles` | 0.354 | 0.139 | 0.214 |
| `per90_xg_shot` | 0.221 | 0.103 | 0.118 |
| `per90_goals` | 0.236 | 0.151 | 0.085 |
| `per90_fouls` | 0.144 | 0.064 | 0.080 |
| `per90_defensive_duels` | 0.130 | 0.057 | 0.073 |

The bias is mostly upward for season variance. This is consistent with a model
trying to approximate low-rank seasonal covariance using diagonal season
variance terms.

## 8. Simulation Validation: Low-Rank Recovery

The low-rank model was intended to recover the true low-rank simulation
parameters. It partially succeeds for some covariance summaries but fails as a
direct loading recovery model under the current specification.

### 8.1 Loading Recovery

Loading recovery metrics:

| parameter | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|
| `Lambda_a` | 0.082 | 0.131 | 0.804 | 0.840 |
| `Lambda_b` | 0.102 | 0.124 | 0.262 | 1.000 |

The player loading matrix has moderate correlation with truth but still has
large errors and poor mixing. The season loading matrix is not recovered well:
the correlation with truth is only 0.262, and coverage is high because
posterior intervals are too wide, not because estimates are precise.

Worst `Lambda_a` errors:

| feature | factor | posterior mean | truth | error | Rhat |
|---|---:|---:|---:|---:|---:|
| `per90_shot_assists` | 1 | -0.003 | 0.491 | -0.494 | 1.736 |
| `per90_fouls` | 1 | -0.001 | 0.352 | -0.352 | 1.736 |
| `per90_shots_blocked` | 1 | -0.001 | 0.340 | -0.342 | 1.735 |
| `per90_opponent_half_recoveries` | 1 | 0.002 | -0.335 | 0.337 | 1.734 |
| `per90_long_passes` | 1 | 0.002 | -0.289 | 0.291 | 1.735 |

Worst `Lambda_b` errors:

| feature | factor | posterior mean | truth | error | Rhat |
|---|---:|---:|---:|---:|---:|
| `rate_successful_sliding_tackles` | 1 | -0.042 | 0.298 | -0.340 | 1.736 |
| `per90_xg_shot` | 1 | -0.008 | 0.230 | -0.238 | 1.737 |
| `per90_passes` | 1 | -0.003 | 0.230 | -0.233 | 1.733 |
| `per90_dangerous_own_half_losses` | 1 | -0.008 | 0.206 | -0.214 | 1.735 |
| `per90_through_passes` | 1 | 0.023 | -0.185 | 0.208 | 1.734 |

The pattern is clear: many free loadings have posterior means near zero,
wide intervals, and very high Rhat. This indicates chains are not agreeing on
the loading structure.

### 8.2 Variance And Scale Recovery

Variance recovery metrics:

| parameter | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|
| `Sigma_a_diag` | 0.009 | 0.013 | 0.997 | 0.811 |
| `Sigma_b_diag` | 0.057 | 0.093 | 0.886 | 0.585 |
| `psi_a_sd` | 0.011 | 0.022 | 0.941 | 0.925 |
| `psi_b_sd` | 0.040 | 0.052 | 0.740 | 0.925 |
| `sigma_epsilon_sd` | 0.004 | 0.005 | 0.998 | 0.925 |

The low-rank model does recover:

- residual scales very well
- player covariance diagonals very well
- player uniqueness scales reasonably well

It does not recover:

- season covariance diagonals precisely
- season uniqueness scales as well as player scales
- loading matrices reliably

Worst `Sigma_b_diag` errors:

| feature | posterior mean | truth | error | Rhat |
|---|---:|---:|---:|---:|
| `rate_successful_sliding_tackles` | 0.572 | 0.139 | 0.432 | 1.051 |
| `per90_goals` | 0.435 | 0.151 | 0.284 | 1.029 |
| `per90_xg_shot` | 0.342 | 0.103 | 0.238 | 1.009 |
| `per90_through_passes` | 0.219 | 0.084 | 0.134 | 1.013 |
| `per90_fouls` | 0.160 | 0.064 | 0.096 | 1.003 |

Season variance is again the weak point.

### 8.3 Why Low-Rank Loading Recovery Failed

The evidence points to identifiability and mixing problems, not HMC geometry.
Likely causes:

1. Only 12 season-level units

   A rank-1 season loading block is being inferred from very few season-level
   units. Even with many player-season rows, the number of independent season
   effects is only 12. That is weak information for a 53-dimensional loading
   vector plus uniqueness scales.

2. Loading and uniqueness tradeoff

   The model can explain covariance through loadings or through uniqueness
   terms. This creates posterior ridges where different combinations of
   `Lambda`, `psi`, and latent scores give similar likelihood values.

3. Anchors are positive but not strong enough

   The pure-anchor constraints impose signs and zero off-factor anchor loadings,
   but they do not force anchors to be strongly informative. In this run,
   `lambda_a_diag[1]` and `lambda_b_diag[1]` have very poor Rhat and low ESS.
   When anchor loading magnitude is weakly identified, the rest of the loading
   column is also unstable.

4. Factor-level multimodality remains

   Even with anchor constraints, chains can occupy different loading/uniqueness
   regimes. The high Rhat values around 1.7 and bulk ESS around 6 are consistent
   with chains stuck in different modes.

5. Direct loading recovery is stricter than covariance recovery

   A covariance matrix can be approximately recovered even when individual
   loadings are not. Loadings are not unique without strong identification,
   while covariance summaries are often more stable.

## 9. PSIS-LOO Results

LOO was computed with `BFA_RUN_LOO=true` after installing the missing `loo`
package. The output is:

- `outputs/tables/08_loo_summary.csv`

Results:

| model | elpd_loo | se_elpd_loo | p_loo | looic | Pareto k > 0.7 | Pareto k > 1 | max Pareto k |
|---|---:|---:|---:|---:|---:|---:|---:|
| real diagonal | -224673.4 | 1007.0 | 49545.7 | 449346.9 | 4354 | 3088 | 3.82 |
| simulated diagonal | -201763.1 | 394.7 | 41871.8 | 403526.1 | 4529 | 2878 | 2.88 |
| simulated low-rank | -196186.8 | 363.3 | 32963.7 | 392373.7 | 4215 | 2104 | 2.08 |

For the simulated data, the low-rank recovery model has better nominal
predictive fit than the diagonal simulated model:

```text
Delta elpd(low-rank - diagonal) = about 5576
```

This direction is expected because the data were generated from a low-rank
process.

However, the PSIS diagnostics are unacceptable:

- thousands of observations have Pareto k above 0.7
- thousands have Pareto k above 1
- max Pareto k is far above 1 for all models

Therefore the LOO estimates should not be used as formal model-selection
evidence. They are computed, but not trustworthy. They are better interpreted
as evidence that many rows are highly influential under conditional
leave-one-row-out scoring.

This is especially unsurprising in a hierarchical player-season model:

- many players appear only once
- player effects can be informed heavily by a single row
- leaving out that row changes the effective posterior substantially
- conditional log-likelihood with latent random effects can be too optimistic
  or unstable for LOO

For real predictive validation, exact K-fold cross-validation is preferable.

## 10. Main Scientific Conclusions

### 10.1 Real Data

The diagonal additive model gives a defensible descriptive result:

- player effects dominate many feature variances
- season effects are generally small but meaningful for some tactical/style
  variables
- residual variation dominates sparse rates and conversion/success variables

This is coherent with football interpretation:

- player role and style are stable
- league/tactical environment changes over seasons, but less than player
  identity
- rate variables are noisy because denominators are often smaller and context
  matters more

The model is acceptable for reporting a variance decomposition, with caveats
about mild convergence warnings.

### 10.2 Simulation

The simulated diagonal model validates the basic additive pipeline. It can
recover player and residual variance components very well, and season variance
moderately well.

The simulated low-rank model does not validate direct loading recovery. It
does validate that covariance diagonals and residual scales can be recovered
better than individual loading parameters.

### 10.3 Model Complexity

The current low-rank model is too ambitious for reliable loading recovery in
this setting, especially for season effects. The real data only have 12
seasons, and the simulation confirms that the season block is fragile even when
truth is known.

## 11. Problems Identified

### Problem 1: Mild Mixing Issues In The Real Diagonal Model

Symptoms:

- max Rhat 1.035
- 13 summarized parameters above 1.01
- minimum bulk ESS 159
- minimum tail ESS 175

This mostly affects player variance parameters such as
`per90_forward_passes`, `per90_sliding_tackles`, and
`per90_progressive_passes`.

Severity:

- mild to moderate
- does not invalidate the whole model
- should be disclosed

Possible fixes to try later:

1. Increase sampling iterations for the real diagonal model.
2. Run more chains or a second independent run to check stability.
3. Use slightly stronger priors on variance scales.
4. Revisit highly correlated features.
5. Consider grouping sparse or noisy rate variables separately.

### Problem 2: Low-Rank Loading Non-Identification

Symptoms:

- max Rhat 1.739
- many loadings have bulk ESS around 6
- many loading posterior means are near zero despite nonzero truth
- loading intervals are wide and cross zero

Severity:

- high
- loading interpretations from `06_sim_lowrank_recovery` are not reliable

Possible fixes to try later:

1. Use truth-informed initialization for simulation recovery experiments.

   This is not a general real-data solution, but it helps distinguish model
   non-identifiability from poor exploration.

2. Use stronger anchor constraints.

   Examples:

   - lower-bound anchor loadings away from zero
   - use stronger priors on anchor loadings
   - select anchor variables with strong expected signal
   - anchor on variables known to have large true simulated loadings in recovery
     tests

3. Separate covariance recovery from loading recovery.

   If the scientific target is covariance, report covariance summaries rather
   than raw loadings. Loadings are more fragile.

4. Add stronger shrinkage to free loadings.

   This may reduce multimodality and loading/uniqueness tradeoffs.

5. Use post-processing alignment for simulation diagnostics.

   Procrustes or sign alignment can help diagnose whether the model is
   recovering the subspace but not the exact anchored orientation. This is a
   diagnostic step, not a substitute for identification if direct loadings are
   the target.

6. Fit lower-rank or simpler season structures.

   With only 12 seasons, a low-rank season loading vector is hard to recover.
   A diagonal season effect or even no season factor may be more defensible.

7. Marginalize latent effects if feasible.

   Sampling many latent effects can create complex posterior geometry. A
   marginal covariance formulation may improve mixing, although crossed player
   and season effects make this computationally harder.

### Problem 3: Season Effects Are Hard To Estimate

Symptoms:

- season variance recovery is weaker than player and residual recovery
- low-rank season loading correlation is poor
- `Sigma_b_diag` has low coverage and upward bias

Cause:

- only 12 seasons
- many season-level covariance parameters
- low-rank season structure inferred from too few season units

Possible fixes to try later:

1. Keep season covariance diagonal in the main thesis model.
2. Use a very low-dimensional pre-specified season score instead of estimating
   a full loading vector.
3. Pool season variance more strongly across features.
4. Use informative priors for season effects.
5. Treat season effects as descriptive random intercept variance rather than
   attempting feature-level factor interpretation.

### Problem 4: PSIS-LOO Fails Pareto-k Diagnostics

Symptoms:

- real diagonal: 4354 rows with k > 0.7, 3088 with k > 1
- simulated diagonal: 4529 rows with k > 0.7, 2878 with k > 1
- simulated low-rank: 4215 rows with k > 0.7, 2104 with k > 1

Severity:

- high for model comparison
- LOO values are computed but not reliable

Likely causes:

1. Hierarchical random effects create influential observations.
2. One-season players produce nearly leave-one-group-out behavior when a row is
   removed.
3. The pointwise log-likelihood is conditional on latent effects, which can be
   problematic for hierarchical LOO.
4. Some player-season rows are high leverage in a 53-dimensional response.

Possible fixes to try later:

1. Use exact K-fold cross-validation.

   Suggested folds:

   - leave-player-out folds for generalizing to unseen players
   - leave-season-out folds for generalizing to new seasons
   - random row folds only if the prediction target is another observed row
     from the same player/season process

2. Use integrated or marginal log-likelihood for LOO if feasible.

   This avoids conditioning too strongly on row-informed random effects.

3. Use `loo` moment matching only as a secondary attempt.

   Because thousands of Pareto-k values are problematic, moment matching may be
   insufficient.

4. Report Pareto-k failures explicitly if LOO is mentioned.

### Problem 5: Memory And Output Size

Symptoms:

- full fit CSVs are very large
- full draw extraction can hit the 16 GB R vector limit
- trace plots from full fits are not practical without targeted extraction

Current mitigation:

- fit objects are lightweight pointers
- CmdStan CSVs are persistent
- posterior summaries are compact
- trace-plot loading was removed from stage `07`
- LOO is opt-in with `BFA_RUN_LOO=true`

Possible fixes to try later:

1. Add a dedicated targeted-draw extractor for a small variable set.
2. Store only required generated quantities for each workflow.
3. Split diagnostic jobs by model.
4. Add optional thinning for saved CSVs if trace-level detail is not needed.
5. Avoid writing `log_lik` unless predictive validation is planned.

## 12. Recommended Next Experiments

These are ideas only. They have not been implemented here.

### 12.1 Strengthen The Real Diagonal Result

Run a longer version of the real diagonal model:

```sh
BFA_RUN_STAN=true BFA_REUSE_FIT=false \
BFA_CHAINS=4 BFA_PARALLEL_CHAINS=4 \
BFA_ITER_WARMUP=1000 BFA_ITER_SAMPLING=2000 \
BFA_ADAPT_DELTA=0.95 BFA_MAX_TREEDEPTH=12 \
Rscript model/run_pipeline.R 03
```

Purpose:

- improve ESS for the few weak variance parameters
- check if Rhat drops below 1.01

This is not urgent unless the thesis needs precise intervals for the worst
parameters.

### 12.2 Build Exact K-Fold Validation

Add a K-fold script rather than relying on PSIS-LOO. Define the prediction
target first:

- unseen player: leave-player-out K-fold
- future season: leave-season-out K-fold
- held-out player-season row: random row K-fold

Given the structure of football panel data, leave-player-out or leave-season-out
is more interpretable than random row folds.

### 12.3 Low-Rank Recovery Debug Run With Truth Initialization

For simulation only, initialize `Lambda_a`, `Lambda_b`, `psi_a`, `psi_b`, and
`sigma_e` near the saved truth. If recovery is good with truth initialization
but poor from diffuse initialization, the issue is posterior exploration. If it
is still poor, the issue is deeper non-identifiability.

### 12.4 Stronger Anchor Experiment

Try a low-rank recovery model with:

- anchor loadings constrained above a positive floor
- stronger priors on anchor loadings
- stronger shrinkage on non-anchor loadings
- anchors chosen from high-signal simulated variables

This would directly test whether current low-rank failure is mostly due to weak
anchors.

### 12.5 Season Simplification

Given only 12 seasons, consider:

- diagonal season effects only
- no low-rank season block
- a single scalar season intensity score
- stronger pooled season prior

The simulation strongly suggests that feature-level season loadings are too
fragile under the current setup.

### 12.6 Feature Sensitivity

Run sensitivity checks with:

- no highly correlated pass-volume duplicate
- separate rate-only and volume-only feature sets
- position-stratified features if position metadata is available
- minimum-minute thresholds above 450

This can show whether variance decomposition conclusions are robust.

## 13. Recommended Thesis Framing

A defensible thesis framing is:

1. The real-data diagonal additive model is the main interpretable model.
2. It shows that stable player effects dominate most volume/style features,
   while rate and sparse outcome variables are much noisier.
3. Season effects exist but are feature-specific and generally smaller, with
   strongest signals in tactical/style variables such as smart passes and
   pressing duels.
4. Simulation confirms that diagonal additive variance recovery works well for
   player and residual components.
5. Low-rank recovery is not yet reliable for direct loading interpretation,
   especially for season effects.
6. PSIS-LOO was computed but fails Pareto-k diagnostics, so it should not be
   used as decisive model-selection evidence.
7. More complex covariance models require stronger identification, more careful
   validation, or exact cross-validation.

## 14. Bottom Line

The pipeline has succeeded for descriptive additive variance decomposition. The
real diagonal model is usable with mild convergence caveats. The simulation
confirms the diagonal model's ability to recover variance components. The
current low-rank recovery model should be treated as a diagnostic negative
result: it samples without divergences but does not identify loadings reliably.

The main next methodological improvements are:

- exact K-fold validation instead of PSIS-LOO
- stronger low-rank identifiability constraints
- simpler or more strongly pooled season structure
- targeted memory-safe diagnostics
- feature sensitivity analyses
