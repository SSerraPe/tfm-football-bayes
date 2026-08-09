# Stage 32 — Player-holdout K-selection (K = 0..8), with minutes-scaling sensitivity
#
# TWO ELPD metrics computed for each K from the same set of posterior draws:
#
#   Phase 1 (base): Normal marginal with constant sigma_e[p] per feature.
#     r_{i,t,p} ~ N(0, Sigma_a,pp × J_T + sigma_e_p^2 × I_T)
#     → uniform-diagonal Woodbury: O(T) per player.
#
#   Phase 2 (mv, phi=0.5): Same draws, minutes-scaled residual variance.
#     Cov_p = Sigma_a,pp × 11' + diag(sigma_e_p^2 × s_t^2)  where s_t = sqrt(m_ref/min_t)
#     → rank-1-plus-diagonal Woodbury: same O(T) cost, non-uniform d_t.
#
# Rationale: the base Phase 1 model is used for fast fitting (marginalized Normal,
# ~1k params). Phase 2 re-evaluates the holdout likelihood under minutes-scaled
# residuals using the same posterior draws as an approximation. If K* is the same
# under both metrics, the K ranking is robust. If it differs, Phase 2 governs
# (minutes scaling is the intended production model).
#
# Stan model: player_season_lowrank_marginalized.stan (Normal, no t-errors).
# K selection from Phase 2 (mv) is written to 32_player_holdout_elpd.csv and
# to 32_recommended_k.txt for consumption by downstream stages.

source("src/bootstrap.R")
check_packages(c(required_base_packages, "cmdstanr", "posterior"))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(ggplot2)
})

K_max      <- as.integer(Sys.getenv("BFA_RANK_MAX",   unset = "8"))
pf_paths   <- as.integer(Sys.getenv("BFA_PF_PATHS",  unset = "4"))
pf_draws   <- as.integer(Sys.getenv("BFA_PF_DRAWS",  unset = "600"))
pf_iters   <- as.integer(Sys.getenv("BFA_PF_ITERS",  unset = "500"))
seed_val <- pipeline$seed + 320L
HOLDOUT_FRAC <- 0.15

mo <- readRDS(file.path(paths$processed, "model_objects.rds"))
Y       <- mo$Y           # N × P
minutes <- mo$metadata$minutes
P <- ncol(Y); N <- nrow(Y); S <- length(mo$season_levels)
m_ref <- median(minutes)
message(sprintf("Loaded data: N=%d P=%d S=%d | m_ref=%.0f min", N, P, S, m_ref))

# ── 1. Holdout split (stratified by recurrence class) ─────────────────────────

recurrence    <- mo$player_recurrence
player_levels <- mo$player_levels
pi_all <- mo$player_index
si_all <- mo$season_index
set.seed(seed_val)

holdout_ids <- recurrence |>
  group_by(recurrence_class) |>
  slice_sample(prop = HOLDOUT_FRAC) |>
  ungroup() |>
  pull(player_id)

is_holdout_player <- player_levels %in% holdout_ids
holdout_row_mask  <- is_holdout_player[pi_all]
train_row_mask    <- !holdout_row_mask

N_train <- sum(train_row_mask); N_hold <- sum(holdout_row_mask)
message(sprintf("Holdout: %d players, %d obs | Training: %d obs",
                length(holdout_ids), N_hold, N_train))

# Reindex training players contiguously
train_player_ids <- player_levels[!is_holdout_player]
I_train <- length(train_player_ids)

Y_train  <- Y[train_row_mask, ]
pi_train <- match(player_levels[pi_all[train_row_mask]], train_player_ids)
si_train <- si_all[train_row_mask]

Y_hold       <- Y[holdout_row_mask, ]
si_hold      <- si_all[holdout_row_mask]
minutes_hold <- minutes[holdout_row_mask]

hold_player_ids <- player_levels[is_holdout_player]
I_hold <- length(hold_player_ids)
pi_hold <- match(player_levels[pi_all[holdout_row_mask]], hold_player_ids)

S_count_train <- as.integer(tabulate(pi_train, nbins = I_train))

# Minutes scale for holdout observations (phi = 0.5)
sigma_scale_hold <- sqrt(m_ref / minutes_hold)

write_csv(
  tibble(player_id = holdout_ids, split = "holdout"),
  file.path(paths$tables, "32_holdout_player_ids.csv")
)

# ── 2. PCA init for marginalized model ────────────────────────────────────────

make_marginalized_init <- function(stan_data) {
  sigma_floor <- stan_data$sigma_floor
  Yd <- stan_data$Y; pid <- stan_data$player_index
  sid <- stan_data$season_index; I <- stan_data$I
  Sv <- stan_data$S; Pd <- stan_data$P; K <- stan_data$Q_a

  sc <- tabulate(sid, nbins = Sv)
  B0 <- sweep(rowsum(Yd, sid), 1, pmax(sc, 1), "/")
  B0 <- sweep(B0, 2, colMeans(B0), "-")
  Y_adj <- Yd - B0[sid, ]
  pc <- tabulate(pid, nbins = I)
  A0 <- sweep(rowsum(Y_adj, pid), 1, pmax(pc, 1), "/")
  A0 <- sweep(A0, 2, colMeans(A0), "-")

  sigma_b0 <- pmax(apply(B0, 2, sd, na.rm = TRUE), sigma_floor)
  z_b0     <- sweep(B0, 2, sigma_b0, "/")
  E0       <- Yd - A0[pid, ] - B0[sid, ]
  sigma_e0 <- pmax(apply(E0, 2, sd, na.rm = TRUE), sigma_floor)

  if (K == 0L) {
    psi_a0 <- pmax(apply(A0, 2, sd, na.rm = TRUE), sigma_floor)
    return(function(chain_id = 1) list(
      psi_a   = psi_a0,
      z_b     = z_b0 + matrix(rnorm(Sv * Pd, 0, 0.05), Sv, Pd),
      sigma_b = sigma_b0, sigma_e = sigma_e0
    ))
  }

  S_A <- crossprod(A0) / I
  eig <- eigen(S_A, symmetric = TRUE)
  vals <- pmax(eig$values, 0); vecs <- eig$vectors
  sig2_iso  <- if (K < Pd) mean(vals[(K + 1):Pd]) else min(vals) * 0.5
  load_vals <- pmax(vals[seq_len(K)] - sig2_iso, 1e-4)
  L0 <- vecs[, seq_len(K), drop = FALSE] %*% diag(sqrt(load_vals), K)
  L0 <- make_plt(L0)
  psi_a0 <- pmax(sqrt(pmax(diag(S_A) - rowSums(L0^2), sigma_floor^2)), sigma_floor)

  function(chain_id = 1) {
    Lj <- L0 + matrix(rnorm(Pd * K, 0, 0.02), Pd, K)
    for (q in seq_len(K)) {
      for (p in seq_len(Pd)) if (q > p) Lj[p, q] <- 0
      Lj[q, q] <- abs(Lj[q, q]) + 0.05
    }
    ld <- numeric(K); lf <- numeric(stan_data$N_lambda_a_free); pos <- 1L
    for (q in seq_len(K)) {
      ld[q] <- Lj[q, q]
      for (p in seq(q + 1L, Pd)) { lf[pos] <- Lj[p, q]; pos <- pos + 1L }
    }
    list(
      lambda_a_diag = ld, lambda_a_free = lf, psi_a = psi_a0,
      z_b     = z_b0 + matrix(rnorm(Sv * Pd, 0, 0.05), Sv, Pd),
      sigma_b = sigma_b0, sigma_e = sigma_e0
    )
  }
}

# ── 3. ELPD helpers ───────────────────────────────────────────────────────────

# Phase 1 (base): uniform residual variance per feature, all obs of same player share same d.
# C_p = Sigma_a,pp × J + sigma_e^2 × I  →  uniform-diagonal Woodbury.
# NB: features are treated as independent (uses diagonal Sigma_a,pp, ignoring ΛΛ' cross-feature
# terms). This underestimates joint likelihood but is consistent across K → ranking is preserved.
compute_holdout_ll_base <- function(Y_hold, pi_hold, si_hold, I_hold, P,
                                    Sigma_a_pp, sigma_e2, B_draw) {
  ll <- 0.0
  for (i in seq_len(I_hold)) {
    idx <- which(pi_hold == i)
    T_i <- length(idx)
    if (T_i == 0L) next
    R      <- Y_hold[idx, , drop = FALSE] - B_draw[si_hold[idx], , drop = FALSE]
    sum_r  <- colSums(R); sum_r2 <- colSums(R^2)
    v <- Sigma_a_pp; d <- sigma_e2
    log_det <- (T_i - 1L) * log(d) + log(v * T_i + d)
    quad    <- sum_r2 / d - v / (d * (v * T_i + d)) * sum_r^2
    ll      <- ll + sum(-T_i / 2 * log(2 * pi) - 0.5 * log_det - 0.5 * quad)
  }
  ll
}

# Phase 2 (mv, fixed phi=0.5): non-uniform diagonal — d_{t,p} = sigma_e^2 × s_t^2.
# Rank-1-plus-diagonal Woodbury:
#   |C_p| = prod(d_{t,p}) × (1 + Sigma_a,pp × sum(1/d_{t,p}))
#   r'C^{-1}r = sum_t r_t^2/d_{t,p} − Sigma_a,pp × (sum_t r_t/d_{t,p})^2 / (1 + sum_t Sigma_a,pp/d_{t,p})
# Vectorised over P features using pre-summed quantities per player.
compute_holdout_ll_mv <- function(Y_hold, pi_hold, si_hold, I_hold, P,
                                  Sigma_a_pp, sigma_e2, B_draw, sigma_scale_hold) {
  ll <- 0.0
  for (i in seq_len(I_hold)) {
    idx <- which(pi_hold == i)
    T_i <- length(idx)
    if (T_i == 0L) next
    R <- Y_hold[idx, , drop = FALSE] - B_draw[si_hold[idx], , drop = FALSE]  # T_i × P

    # Observation weights: 1/s_t^2 (same for all features for a given observation)
    inv_s2 <- 1 / sigma_scale_hold[idx]^2    # T_i-vector

    # Per-feature weighted sums (vectorised over P via matrix-vector products)
    # w_t,p = inv_s2[t] / sigma_e2[p]
    # sum_t w_t,p  = sum(inv_s2) / sigma_e2[p]  (scalar × P-vector)
    # sum_t w_t,p × r_t,p = (R' inv_s2)[p] / sigma_e2[p]  (P-vector)
    # sum_t w_t,p × r_t,p^2 = (R^2' inv_s2)[p] / sigma_e2[p]  (P-vector)

    wsum_scalar <- sum(inv_s2)               # scalar
    w_dot_r     <- as.vector(crossprod(R, inv_s2))   / sigma_e2   # P-vector
    w_dot_r2    <- as.vector(crossprod(R^2, inv_s2)) / sigma_e2   # P-vector
    wsum_p      <- wsum_scalar / sigma_e2    # P-vector: sum_t(1/d_{t,p}) = sum(inv_s2)/sigma_e2

    # log|C_p| = T*log(sigma_e2) + 2*sum(log(s_t)) + log(1 + Sigma_a,pp × wsum_p)
    log_s_sum <- sum(log(sigma_scale_hold[idx]))   # scalar = 0.5 sum log(m_ref/min_t)
    log_det   <- T_i * log(sigma_e2) + 2 * log_s_sum + log(1 + Sigma_a_pp * wsum_p)

    # r'C^{-1}r = w_dot_r2 − Sigma_a,pp × w_dot_r^2 / (1 + Sigma_a,pp × wsum_p)
    denom <- 1 + Sigma_a_pp * wsum_p
    quad  <- w_dot_r2 - Sigma_a_pp * w_dot_r^2 / denom

    ll <- ll + sum(-T_i / 2 * log(2 * pi) - 0.5 * log_det - 0.5 * quad)
  }
  ll
}

lse_avg <- function(ll_vec) {
  lse <- max(ll_vec)
  lse + log(mean(exp(ll_vec - lse)))
}

# ── 4. Main K loop ────────────────────────────────────────────────────────────

stan_file_marg <- file.path(paths$stan, "player_season_lowrank_marginalized.stan")
if (!file.exists(stan_file_marg))
  stop("Marginalized Stan model not found: ", stan_file_marg)

model_compiled <- cmdstanr::cmdstan_model(stan_file_marg)
z_b_vars <- as.vector(outer(seq_len(S), seq_len(P),
                            function(s, p) sprintf("z_b[%d,%d]", s, p)))

elpd_results <- vector("list", K_max)

for (K in seq(1L, K_max)) {
  message("\n=== K = ", K, " ===")
  N_free <- lower_tri_free_count(P, K)

  stan_data <- list(
    N = N_train, P = P, I = I_train, S = S,
    Q_a = K, N_lambda_a_free = N_free,
    Y = unname(Y_train),
    player_index = as.integer(pi_train),
    season_index = as.integer(si_train),
    S_count = S_count_train,
    sigma_floor = pipeline$sigma_floor
  )

  init_fn <- make_marginalized_init(stan_data)
  fit_id    <- sprintf("32_marg_k%d", K)
  draws_rds <- file.path(paths$fits, paste0(fit_id, "_draws.rds"))

  draws_cache <- NULL
  if (sampler$reuse_fit && file.exists(draws_rds)) {
    draws_cache <- readRDS(draws_rds)
    message("Reusing Pathfinder draws: ", fit_id)
  }

  if (is.null(draws_cache)) {
    if (!sampler$run_stan) {
      message("BFA_RUN_STAN=false — skipping K=", K)
      elpd_results[[K]] <- tibble(
        K = K, elpd_base = NA_real_, elpd_base_se = NA_real_,
        elpd_mv = NA_real_, elpd_mv_se = NA_real_,
        n_holdout_obs = N_hold, n_holdout_players = I_hold
      )
      next
    }

    message("Running Pathfinder (num_paths=", pf_paths, ", draws=", pf_draws, ") ...")
    fit_pf <- model_compiled$pathfinder(
      data          = stan_data,
      seed          = seed_val + K,
      num_paths     = pf_paths,
      draws         = pf_draws,
      max_lbfgs_iters = pf_iters,
      refresh       = 50,
      init          = lapply(seq_len(pf_paths), init_fn)
    )

    draws_cache <- list(
      sa = posterior::as_draws_matrix(fit_pf$draws("Sigma_a_diag")),
      se = posterior::as_draws_matrix(fit_pf$draws("sigma_e")),
      sb = posterior::as_draws_matrix(fit_pf$draws("sigma_b")),
      zb = posterior::as_draws_matrix(fit_pf$draws("z_b"))
    )
    saveRDS(draws_cache, draws_rds)
    message("Pathfinder complete — draws cached to ", basename(draws_rds))
  }

  # ── Extract draws ─────────────────────────────────────────────────────────
  message("Extracting draws ...")
  draws_sa <- draws_cache$sa   # D × P
  draws_se <- draws_cache$se   # D × P
  draws_sb <- draws_cache$sb   # D × P
  draws_zb <- draws_cache$zb   # D × (S×P)
  D <- nrow(draws_sa)
  message("D = ", D, " draws")

  # ── ELPD: both base and mv per draw ───────────────────────────────────────
  ll_base_vec <- numeric(D)
  ll_mv_vec   <- numeric(D)

  for (d in seq_len(D)) {
    Sigma_a_pp <- as.numeric(draws_sa[d, ])
    sigma_e2   <- as.numeric(draws_se[d, ])^2
    sigma_b_d  <- as.numeric(draws_sb[d, ])
    z_b_mat    <- matrix(as.numeric(draws_zb[d, z_b_vars]), S, P)
    z_b_c      <- sweep(z_b_mat, 2, colMeans(z_b_mat), "-")
    B_draw     <- sweep(z_b_c, 2, sigma_b_d, "*")  # S × P

    ll_base_vec[d] <- compute_holdout_ll_base(
      Y_hold, pi_hold, si_hold, I_hold, P, Sigma_a_pp, sigma_e2, B_draw)
    ll_mv_vec[d] <- compute_holdout_ll_mv(
      Y_hold, pi_hold, si_hold, I_hold, P, Sigma_a_pp, sigma_e2, B_draw,
      sigma_scale_hold)
  }

  elpd_base <- lse_avg(ll_base_vec);  elpd_base_se <- sd(ll_base_vec) / sqrt(D)
  elpd_mv   <- lse_avg(ll_mv_vec);    elpd_mv_se   <- sd(ll_mv_vec)   / sqrt(D)

  message(sprintf("K=%d  ELPD_base=%.1f (SE=%.1f)  ELPD_mv=%.1f (SE=%.1f)",
                  K, elpd_base, elpd_base_se, elpd_mv, elpd_mv_se))

  elpd_results[[K]] <- tibble(
    K = K,
    elpd_base = elpd_base, elpd_base_se = elpd_base_se,
    elpd_mv   = elpd_mv,   elpd_mv_se   = elpd_mv_se,
    n_holdout_obs = N_hold, n_holdout_players = I_hold
  )
}

# ── 5. Output ─────────────────────────────────────────────────────────────────

elpd_tbl <- bind_rows(elpd_results)
write_csv(elpd_tbl, file.path(paths$tables, "32_player_holdout_elpd.csv"))
message("\nPlayer-holdout ELPD table:\n"); print(elpd_tbl)

complete <- elpd_tbl |> filter(!is.na(elpd_mv))

if (nrow(complete) > 0L) {
  K_base_star <- complete$K[which.max(complete$elpd_base)]
  K_mv_star   <- complete$K[which.max(complete$elpd_mv)]
  message(sprintf("\nK* (base model):           K = %d", K_base_star))
  message(sprintf("K* (minutes-scaled phi=0.5): K = %d", K_mv_star))
  if (K_base_star != K_mv_star)
    message("NOTE: K* differs between base and mv — use K_mv* for stage 34.")
  else
    message("K* is consistent across both ELPD metrics — robust selection.")

  writeLines(as.character(K_mv_star),
             file.path(paths$notes, "32_recommended_k.txt"))
  message("Recommended K written to outputs/notes/32_recommended_k.txt")

  # Combined plot
  plot_df <- bind_rows(
    complete |> select(K, elpd = elpd_base, elpd_se = elpd_base_se) |>
      mutate(metric = "Base (constant σ_e)"),
    complete |> select(K, elpd = elpd_mv,   elpd_se = elpd_mv_se)   |>
      mutate(metric = "Minutes-scaled (φ=0.5)")
  )

  p_elpd <- ggplot(plot_df, aes(x = K, y = elpd, colour = metric, group = metric)) +
    geom_errorbar(aes(ymin = elpd - 2 * elpd_se, ymax = elpd + 2 * elpd_se),
                  width = 0.2, position = position_dodge(0.3)) +
    geom_line(linewidth = 0.9, position = position_dodge(0.3)) +
    geom_point(size = 3, position = position_dodge(0.3)) +
    scale_x_continuous(breaks = 0:K_max) +
    scale_colour_manual(values = c("Base (constant σ_e)" = "#2166ac",
                                   "Minutes-scaled (φ=0.5)" = "#d73027"),
                        name = "ELPD metric") +
    geom_vline(xintercept = K_mv_star, linetype = "dashed",
               colour = "#d73027", linewidth = 0.5) +
    labs(
      title    = "Player-holdout ELPD vs number of factors K",
      subtitle = sprintf(
        "%d held-out players (%d obs). Dashed line = recommended K*=%d (minutes-scaled).",
        I_hold, N_hold, K_mv_star),
      x = "K (number of latent factors)",
      y = "Holdout ELPD (log marginal likelihood)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  ggsave(file.path(paths$figures, "32_player_holdout_elpd.png"),
         p_elpd, width = 9, height = 5.5, dpi = 150)
  message("Saved: 32_player_holdout_elpd.png")
}

message("\nStage 32 complete.")
