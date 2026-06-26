# -----------------------------
# PCA-like postprocessing of posterior draws
# -----------------------------

lambda_arr = fit$draws(variables = "Lambda", format = "draws_array")
F_arr = fit$draws(variables = "F", format = "draws_array")

n_iter = dim(lambda_arr)[1]
n_chains = dim(lambda_arr)[2]

all_draws = expand.grid(iter = 1:n_iter, chain = 1:n_chains)

set.seed(123)
n_use = min(1000, nrow(all_draws))
use_draws = all_draws[sample(seq_len(nrow(all_draws)), n_use), ]

Lambda_vars = as.vector(
  outer(
    1:P,
    1:K_fit,
    function(j, h) paste0("Lambda[", j, ",", h, "]")
  )
)

F_vars = as.vector(
  outer(
    1:I,
    1:K_fit,
    function(i, h) paste0("F[", i, ",", h, "]")
  )
)

# -----------------------------
# Helper functions
# -----------------------------

fix_signs = function(L) {
  for (h in 1:ncol(L)) {
    jmax = which.max(abs(L[, h]))
    if (L[jmax, h]<0) {
      L[, h] = -L[, h]
    }
  }
  L
}

pca_from_B = function(B, K) {
  eg = eigen(B, symmetric = TRUE)
  vals = pmax(eg$values[1:K], 0)
  U = eg$vectors[, 1:K, drop = FALSE]
  L = U%*%diag(sqrt(vals), K)
  list(U = U, values = vals, loadings = L)
}

align_to_reference = function(L, G, L_ref) {
  for (h in 1:ncol(L)) {
    if (sum(L[, h]*L_ref[, h])<0) {
      L[, h] = -L[, h]
      G[, h] = -G[, h]
    }
  }
  list(loadings = L, scores = G)
}

# -----------------------------
# Reference orientation from posterior mean B
# -----------------------------

B_ref = matrix(0, P, P)

for (s in 1:n_use) {
  it = use_draws$iter[s]
  ch = use_draws$chain[s]
  
  L_s = matrix(as.numeric(lambda_arr[it, ch, Lambda_vars]), P, K_fit)
  B_ref = B_ref+(L_s%*%t(L_s))/n_use
}

pca_ref = pca_from_B(B_ref, K_fit)
L_ref = fix_signs(pca_ref$loadings)

# -----------------------------
# Draw-wise PCA loadings and scores
# -----------------------------

L_pca_draws = array(NA, dim = c(n_use, P, K_fit))
G_pca_draws = array(NA, dim = c(n_use, I, K_fit))
eig_draws = matrix(NA, n_use, K_fit)

for (s in 1:n_use) {
  it = use_draws$iter[s]
  ch = use_draws$chain[s]
  
  L_tri = matrix(as.numeric(lambda_arr[it, ch, Lambda_vars]), P, K_fit)
  F_s = matrix(as.numeric(F_arr[it, ch, F_vars]), I, K_fit)
  
  B_s = L_tri%*%t(L_tri)
  pca_s = pca_from_B(B_s, K_fit)
  
  U_s = pca_s$U
  vals_s = pca_s$values
  L_pca_s = pca_s$loadings
  
  C_s = F_s%*%t(L_tri)
  G_s = C_s%*%U_s%*%diag(1/sqrt(pmax(vals_s, 1e-10)), K_fit)
  
  aligned = align_to_reference(L_pca_s, G_s, L_ref)
  
  L_pca_draws[s,,] = aligned$loadings
  G_pca_draws[s,,] = aligned$scores
  eig_draws[s, ] = vals_s
}

L_pca_mean = apply(L_pca_draws, c(2, 3), mean)
L_pca_q05 = apply(L_pca_draws, c(2, 3), quantile, probs = 0.05)
L_pca_q95 = apply(L_pca_draws, c(2, 3), quantile, probs = 0.95)

G_pca_mean = apply(G_pca_draws, c(2, 3), mean)
G_pca_q05 = apply(G_pca_draws, c(2, 3), quantile, probs = 0.05)
G_pca_q95 = apply(G_pca_draws, c(2, 3), quantile, probs = 0.95)

prop_common = eig_draws/rowSums(eig_draws)

# -----------------------------
# Truth in PCA orientation, for simulation checks
# -----------------------------

pca_true = pca_from_B(B_a_true, K_fit)
L_pca_true = pca_true$loadings

# Align truth to posterior reference.
dummy_scores = matrix(0, I, K_fit)
aligned_true = align_to_reference(L_pca_true, dummy_scores, L_ref)
L_pca_true = aligned_true$loadings

C_true = F_true%*%t(Lambda_true_std)
G_pca_true = C_true%*%pca_true$U%*%diag(1/sqrt(pmax(pca_true$values, 1e-10)), K_fit)

for (h in 1:K_fit) {
  if (sum(L_pca_true[, h]*L_ref[, h])<0) {
    G_pca_true[, h] = -G_pca_true[, h]
  }
}

prop_common_true = pca_true$values/sum(pca_true$values)

# -----------------------------
# Numerical summaries
# -----------------------------

cat("\nPCA-oriented loading recovery:\n")
print(c(
  cor_pca_loadings = cor(as.vector(L_pca_true), as.vector(L_pca_mean)),
  rmse_pca_loadings = sqrt(mean((L_pca_mean-L_pca_true)^2))
))

cat("\nPCA-oriented score recovery:\n")
for (h in 1:K_fit) {
  cat(
    "factor", h,
    ": cor =",
    round(cor(G_pca_true[, h], G_pca_mean[, h]), 3),
    ", RMSE =",
    round(sqrt(mean((G_pca_mean[, h]-G_pca_true[, h])^2)), 3),
    "\n"
  )
}

cat("\nProportion of common player variance explained:\n")
print(apply(prop_common, 2, quantile, probs = c(0.05, 0.5, 0.95)))

# -----------------------------
# Plots: PCA-oriented loadings
# -----------------------------

par(mfrow = c(1, K_fit), mar = c(5, 5, 3, 1))

for (h in 1:K_fit) {
  ylim_h = range(L_pca_true[, h], L_pca_q05[, h], L_pca_q95[, h])
  
  plot(
    1:P,
    L_pca_true[, h],
    type = "b",
    pch = 19,
    ylim = ylim_h,
    xlab = "Variable",
    ylab = "Loading",
    main = paste("PCA factor", h)
  )
  
  lines(1:P, L_pca_mean[, h], type = "b", pch = 19, lty = 2)
  segments(1:P, L_pca_q05[, h], 1:P, L_pca_q95[, h])
  abline(h = 0, lty = 3)
  
  legend(
    "topright",
    legend = c("Truth", "Posterior mean"),
    lty = c(1, 2),
    pch = 19,
    bty = "n"
  )
}

par(mfrow = c(1, 1))

plot(
  as.vector(L_pca_true),
  as.vector(L_pca_mean),
  pch = 19,
  xlab = "True PCA-oriented loading",
  ylab = "Posterior mean PCA-oriented loading",
  main = "PCA-oriented loading recovery"
)
abline(0, 1, lty = 2)

# -----------------------------
# Plots: common variance explained
# -----------------------------

boxplot(
  prop_common,
  xlab = "PCA factor",
  ylab = "Proportion",
  main = "Common player variance explained"
)
points(1:K_fit, prop_common_true, pch = 19)

# -----------------------------
# Plots: PCA-oriented player scores
# -----------------------------

par(mfrow = c(1, K_fit), mar = c(5, 5, 3, 1))

for (h in 1:K_fit) {
  plot(
    G_pca_true[, h],
    G_pca_mean[, h],
    pch = 19,
    xlab = paste("True PCA score", h),
    ylab = paste("Posterior mean PCA score", h),
    main = paste("Score recovery, factor", h)
  )
  abline(0, 1, lty = 2)
}

par(mfrow = c(1, 1))

if (K_fit>=2) {
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
  
  plot(
    G_pca_true[, 1],
    G_pca_true[, 2],
    pch = 19,
    xlab = "True PCA score 1",
    ylab = "True PCA score 2",
    main = "True player scores"
  )
  
  plot(
    G_pca_mean[, 1],
    G_pca_mean[, 2],
    pch = 19,
    xlab = "Posterior mean PCA score 1",
    ylab = "Posterior mean PCA score 2",
    main = "Estimated player scores"
  )
  
  par(mfrow = c(1, 1))
}