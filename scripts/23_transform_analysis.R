# Stage 23: Variable transform analysis (Issue 8, transform half).
#
# Goal: identify which per90_* and rate_* features have heavy right skew in the
# pre-scaled data, and assess which variance-stabilising transform (log1p / sqrt /
# asinh) would best compress the tail before Z-scaling.
#
# Mechanism (why transforms help):
#   Many per90_* features are per-90-minute counts of rare events (goals, yellow
#   cards, tackles). Their raw distributions are right-skewed (Poisson-like or
#   zero-inflated). When Z-scaled, the long right tail becomes extreme positive
#   values in Y that the model must absorb with heavy-tailed residuals or large
#   factor weights. A variance-stabilising transform applied BEFORE Z-scaling
#   compresses the tail, reducing the magnitude of these extreme Z-scores.
#   The expected consequence (if the transform is appropriate) is that the
#   standardised residuals become less extreme, and the model's degrees-of-
#   freedom parameter ν may shift upward — i.e., the t-distribution's tail
#   was partly compensating for distributional misspecification, not pure
#   observation noise. Whether ν actually rises requires a refit to confirm.
#
# REFIT GATE: applying transforms requires re-running stages 01+ (feature
# engineering with transforms applied before Z-scaling). This script is ANALYSIS
# ONLY. Report findings to professor; get a go-ahead before rebuilding.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/23_transform_analysis.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Base R skewness (type 1 moment estimator) and excess kurtosis
skewness_r <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  m3 <- sum((x - m)^3) / n
  m2 <- sum((x - m)^2) / n
  m3 / m2^(3/2)
}

kurtosis_r <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  m4 <- mean((x - m)^4)
  m2 <- mean((x - m)^2)
  m4 / m2^2 - 3   # excess kurtosis
}

message("=== Stage 23: Variable transform analysis ===")

# ---- Load raw (pre-scaled) data -------------------------------------------------

model_objects <- readRDS(file.path(paths$processed, "model_objects.rds"))

# model_objects$Y_raw is the unscaled feature matrix (N × P)
Y_raw <- model_objects$Y_raw
feat_names <- model_objects$variable_names
colnames(Y_raw) <- feat_names

message("Loaded Y_raw: ", nrow(Y_raw), " rows × ", ncol(Y_raw), " features")

# Fallback: if Y_raw is absent, load from model_feature_matrix_unscaled.csv
if (is.null(Y_raw) || nrow(Y_raw) == 0) {
  message("Y_raw not in model_objects; loading from model_feature_matrix_unscaled.csv")
  Y_raw <- as.matrix(
    read_csv(file.path(paths$processed, "model_feature_matrix_unscaled.csv"),
             show_col_types = FALSE)
  )
  colnames(Y_raw) <- feat_names
}

# ---- Compute skewness, kurtosis, and minimum value per feature ------------------

compute_stats <- function(x) {
  x <- x[is.finite(x)]
  data.frame(
    n          = length(x),
    mean       = mean(x),
    sd         = sd(x),
    min_val    = min(x),
    p01        = quantile(x, 0.01),
    p25        = quantile(x, 0.25),
    median_val = median(x),
    p75        = quantile(x, 0.75),
    p99        = quantile(x, 0.99),
    max_val    = max(x),
    skewness   = skewness_r(x),
    kurtosis   = kurtosis_r(x)
  )
}

stats_list <- lapply(feat_names, function(f) {
  df <- compute_stats(Y_raw[, f])
  df$feature <- f
  df
})
feat_stats <- bind_rows(stats_list) |>
  mutate(
    feature_type = ifelse(grepl("^rate_", feature), "rate", "per90"),
    # Recommend transform based on skewness and data range:
    #   log1p  — right-skewed count features with min ≥ 0
    #   asinh  — right-skewed but can include 0 or small negatives after rounding
    #   sqrt   — moderate skew (1–2), min ≥ 0
    #   none   — already near-symmetric (|skewness| < 1)
    recommended_transform = dplyr::case_when(
      abs(skewness) < 1                          ~ "none",
      skewness >= 1 & skewness < 2 & min_val >= 0 ~ "sqrt",
      skewness >= 2 & min_val >= 0               ~ "log1p",
      skewness >= 1 & min_val < 0               ~ "asinh",
      TRUE                                       ~ "none"
    ),
    mechanism = dplyr::case_when(
      recommended_transform == "log1p"  ~
        "Sparse count; log1p compresses right tail, equidistances small counts",
      recommended_transform == "sqrt"   ~
        "Moderate count skew; sqrt is variance-stabilising for Poisson-like distributions",
      recommended_transform == "asinh"  ~
        "Right-skewed but can include near-zero or small negatives; asinh handles both",
      TRUE                              ~ "Distribution already near-symmetric; no transform needed"
    )
  ) |>
  select(feature, feature_type, skewness, kurtosis, min_val, max_val,
         recommended_transform, mechanism) |>
  arrange(desc(skewness))

write_csv(feat_stats, file.path(paths$tables, "23_feature_skewness.csv"))
message("Saved 23_feature_skewness.csv")

# ---- Print summary --------------------------------------------------------------

message("\n--- Features with skewness > 2 (strong right skew) ---")
print(as.data.frame(feat_stats |> filter(skewness > 2) |>
                      select(feature, skewness, min_val, recommended_transform, mechanism)))

message("\n--- Features with skewness 1–2 (moderate skew) ---")
print(as.data.frame(feat_stats |> filter(skewness >= 1, skewness <= 2) |>
                      select(feature, skewness, min_val, recommended_transform)))

message("\n--- Transform recommendation counts ---")
print(table(feat_stats$recommended_transform))

# ---- Visual comparison: pre- and post-transform distributions -------------------
# Only for the 6 highest-skew features (to keep the output readable).

high_skew <- feat_stats |> filter(skewness > 1) |> head(6) |> pull(feature)

if (length(high_skew) > 0) {
  apply_transform <- function(x, transform) {
    switch(transform,
      "log1p" = log1p(x),
      "sqrt"  = sqrt(pmax(x, 0)),
      "asinh" = asinh(x),
      x
    )
  }

  plot_data <- lapply(high_skew, function(f) {
    t <- feat_stats$recommended_transform[feat_stats$feature == f]
    x_raw <- Y_raw[, f]
    x_trans <- apply_transform(x_raw, t)
    data.frame(
      feature   = f,
      transform = t,
      raw       = x_raw,
      transformed = x_trans
    )
  }) |> bind_rows()

  # Long format for faceting
  plot_long <- plot_data |>
    pivot_longer(c(raw, transformed), names_to = "version", values_to = "value") |>
    mutate(
      version = factor(version, levels = c("raw", "transformed"),
                       labels = c("Raw", "After transform"))
    )

  p <- ggplot(plot_long, aes(x = value)) +
    geom_histogram(bins = 40, fill = "steelblue", colour = "white", alpha = 0.8) +
    facet_grid(feature ~ version, scales = "free") +
    labs(
      title    = "Distribution before and after recommended transform",
      subtitle = "Transforms applied before Z-scaling (analysis only — no refit yet)",
      x = "Feature value", y = "Count"
    ) +
    theme_minimal(base_size = 9) +
    theme(strip.text.y = element_text(size = 7, angle = 0))

  fig_path <- file.path(paths$figures, "23_transform_examples.png")
  ggsave(fig_path, p, width = 10, height = max(4, length(high_skew) * 1.5),
         dpi = 150)
  message("Saved 23_transform_examples.png")
} else {
  message("No features with skewness > 1; no transform plot produced.")
}

message("\nREFIT GATE: Applying these transforms requires re-running stage 01+ with")
message("the transforms applied before Z-scaling in 01_read_clean_feature_engineering.R.")
message("This changes Y_scaled.csv and every downstream stage. Confirm before running.")

message("\nStage 23 complete.")
