# Shared feature engineering helpers.

safe_divide <- function(num, den, multiplier = 1) {
  ifelse(is.na(den) | den <= 0, NA_real_, multiplier * num / den)
}

weighted_mean_minutes <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  if (!any(keep)) return(NA_real_)
  stats::weighted.mean(x[keep], w[keep])
}

first_non_missing_character <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

median_impute <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  if (is.na(med)) med <- 0
  ifelse(is.na(x), med, x)
}

anchor_free_count <- function(P, Q) {
  if (Q < 1L) return(0L)
  if (Q > P) {
    stop("Cannot identify rank ", Q, " with only ", P, " observed variables.",
         call. = FALSE)
  }
  as.integer((P - Q) * Q)
}

default_anchor_index <- function(P, Q) {
  if (Q > P) {
    stop("Rank cannot exceed the number of observed variables.", call. = FALSE)
  }
  seq_len(Q)
}

audit_scaled_feature_matrix <- function(Y, variable_names,
                                        exact_correlation_threshold = 0.999,
                                        warning_correlation_threshold = 0.95) {
  P <- ncol(Y)
  if (P != length(variable_names)) {
    stop("Feature matrix columns do not match variable_names.", call. = FALSE)
  }

  empirical_rank <- qr(Y, tol = 1e-8)$rank
  if (empirical_rank < P) {
    stop("Scaled feature matrix is rank deficient: P = ", P,
         ", empirical rank = ", empirical_rank, ".", call. = FALSE)
  }

  correlation_matrix <- stats::cor(Y)
  pair_index <- which(upper.tri(correlation_matrix), arr.ind = TRUE)
  pair_table <- tibble::tibble(
    variable_a = variable_names[pair_index[, "row"]],
    variable_b = variable_names[pair_index[, "col"]],
    correlation = correlation_matrix[pair_index],
    abs_correlation = abs(correlation)
  ) |>
    dplyr::arrange(dplyr::desc(abs_correlation), variable_a, variable_b)

  exact_pairs <- pair_table |>
    dplyr::filter(abs_correlation > exact_correlation_threshold)
  if (nrow(exact_pairs) > 0L) {
    stop("Scaled feature matrix has near-exact correlations above ",
         exact_correlation_threshold, ": ",
         paste(paste0(exact_pairs$variable_a, " vs ", exact_pairs$variable_b,
                      " (r = ", signif(exact_pairs$correlation, 6), ")"),
               collapse = "; "),
         call. = FALSE)
  }

  high_pairs <- pair_table |>
    dplyr::filter(abs_correlation >= warning_correlation_threshold)
  if (nrow(high_pairs) > 0L) {
    warning("High feature correlations retained for audit: ",
            paste(paste0(high_pairs$variable_a, " vs ", high_pairs$variable_b,
                         " (r = ", signif(high_pairs$correlation, 4), ")"),
                  collapse = "; "),
            call. = FALSE)
  }

  list(
    summary = tibble::tibble(
      n_observations = nrow(Y),
      n_variables = P,
      empirical_rank = empirical_rank,
      max_abs_correlation = max(pair_table$abs_correlation, na.rm = TRUE),
      exact_correlation_threshold = exact_correlation_threshold,
      warning_correlation_threshold = warning_correlation_threshold,
      high_correlation_pairs = nrow(high_pairs)
    ),
    high_correlation_pairs = high_pairs
  )
}
