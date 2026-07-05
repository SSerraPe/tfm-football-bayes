suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

fs  <- read.csv("outputs/tables/29_residual_feature_summary.csv",
                stringsAsFactors = FALSE)
tri <- read.csv("outputs/tables/29_residual_skew_triangulation.csv",
                stringsAsFactors = FALSE)

# Defensive: drop any row where verdict is NA before plotting
tri_plot <- tri |>
  filter(!is.na(verdict)) |>
  mutate(feature = factor(feature, levels = rev(fs$feature[order(-fs$p99_abs)])))

cat("Verdict distribution used for plot:\n")
print(table(tri_plot$verdict))
cat("Any NA in feature factor:", any(is.na(tri_plot$feature)), "\n")

tau_99 <- tri$tau_99_val[1]
nu_hat <- 4.9023   # from 28_diagnostics_cache.csv

verdict_colours <- c(
  "contamination-driven" = "#d73027",
  "season-driven"        = "#fc8d59",
  "skew-residual"        = "#fee090",
  "mixed"                = "#abd9e9",
  "unremarkable_p99"     = "#d0d0d0"
)

p_bar <- ggplot(tri_plot, aes(x = p99_abs, y = feature, fill = verdict)) +
  geom_col() +
  geom_vline(xintercept = tau_99, linetype = "dashed", colour = "black",
             linewidth = 0.8) +
  annotate("text", x = tau_99 + 0.02, y = 1, hjust = 0, size = 3,
           label = sprintf("τₙ₉ = %.3f", tau_99)) +
  scale_fill_manual(values = verdict_colours, name = "Verdict",
                    na.value = NA, drop = TRUE) +
  labs(
    title    = "Stage 28 — Residual tail by feature (P=48, K=3, t-model)",
    subtitle = sprintf(
      "Threshold τₙₙ = %.4f at ν = %.3f; above dashed line = heavier than model predicts",
      tau_99, nu_hat),
    x = "p99 of |zₛₜₑ|", y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))

ggplot2::ggsave("outputs/figures/29_residual_p99_barchart.png", p_bar,
                width = 12, height = 10, dpi = 150)
message("Saved: outputs/figures/29_residual_p99_barchart.png")
