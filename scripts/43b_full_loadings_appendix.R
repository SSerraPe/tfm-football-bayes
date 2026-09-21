# Stage 43b -- Full canonical-loading appendix table (Results-chapter §6 revision).
#
# Table 5 in results.tex (tab:pc-comparison) shows only the four largest reliable loadings
# per pole per PC, for readability in the main text. This script builds the exhaustive
# companion: every reliable (90% CI excludes zero) loading on every PC, split by sign and
# grouped by feature category, as a bulleted longtable for the appendix.
#
# Reads:  outputs/tables/31_pca_loading_ci.csv (pc, feature, group, loading, ci_lo, ci_hi, reliable)
# Writes: tfm_mesio_latex_2026/tables/43b_full_loadings_appendix.tex (a \begin{longtable}...
#         \end{longtable} fragment, \input from the appendix in main.tex)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/43b_full_loadings_appendix.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
source(file.path(model_root, "src", "loading_visualization.R"))

suppressPackageStartupMessages({ library(readr); library(dplyr) })

d <- read_csv(file.path(paths$tables, "31_pca_loading_ci.csv"), show_col_types = FALSE)

clean_feature <- function(x) gsub("_", " ", gsub("^per90_|^rate_", "", x))
esc <- function(x) gsub("_", "\\\\_", x)  # only underscores appear in this content

GROUP_ORDER <- c("Discipline", "Defending", "Attacking output", "Chance creation",
                  "Progressive actions", "Passing volume", "Passing quality", "Dueling")

PC_NAMES <- c(PC1 = "Defensive engagement", PC2 = "Technical orientation", PC3 = "Possession retention")
PC_POLES <- list(
  PC1 = c(neg = "attacking",  pos = "defensive"),
  PC2 = c(neg = "aerial",     pos = "technical"),
  PC3 = c(neg = "crossing",   pos = "retention")
)

d <- d |>
  mutate(group = factor(group, levels = GROUP_ORDER), feature_clean = clean_feature(feature))

# One row PER GROUP (not one row per PC/sign with every group crammed into a single cell):
# longtable can only break between rows, never inside one, and a sign's full group list
# (up to 8 groups, 33 features for PC2's positive pole) would risk overflowing a page if it
# were a single cell/row.
group_rows <- function(sub_pc, sign) {
  rows <- if (sign == "pos") sub_pc |> filter(loading > 0) else sub_pc |> filter(loading < 0)
  if (nrow(rows) == 0) return(character(0))
  out <- character(0)
  for (g in GROUP_ORDER) {
    gp <- rows |> filter(group == g) |> arrange(desc(abs(loading)))
    if (nrow(gp) == 0) next
    items <- paste0("\\item ", esc(gp$feature_clean), " (", sprintf("%.2f", gp$loading), ")", collapse = "\n    ")
    out <- c(out, sprintf("%s & \\begin{itemize}[nosep,leftmargin=1.1em,topsep=1pt]\n    %s\n  \\end{itemize}\\\\",
                           g, items))
  }
  out
}

out_lines <- c(
  "\\begin{longtable}{p{2.6cm}p{13cm}}",
  "\\caption{Every reliable (90\\% credible interval excluding zero) posterior-mean loading on each canonical",
  "factor, by sign and feature group. Table~\\ref{tab:pc-comparison} in the main text shows only the four",
  "largest per pole; this is the full list underlying Figure~\\ref{fig:loading-heatmap}.}",
  "\\label{tab:full-loadings}\\\\",
  "\\toprule",
  "Group & Loadings \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{2}{l}{\\small\\itshape (Table~\\ref{tab:full-loadings} continued from previous page)}\\\\",
  "\\toprule",
  "Group & Loadings \\\\",
  "\\midrule",
  "\\endhead",
  "\\midrule",
  "\\multicolumn{2}{r}{\\small\\itshape continued on next page}\\\\",
  "\\endfoot",
  "\\bottomrule",
  "\\endlastfoot"
)

for (pcx in c("PC1", "PC2", "PC3")) {
  sub_pc <- d |> filter(pc == pcx, reliable)
  poles <- PC_POLES[[pcx]]
  pos_rows <- group_rows(sub_pc, "pos")
  neg_rows <- group_rows(sub_pc, "neg")
  # add extra trailing space to the last row of each block via \\[4pt] instead of a plain \\
  if (length(pos_rows) > 0) pos_rows[length(pos_rows)] <- sub("\\\\\\\\$", "\\\\\\\\[4pt]", pos_rows[length(pos_rows)])
  if (length(neg_rows) > 0) neg_rows[length(neg_rows)] <- sub("\\\\\\\\$", "\\\\\\\\[8pt]", neg_rows[length(neg_rows)])
  # Section headers as plain 2-column rows (left cell blank) rather than \multicolumn --
  # \multicolumn inside a table body has been observed elsewhere in this document to corrupt
  # main.aux (a hyperref/aux-writing interaction, not specific to this content), so it is
  # avoided here even though it would give a cleaner full-width visual span.
  out_lines <- c(out_lines,
    sprintf(" & \\rule{0pt}{2.2ex}\\textbf{%s: %s}\\\\", pcx, PC_NAMES[[pcx]]),
    sprintf(" & \\textit{Positive loadings (%s pole)}\\\\", poles[["pos"]]),
    pos_rows,
    sprintf(" & \\textit{Negative loadings (%s pole)}\\\\", poles[["neg"]]),
    neg_rows
  )
}

out_lines <- c(out_lines, "\\end{longtable}")

out_dir <- file.path(paths$scripts, "..", "tfm_mesio_latex_2026", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "43b_full_loadings_appendix.tex")
writeLines(out_lines, out_path)
message("Wrote: ", out_path)
