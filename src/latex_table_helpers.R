# Shared helper for emitting native LaTeX table fragments (Chapter 5 / results.tex).
#
# Produces just the `tabular` body (via knitr::kable, booktabs style) -- NOT a full
# `table` float with caption/label. results.tex wraps each fragment in its own
# \begin{table}[h] ... \input{tables/<name>.tex} ... \caption{...} \label{...}
# \end{table}, matching the hand-written table convention already used in
# methodology.tex (tab:notation_model), rather than letting kable generate its own
# float (keeps caption wording/positioning under prose control).
#
# Fragments are written to tfm_latex/tables/<name>.tex.

write_latex_table <- function(df, name, col_names = NULL, digits = 3, align = NULL, longtable = FALSE) {
  check_packages(c("knitr"))
  if (!is.null(col_names)) colnames(df) <- col_names
  out_dir <- file.path(paths$scripts, "..", "tfm_latex", "tables")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, paste0(name, ".tex"))
  tex <- knitr::kable(
    df, format = "latex", booktabs = TRUE, digits = digits,
    align = align, row.names = FALSE, linesep = "", longtable = longtable
  )
  writeLines(as.character(tex), out_path)
  message("Wrote LaTeX table fragment: ", out_path)
  invisible(out_path)
}
