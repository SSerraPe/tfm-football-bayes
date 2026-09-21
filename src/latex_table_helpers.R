# Shared helper for emitting native LaTeX table fragments (Chapter 5 / results.tex).
#
# Produces just the `tabular` body (via knitr::kable, booktabs style) -- NOT a full
# `table` float with caption/label. results.tex wraps each fragment in its own
# \begin{table}[h] ... \input{tables/<name>.tex} ... \caption{...} \label{...}
# \end{table}, matching the hand-written table convention already used in
# methodology.tex (tab:notation_model), rather than letting kable generate its own
# float (keeps caption wording/positioning under prose control).
#
# Fragments are written to tfm_mesio_latex_2026/tables/<name>.tex.
# (Fixed 2026-09-20, revision pass 3: this pointed at the old, superseded tfm_latex/
# folder -- removed from git and from GitHub, but still present on local disk -- so every
# write_latex_table() call was silently writing to a defunct location while results.tex's
# \input calls kept reading whatever stale content happened to already be sitting in
# tfm_mesio_latex_2026/tables/. Any table fragment not otherwise re-verified in this pass
# should be treated as unconfirmed until regenerated post-fix.)

write_latex_table <- function(df, name, col_names = NULL, digits = 3, align = NULL,
                               longtable = FALSE, escape = TRUE) {
  check_packages(c("knitr"))
  if (!is.null(col_names)) colnames(df) <- col_names
  out_dir <- file.path(paths$scripts, "..", "tfm_mesio_latex_2026", "tables")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, paste0(name, ".tex"))
  tex <- knitr::kable(
    df, format = "latex", booktabs = TRUE, digits = digits,
    align = align, row.names = FALSE, linesep = "", longtable = longtable, escape = escape
  )
  writeLines(as.character(tex), out_path)
  message("Wrote LaTeX table fragment: ", out_path)
  invisible(out_path)
}
