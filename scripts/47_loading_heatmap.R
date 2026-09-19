# Stage 47 -- Grouped canonical-loading heatmap (Results-chapter rewrite).
#
# Replaces the crowded loading-CI dot plot (stage 31's 31_pca_loading_ci_plot.png) and the
# three separate top-loader tables with one figure: one panel per feature group, 3 columns
# (PC1/PC2/PC3) x n-feature rows, square cells, posterior-mean loading on a diverging scale,
# a solid group-colour strip in place of per-row colour, no in-cell numbers, one shared
# colour bar, a legend below.
#
# Group scheme: the REAL 7-group feature_group_lookup() (src/loading_visualization.R),
# restricted to the 48 features actually in the model -- NOT the mockup brief's synthetic
# 8-group scheme ("Shooting rates", the 8th category, has zero members among the 48 real
# features; both its features were dropped in the 53->48 rebuild). Groups ordered so PC1's
# sign runs monotonically through the reading order (down column 1, then column 2), split
# 27/21 rows across the two columns (Defending+Passing quality+Passing volume | Dueling+
# Discipline+Progressive actions+Attacking output) rather than an even 4/3 group-count
# split, which would have put 34 rows in one column and 14 in the other.
#
# Reliability encoding: a thin dark outline on cells whose 90% CI excludes zero (137/144
# cells, 95.1%), no outline on the rest -- outlining reliable cells (the large majority)
# makes the handful of unreliable ones the visible exception, rather than dimming the
# unreliable minority and fighting the diverging colour scale (a faded cell already reads
# as "near zero," conflating magnitude with reliability).
#
# Reads:  outputs/tables/31_pca_loading_ci.csv (pc, feature, group, loading, ci_lo, ci_hi, reliable)
# Writes: outputs/figures/47_loading_heatmap.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/47_loading_heatmap.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
source(file.path(model_root, "src", "loading_visualization.R"))
check_packages(c(required_base_packages, "grid"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(grid) })

d <- read_csv(file.path(paths$tables, "31_pca_loading_ci.csv"), show_col_types = FALSE)

# ---- group order (real 7-group scheme, PC1-monotonic) -----------------------------------
GROUP_ORDER_COL1 <- c("Defending", "Passing quality", "Passing volume")
GROUP_ORDER_COL2 <- c("Dueling", "Discipline", "Progressive actions", "Attacking output")
ALL_GROUPS <- c(GROUP_ORDER_COL1, GROUP_ORDER_COL2)

gcols <- group_colours()
stopifnot(all(ALL_GROUPS %in% names(gcols)))

clean_feature <- function(x) gsub("_", " ", gsub("^per90_|^rate_", "", x))

vmax <- max(abs(d$loading))
message(sprintf("Colour scale pinned to real max |loading| = %.4f (feature: %s)",
                 vmax, d$feature[which.max(abs(d$loading))]))

# ---- build per-group ordered feature lists + matrices ------------------------------------
group_data <- lapply(ALL_GROUPS, function(g) {
  sub <- d |> filter(group == g)
  # order within panel by |PC1 loading| descending (consistent visual anchor: the axis
  # driving the group-ordering decision also orders each panel's rows)
  pc1_ord <- sub |> filter(pc == "PC1") |> arrange(desc(abs(loading))) |> pull(feature)
  feats <- pc1_ord
  M <- matrix(NA_real_, nrow = length(feats), ncol = 3, dimnames = list(feats, c("PC1", "PC2", "PC3")))
  R <- matrix(FALSE, nrow = length(feats), ncol = 3, dimnames = list(feats, c("PC1", "PC2", "PC3")))
  for (pcx in c("PC1", "PC2", "PC3")) {
    row <- sub |> filter(pc == pcx)
    M[row$feature, pcx] <- row$loading
    R[row$feature, pcx] <- row$reliable
  }
  list(group = g, features = feats, labels = clean_feature(feats), M = M, R = R, n = length(feats))
})
names(group_data) <- ALL_GROUPS

n_col1 <- sum(vapply(group_data[GROUP_ORDER_COL1], `[[`, integer(1), "n"))
n_col2 <- sum(vapply(group_data[GROUP_ORDER_COL2], `[[`, integer(1), "n"))
message(sprintf("Column split: column 1 = %d rows (%s), column 2 = %d rows (%s)",
                 n_col1, paste(GROUP_ORDER_COL1, collapse = "+"),
                 n_col2, paste(GROUP_ORDER_COL2, collapse = "+")))

# ---- colour ramp: RdBu_r-equivalent (blue = negative, white = zero, red = positive) ------
rdbu_r_anchors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                     "#F7F7F7",
                     "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
ramp <- colorRampPalette(rdbu_r_anchors, space = "Lab")
loading_to_colour <- function(x, vmax_) {
  t <- pmin(pmax((x / vmax_ + 1) / 2, 0), 1)
  ramp(256)[pmax(1, pmin(256, round(t * 255) + 1))]
}

# ---- geometry (inches) --------------------------------------------------------------------
# Canvas widened from the spec's 6.4in to 7.3in so the colour bar and the two longest
# legend entries ("Progressive actions", "Attacking output") have room -- cell size, panel
# heights, and the square-cell constraint are unaffected; only the right-hand margin grows.
FW <- 7.3; FH <- 8.6
cell <- 0.26            # square cell, in
gap  <- 0.09            # gap between stacked panels within a column, in
top_margin <- 0.40      # space for the title
bottom_margin <- 0.85   # space for the legend
strip_w <- 0.10         # colour-strip width, in
label_gap <- 0.06       # gap between labels and strip

avail_h <- FH - top_margin - bottom_margin
col1_h <- n_col1 * cell + (length(GROUP_ORDER_COL1) - 1) * gap
col2_h <- n_col2 * cell + (length(GROUP_ORDER_COL2) - 1) * gap
stopifnot(col1_h <= avail_h + 1e-6)  # sanity: must fit

# Horizontal layout: label zone | strip | heatmap(3 cells), for each of 2 columns, plus a
# right-hand colour-bar zone.
label_w <- 1.55
hm_w <- 3 * cell
col_block_w <- label_w + label_gap + strip_w + hm_w
col1_x0 <- 0.12
col2_x0 <- col1_x0 + col_block_w + 0.55
cbar_x0 <- col2_x0 + col_block_w + 0.35
stopifnot(cbar_x0 + 0.22 <= FW)

out_path <- file.path(paths$figures, "47_loading_heatmap.png")
png(out_path, width = FW, height = FH, units = "in", res = 220)
grid.newpage()

draw_panel <- function(gname, x0, y_top) {
  gd <- group_data[[gname]]
  n  <- gd$n
  h  <- n * cell
  y0 <- y_top - h  # bottom of this panel, inches from bottom

  # colour strip
  strip_vp <- viewport(x = unit(x0 + label_w + label_gap, "in"), y = unit(y0, "in"),
                        width = unit(strip_w, "in"), height = unit(h, "in"),
                        just = c("left", "bottom"))
  pushViewport(strip_vp)
  grid.rect(gp = gpar(fill = gcols[[gname]], col = NA))
  popViewport()

  # heatmap cells
  hm_x0 <- x0 + label_w + label_gap + strip_w
  for (i in seq_len(n)) {
    row_y0 <- y0 + h - i * cell   # top row (i=1) drawn at the top
    for (j in 1:3) {
      cell_x0 <- hm_x0 + (j - 1) * cell
      col <- loading_to_colour(gd$M[i, j], vmax)
      vp <- viewport(x = unit(cell_x0, "in"), y = unit(row_y0, "in"),
                      width = unit(cell, "in"), height = unit(cell, "in"),
                      just = c("left", "bottom"))
      pushViewport(vp)
      grid.rect(gp = gpar(fill = col, col = "grey92", lwd = 0.3))
      if (isTRUE(gd$R[i, j])) {
        grid.rect(width = unit(1, "npc") - unit(1.1, "pt"), height = unit(1, "npc") - unit(1.1, "pt"),
                   gp = gpar(fill = NA, col = "black", lwd = 1.0))
      }
      popViewport()
    }
    # feature label, right-aligned, ending just left of the strip
    lab_vp <- viewport(x = unit(x0, "in"), y = unit(row_y0, "in"),
                        width = unit(label_w, "in"), height = unit(cell, "in"),
                        just = c("left", "bottom"))
    pushViewport(lab_vp)
    grid.text(gd$labels[i], x = unit(1, "npc"), y = unit(0.5, "npc"),
               just = c("right", "center"), gp = gpar(fontsize = 6.8, fontfamily = "serif"))
    popViewport()
  }

  invisible(y0)
}

draw_column <- function(groups, x0, top_y, is_first_col) {
  y_cursor <- top_y
  for (k in seq_along(groups)) {
    gname <- groups[k]
    y0 <- draw_panel(gname, x0, y_cursor)
    if (k == 1) {
      # PC1/PC2/PC3 column headers on top of the first panel in this column
      hm_x0 <- x0 + label_w + label_gap + strip_w
      for (j in 1:3) {
        vp <- viewport(x = unit(hm_x0 + (j - 1) * cell, "in"), y = unit(y_cursor, "in"),
                        width = unit(cell, "in"), height = unit(0.16, "in"),
                        just = c("left", "bottom"))
        pushViewport(vp)
        grid.text(paste0("PC", j), y = unit(0.15, "npc"), gp = gpar(fontsize = 7, fontfamily = "serif"))
        popViewport()
      }
    }
    y_cursor <- y0 - gap
  }
}

top_y <- FH - top_margin
draw_column(GROUP_ORDER_COL1, col1_x0, top_y, TRUE)
draw_column(GROUP_ORDER_COL2, col2_x0, top_y, FALSE)

# ---- shared colour bar --------------------------------------------------------------------
cbar_h <- min(col1_h, 3.2)
cbar_y0 <- (FH - bottom_margin - top_margin - cbar_h) / 2 + bottom_margin - 0.15
n_grad <- 256
grad_cols <- loading_to_colour(seq(-vmax, vmax, length.out = n_grad), vmax)
cbar_vp <- viewport(x = unit(cbar_x0, "in"), y = unit(cbar_y0, "in"),
                     width = unit(0.16, "in"), height = unit(cbar_h, "in"), just = c("left", "bottom"))
pushViewport(cbar_vp)
grid.raster(matrix(rev(grad_cols), ncol = 1), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
grid.rect(gp = gpar(fill = NA, col = "grey40", lwd = 0.4))
grid.text(sprintf("+%.2f", vmax), x = unit(1, "npc") + unit(0.05, "in"), y = unit(1, "npc"),
           just = c("left", "top"), gp = gpar(fontsize = 6.5, fontfamily = "serif"))
grid.text(sprintf("%.2f", -vmax), x = unit(1, "npc") + unit(0.05, "in"), y = unit(0, "npc"),
           just = c("left", "bottom"), gp = gpar(fontsize = 6.5, fontfamily = "serif"))
grid.text("0", x = unit(1, "npc") + unit(0.05, "in"), y = unit(0.5, "npc"),
           just = c("left", "center"), gp = gpar(fontsize = 6.5, fontfamily = "serif"))
popViewport()
grid.text("posterior-mean\nloading", x = unit(cbar_x0, "in"), y = unit(cbar_y0 + cbar_h + 0.12, "in"),
           just = c("left", "bottom"), gp = gpar(fontsize = 7.2, fontfamily = "serif", lineheight = 0.95))

# ---- legend (below the panels; only place group names appear) -----------------------------
# Four groups on the first row, three on the second -- x-positions computed from measured
# text width (in inches) rather than an even grid, since group names vary a lot in length
# ("Dueling" vs "Progressive actions") and an even grid either overlaps or wastes space.
leg_font <- gpar(fontsize = 7.6, fontfamily = "serif")
swatch_w <- 0.15; swatch_gap <- 0.07; entry_gap <- 0.30
leg_rows <- list(ALL_GROUPS[1:4], ALL_GROUPS[5:7])
row_y <- c(bottom_margin - 0.28, bottom_margin - 0.58)
for (r in 1:2) {
  gs <- leg_rows[[r]]
  widths <- vapply(gs, function(g) convertWidth(grobWidth(textGrob(g, gp = leg_font)), "in", valueOnly = TRUE), numeric(1))
  entry_w <- swatch_w + swatch_gap + widths
  total_w <- sum(entry_w) + entry_gap * (length(gs) - 1)
  x0 <- (FW - total_w) / 2
  xcur <- x0
  for (k in seq_along(gs)) {
    sw_vp <- viewport(x = unit(xcur, "in"), y = unit(row_y[r], "in"),
                        width = unit(swatch_w, "in"), height = unit(swatch_w, "in"), just = c("left", "center"))
    pushViewport(sw_vp)
    grid.rect(gp = gpar(fill = gcols[[gs[k]]], col = NA))
    popViewport()
    grid.text(gs[k], x = unit(xcur + swatch_w + swatch_gap, "in"), y = unit(row_y[r], "in"),
               just = c("left", "center"), gp = leg_font)
    xcur <- xcur + entry_w[k] + entry_gap
  }
}

# ---- title ----------------------------------------------------------------------------------
grid.text("Canonical factor loadings by feature group", x = unit(FW / 2, "in"), y = unit(FH - 0.20, "in"),
           gp = gpar(fontsize = 11.5, fontface = "bold", fontfamily = "serif"))

dev.off()
message("Stage 47 complete: ", out_path)
