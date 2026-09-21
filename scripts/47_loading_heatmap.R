# Stage 47 -- Grouped canonical-loading heatmap (Results-chapter rewrite, §6 revision pass).
#
# One panel per feature group, 3 columns (PC1/PC2/PC3) x n-feature rows, square cells,
# posterior-mean loading on a diverging scale, a solid group-colour strip in place of
# per-row colour, no in-cell numbers, one shared colour bar, a titled legend below.
#
# Group scheme (revision pass, 2026-09): the REAL 8-group feature_group_lookup()
# (src/loading_visualization.R), restricted to the 48 features actually in the model --
# "Shooting rates" still has zero members among the 48 real features (both its features
# were dropped in the 53->48 rebuild) and is excluded here as before. "Passing volume" was
# split this revision into "Passing volume" (bulk/directional passing) and "Chance creation"
# (through balls, smart passes, crosses, key passes, shot assists -- the creative/attacking-
# associated passes, all PC1-negative) so no single panel dominates the figure's height and
# season-trend facets (stage 33) get a clean 2x4 grid instead of one empty cell.
#
# Layout: THREE panel columns (was two), since 8 groups in 2 columns would leave one column
# very tall. Column split chosen to balance row counts as evenly as an 8-group/3-column split
# allows, grouped thematically rather than strictly by PC1 sign: Discipline+Defending (15
# rows, "defensive/discipline") | Attacking output+Chance creation+Progressive actions (16
# rows, "attacking") | Passing volume+Passing quality+Dueling (17 rows, "passing/dueling") --
# ordered shortest-to-longest left to right. Each column carries its own PC1/PC2/PC3 header
# (previously only the first column did, since there was one shared visual block); headers
# get their own dedicated vertical band above the panels rather than sharing space with them.
#
# Colour palette: group_colours() was redesigned this revision to avoid red and blue
# entirely (both are used by this heatmap's own diverging loading scale below), and to keep
# every colour mutually distinct (an earlier palette had three different greens).
#
# Reliability: no longer marked on the heatmap itself (a per-cell asterisk was tried and
# dropped as not adding enough to be worth the clutter) -- the 90% CI bounds are still in
# outputs/tables/31_pca_loading_ci.csv (and Table 5 in results.tex) for anyone who needs them.
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

# ---- group order (real 8-group scheme, thematic split, balanced to 15/16/17 rows,
# ---- shortest column first) ---------------------------------------------------------------
GROUP_ORDER_COL1 <- c("Discipline", "Defending")
GROUP_ORDER_COL2 <- c("Attacking output", "Chance creation", "Progressive actions")
GROUP_ORDER_COL3 <- c("Passing volume", "Passing quality", "Dueling")
ALL_GROUPS <- c(GROUP_ORDER_COL1, GROUP_ORDER_COL2, GROUP_ORDER_COL3)

gcols <- group_colours()
stopifnot(all(ALL_GROUPS %in% names(gcols)))

clean_feature <- function(x) gsub("_", " ", gsub("^per90_|^rate_", "", x))

# Fixed +/-1 colour scale (not data-driven) per user request: an interpretable, stable range
# for a correlation-like loading rather than one that rescales itself around whatever the
# current real max happens to be.
vmax <- 1.0
message(sprintf("Colour scale fixed at +/-%.2f (actual max |loading| in data: %.4f, feature: %s)",
                 vmax, max(abs(d$loading)), d$feature[which.max(abs(d$loading))]))

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
n_col3 <- sum(vapply(group_data[GROUP_ORDER_COL3], `[[`, integer(1), "n"))
message(sprintf("Column split: col1 = %d rows (%s), col2 = %d rows (%s), col3 = %d rows (%s)",
                 n_col1, paste(GROUP_ORDER_COL1, collapse = "+"),
                 n_col2, paste(GROUP_ORDER_COL2, collapse = "+"),
                 n_col3, paste(GROUP_ORDER_COL3, collapse = "+")))

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
# Three label zones (one per column) are the dominant width cost of moving from 2 to 3
# columns -- each of the 8 panels can carry a feature as long as "dangerous opponent half
# recoveries", so label_w can't shrink much without clipping. Font sizes trimmed slightly
# (label 6.8->6.3pt, legend 7.6->7.3pt) to claw back some width; the canvas grows from the
# previous 7.3in to 8.4in to fit a third label+strip+heatmap block, and the LaTeX
# \includegraphics width (results.tex) is bumped from 0.93 to 0.97\textwidth to keep the
# on-page font size close to the previous revision's.
cell <- 0.22            # square cell, in
gap  <- 0.08            # gap between stacked panels within a column, in
top_margin <- 0.15      # whitespace above the header band
header_h   <- 0.24      # dedicated header band per column (PC1/PC2/PC3 labels) -- no
                         # longer shares space with the first panel, fixing the collision
strip_w <- 0.09         # colour-strip width, in
label_gap <- 0.055      # gap between labels and strip
label_w <- 1.45         # feature-label column width, in
label_fontsize  <- 6.3
header_fontsize <- 7.3

# Legend box geometry (bottom-left corner, boxed, 2 columns x 4 rows -- see the drawing code
# below). Only the fixed constants that determine its HEIGHT are needed this early: bottom_
# margin has to reserve enough room for the box before FH is fixed and the device opened;
# the box's WIDTH depends on measured label text and is resolved later, once the device is open.
leg_box_pad   <- 0.09
leg_title_h   <- 0.16
leg_title_gap <- 0.06
leg_row_h     <- 0.17
leg_row_gap   <- 0.02
leg_box_h     <- 2 * leg_box_pad + leg_title_h + leg_title_gap + 4 * leg_row_h + 3 * leg_row_gap
leg_box_clearance <- 0.14  # buffer above and below the box within bottom_margin

bottom_margin <- leg_box_h + 2 * leg_box_clearance

col_h <- function(groups) {
  n <- sum(vapply(group_data[groups], `[[`, integer(1), "n"))
  n * cell + (length(groups) - 1) * gap
}
col1_h <- col_h(GROUP_ORDER_COL1)
col2_h <- col_h(GROUP_ORDER_COL2)
col3_h <- col_h(GROUP_ORDER_COL3)
max_col_h <- max(col1_h, col2_h, col3_h)

FH <- top_margin + header_h + max_col_h + bottom_margin

# Horizontal layout: label zone | strip | heatmap(3 cells), for each of 3 columns. No more
# right-hand colour-bar zone needed -- the colour bar moved below columns 2-3 (horizontal,
# drawn later) instead of occupying a fourth vertical zone to the right. inter_col_gap
# trimmed (was 0.30) per user feedback that the space between columns looked a little wide.
hm_w <- 3 * cell
col_block_w <- label_w + label_gap + strip_w + hm_w
inter_col_gap <- 0.22
col1_x0 <- 0.12
col2_x0 <- col1_x0 + col_block_w + inter_col_gap
col3_x0 <- col2_x0 + col_block_w + inter_col_gap

FW <- col3_x0 + col_block_w + 0.15  # just enough right-hand breathing room, no bar zone

out_path <- file.path(paths$figures, "47_loading_heatmap.png")
png(out_path, width = FW, height = FH, units = "in", res = 220)
grid.newpage()

draw_panel <- function(gname, x0, y_top) {
  gd <- group_data[[gname]]
  n  <- gd$n
  h  <- n * cell
  y0 <- y_top - h  # bottom of this panel, inches from bottom

  # colour strip -- black border added around the swatch so every group colour reads as a
  # bordered block, same treatment as the legend swatches below
  strip_vp <- viewport(x = unit(x0 + label_w + label_gap, "in"), y = unit(y0, "in"),
                        width = unit(strip_w, "in"), height = unit(h, "in"),
                        just = c("left", "bottom"))
  pushViewport(strip_vp)
  grid.rect(gp = gpar(fill = gcols[[gname]], col = "black", lwd = 0.6))
  popViewport()

  # heatmap cells -- contiguous (no gap), separated only by a thin white border between
  # cells, and one black border wraps the whole n x 3 block once, below.
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
      grid.rect(gp = gpar(fill = col, col = "white", lwd = 1.1))
      popViewport()
    }
    # feature label, right-aligned, ending just left of the strip
    lab_vp <- viewport(x = unit(x0, "in"), y = unit(row_y0, "in"),
                        width = unit(label_w, "in"), height = unit(cell, "in"),
                        just = c("left", "bottom"))
    pushViewport(lab_vp)
    grid.text(gd$labels[i], x = unit(1, "npc"), y = unit(0.5, "npc"),
               just = c("right", "center"), gp = gpar(fontsize = label_fontsize, fontfamily = "serif"))
    popViewport()
  }

  # one black border around the whole group's heatmap block (replaces per-cell borders as
  # the "this is one group's block" cue)
  block_vp <- viewport(x = unit(hm_x0, "in"), y = unit(y0, "in"),
                        width = unit(3 * cell, "in"), height = unit(h, "in"), just = c("left", "bottom"))
  pushViewport(block_vp)
  grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1.0))
  popViewport()

  invisible(y0)
}

draw_column <- function(groups, x0, panel_top) {
  # dedicated header band, drawn once per column, directly above that column's own panels
  hm_x0 <- x0 + label_w + label_gap + strip_w
  for (j in 1:3) {
    vp <- viewport(x = unit(hm_x0 + (j - 1) * cell, "in"), y = unit(panel_top, "in"),
                    width = unit(cell, "in"), height = unit(header_h, "in"),
                    just = c("left", "bottom"))
    pushViewport(vp)
    grid.text(paste0("PC", j), y = unit(0.5, "npc"),
               gp = gpar(fontsize = header_fontsize, fontfamily = "serif", fontface = "bold"))
    popViewport()
  }

  y_cursor <- panel_top
  for (k in seq_along(groups)) {
    gname <- groups[k]
    y0 <- draw_panel(gname, x0, y_cursor)
    y_cursor <- y0 - gap
  }
}

panel_top <- FH - top_margin - header_h
draw_column(GROUP_ORDER_COL1, col1_x0, panel_top)
draw_column(GROUP_ORDER_COL2, col2_x0, panel_top)
draw_column(GROUP_ORDER_COL3, col3_x0, panel_top)

# ---- shared colour bar: small, horizontal, centred below columns 2 and 3 -----------------
# Moved off the right margin (was a tall vertical bar) and shrunk twice: first from a full
# right-hand vertical bar to a horizontal one spanning columns 2-3, then to this small fixed-
# width bar centred under that same span, per user feedback that both earlier versions were
# too large. Sits in the same bottom band as the legend box (below column 1).
cbar_bar_h    <- 0.13
cbar_title_h  <- 0.16
cbar_title_gap<- 0.05
cbar_label_gap<- 0.04
cbar_label_h  <- 0.13

cbar_span   <- 1.9  # fixed width -- bigger than the first small pass, still well short of
                     # the full columns 2-3 span it used to stretch across
cbar_centre <- (col2_x0 + col3_x0 + col_block_w) / 2
cbar_x0     <- cbar_centre - cbar_span / 2
cbar_x1     <- cbar_x0 + cbar_span

cbar_bar_y0 <- leg_box_clearance + cbar_label_h + cbar_label_gap
cbar_bar_top<- cbar_bar_y0 + cbar_bar_h
cbar_title_y<- cbar_bar_top + cbar_title_gap + cbar_title_h / 2
cbar_label_y<- leg_box_clearance + cbar_label_h / 2

n_grad <- 256
grad_cols <- loading_to_colour(seq(-vmax, vmax, length.out = n_grad), vmax)
cbar_vp <- viewport(x = unit(cbar_x0, "in"), y = unit(cbar_bar_y0, "in"),
                     width = unit(cbar_span, "in"), height = unit(cbar_bar_h, "in"), just = c("left", "bottom"))
pushViewport(cbar_vp)
grid.raster(matrix(grad_cols, nrow = 1), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
grid.rect(gp = gpar(fill = NA, col = "grey40", lwd = 0.4))
popViewport()

grid.text("posterior-mean loading", x = unit(cbar_x0 + cbar_span / 2, "in"), y = unit(cbar_title_y, "in"),
           just = c("center", "center"), gp = gpar(fontsize = 6.3, fontfamily = "serif"))
grid.text(sprintf("%.0f", -vmax), x = unit(cbar_x0, "in"), y = unit(cbar_label_y, "in"),
           just = c("center", "center"), gp = gpar(fontsize = 5.8, fontfamily = "serif"))
grid.text("0", x = unit(cbar_x0 + cbar_span / 2, "in"), y = unit(cbar_label_y, "in"),
           just = c("center", "center"), gp = gpar(fontsize = 5.8, fontfamily = "serif"))
grid.text(sprintf("+%.0f", vmax), x = unit(cbar_x1, "in"), y = unit(cbar_label_y, "in"),
           just = c("center", "center"), gp = gpar(fontsize = 5.8, fontfamily = "serif"))

# ---- legend (boxed, bottom-left corner, entries aligned in a 2-column x 4-row grid) --------
# Moved off centre-bottom and into an enclosed box anchored under column 1, per user request;
# entries line up in two aligned sub-columns rather than the previous variable-width row flow.
leg_title_font <- gpar(fontsize = 8.0, fontfamily = "serif", fontface = "bold")
leg_font <- gpar(fontsize = 7.0, fontfamily = "serif")
swatch_w <- 0.13; swatch_gap <- 0.06; subcol_gap <- 0.22

leg_pairs <- list(ALL_GROUPS[1:2], ALL_GROUPS[3:4], ALL_GROUPS[5:6], ALL_GROUPS[7:8])
subcol1_groups <- vapply(leg_pairs, `[[`, character(1), 1)
subcol2_groups <- vapply(leg_pairs, `[[`, character(1), 2)
label_w1 <- max(vapply(subcol1_groups, function(g)
  convertWidth(grobWidth(textGrob(g, gp = leg_font)), "in", valueOnly = TRUE), numeric(1)))
label_w2 <- max(vapply(subcol2_groups, function(g)
  convertWidth(grobWidth(textGrob(g, gp = leg_font)), "in", valueOnly = TRUE), numeric(1)))
subcol1_w <- swatch_w + swatch_gap + label_w1
subcol2_w <- swatch_w + swatch_gap + label_w2
box_w <- 2 * leg_box_pad + subcol1_w + subcol_gap + subcol2_w

box_x0 <- col1_x0 + 0.4  # nudged right of column 1's own left edge, per user feedback
box_y0 <- leg_box_clearance
box_top <- box_y0 + leg_box_h
stopifnot(box_top + leg_box_clearance <= bottom_margin + 1e-9)  # must clear the panels above

grid.rect(x = unit(box_x0, "in"), y = unit(box_y0, "in"), width = unit(box_w, "in"), height = unit(leg_box_h, "in"),
          just = c("left", "bottom"), gp = gpar(fill = "white", col = "grey40", lwd = 0.5))

grid.text("Feature group", x = unit(box_x0 + box_w / 2, "in"),
           y = unit(box_top - leg_box_pad - leg_title_h / 2, "in"),
           just = c("center", "center"), gp = leg_title_font)

subcol1_x0 <- box_x0 + leg_box_pad
subcol2_x0 <- subcol1_x0 + subcol1_w + subcol_gap
row_top <- box_top - leg_box_pad - leg_title_h - leg_title_gap
for (r in seq_along(leg_pairs)) {
  ry <- row_top - (r - 1) * (leg_row_h + leg_row_gap) - leg_row_h / 2
  pair <- leg_pairs[[r]]
  for (side in 1:2) {
    gx0 <- if (side == 1) subcol1_x0 else subcol2_x0
    gname <- pair[side]
    sw_vp <- viewport(x = unit(gx0, "in"), y = unit(ry, "in"),
                        width = unit(swatch_w, "in"), height = unit(swatch_w, "in"), just = c("left", "center"))
    pushViewport(sw_vp)
    grid.rect(gp = gpar(fill = gcols[[gname]], col = "black", lwd = 0.6))
    popViewport()
    grid.text(gname, x = unit(gx0 + swatch_w + swatch_gap, "in"), y = unit(ry, "in"),
               just = c("left", "center"), gp = leg_font)
  }
}

# Title deliberately omitted (G2): the LaTeX caption and section heading already name this
# figure ("Canonical factor loadings by feature group") -- a burnt-in title would be redundant
# and taller-than-necessary float boxes are exactly what G1 is fixing.

dev.off()
message("Stage 47 complete: ", out_path)
