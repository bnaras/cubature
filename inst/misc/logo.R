## data-raw/logo.R ---------------------------------------------------
## Hex sticker for the cubature package.
##
## Design: an h-adaptive quadtree subdivision of a 2D peaked integrand,
## clipped to the hex boundary, with warm colors showing cell-averaged
## function value. This visually illustrates what hcubature actually
## does: adaptively refine cells where the integrand has most of its
## mass, coarsely sample regions where the function is near zero.
##
## To regenerate: Rscript data-raw/logo.R
##
## Output: man/figures/logo.png (pkgdown expects the logo there).

library(hexSticker)
library(ggplot2)
library(polyclip)
library(showtext)

font_add_google("Oxanium", "oxanium")
showtext_auto()

## Target: a single centered Gaussian bump. We deliberately use a
## unimodal peak so the quadtree shows a clean radial refinement
## pattern and the upper part of the hex stays dark behind the text.
target <- function(x, y) {
  exp(-30 * ((x - 0.5)^2 + (y - 0.5)^2))
}

## Recursively subdivide [0,1]^2 into cells until each cell is either
## small enough or contributes little per unit volume. This mimics the
## refinement decisions that hcubature makes internally.
subdivide_if_hot <- function(cells, depth, max_depth = 6) {
  if (depth >= max_depth) return(cells)
  out <- list()
  for (c in cells) {
    mid_x <- 0.5 * (c$x0 + c$x1); mid_y <- 0.5 * (c$y0 + c$y1)
    size  <- (c$x1 - c$x0)
    if (target(mid_x, mid_y) * size^2 > 0.0004 && size > 1/64) {
      out[[length(out) + 1]] <- list(x0 = c$x0, x1 = mid_x, y0 = c$y0, y1 = mid_y)
      out[[length(out) + 1]] <- list(x0 = mid_x, x1 = c$x1, y0 = c$y0, y1 = mid_y)
      out[[length(out) + 1]] <- list(x0 = c$x0, x1 = mid_x, y0 = mid_y, y1 = c$y1)
      out[[length(out) + 1]] <- list(x0 = mid_x, x1 = c$x1, y0 = mid_y, y1 = c$y1)
    } else {
      out[[length(out) + 1]] <- c
    }
  }
  subdivide_if_hot(out, depth + 1, max_depth)
}

cells <- subdivide_if_hot(list(list(x0 = 0, x1 = 1, y0 = 0, y1 = 1)), 0)
df <- do.call(rbind, lapply(cells, function(c) {
  cx <- 0.5*(c$x0 + c$x1); cy <- 0.5*(c$y0 + c$y1)
  data.frame(x0 = c$x0, x1 = c$x1, y0 = c$y0, y1 = c$y1,
             val = target(cx, cy))
}))

## Clip every cell against a hex polygon inscribed in the unit square.
## Cells straddling the hex edge are cropped; cells outside are
## dropped. This is what keeps the subdivision inside the sticker.
hex_cx <- 0.5; hex_cy <- 0.5; hex_r <- 0.48
ang <- seq(0, 2*pi, length.out = 7) + pi/2   # pointy-top
hex_path <- list(list(x = hex_cx + hex_r * cos(ang),
                      y = hex_cy + hex_r * sin(ang)))

cell_polys <- list(); pid <- 1
for (k in seq_len(nrow(df))) {
  rect <- list(list(
    x = c(df$x0[k], df$x1[k], df$x1[k], df$x0[k]),
    y = c(df$y0[k], df$y0[k], df$y1[k], df$y1[k])))
  clipped <- polyclip(rect, hex_path, op = "intersection")
  if (length(clipped) == 0) next
  for (piece in clipped) {
    cell_polys[[length(cell_polys) + 1]] <- data.frame(
      x = piece$x, y = piece$y, id = pid, val = df$val[k])
    pid <- pid + 1
  }
}
cells_clipped <- do.call(rbind, cell_polys)

p <- ggplot(cells_clipped) +
  geom_polygon(aes(x = x, y = y, group = id, fill = val),
               color = "#0d0013", linewidth = 0.08) +
  scale_fill_gradientn(
    colors = c("#1a0025", "#3b0f4b", "#7d1e6a", "#c04a45", "#e88b2a", "#f7d148"),
    guide = "none") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(legend.position  = "none",
        plot.background  = element_rect(fill = NA, color = NA),
        panel.background = element_rect(fill = NA, color = NA),
        plot.margin      = margin(0, 0, 0, 0))

outfile <- "man/figures/logo.png"
sticker(p,
        package = "cubature",
        s_x = 1, s_y = 0.9, s_width = 1.75, s_height = 1.75,
        p_size = 24, p_y = 1.52, p_color = "#f7d148",
        p_family = "oxanium",
        h_fill = "#0d0013", h_color = "#f7d148", h_size = 1.3,
        filename = outfile)

cat("Wrote", outfile, "\n")
