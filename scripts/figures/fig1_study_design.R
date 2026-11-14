#!/usr/bin/env Rscript
# ============================================================
#  Figure 1 — study design
#  Card layout and typography after the reference figure.
# ============================================================
suppressMessages({library(grid)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

INK  <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
POS  <- "#8E6FA8"; NEG  <- "#5B8C8A"; AMB <- "#B08A5E"; BLU <- "#6E8CA8"
C1 <- list(bg = "#FFFFFF", bar = INK3, sum = "#EDEFF1")
C2 <- list(bg = "#F4F8FB", bar = BLU,  sum = "#DCE8F2")
C3 <- list(bg = "#F8F5FA", bar = POS,  sum = "#E7DDF0")
C4 <- list(bg = "#F2F8F7", bar = NEG,  sum = "#D7EAE7")
C5 <- list(bg = "#F4F8FB", bar = BLU,  sum = "#DCE8F2")   # conclusion layer
C6 <- list(bg = "#F6F7F8", bar = INK3, sum = "#E8EAEC")   # open questions: neutral grey
BRD <- "#C9CDD1"

W <- 170/25.4; Hh <- 192.37/25.4    # 170 x 192.4 mm
OUT <- file.path(ROOT, "Fig1_study_design.pdf")

# design coordinates (0.150-0.978) -> canvas coordinates (0.025-0.985)
CY <- function(y) 0.025 + (y - 0.150) * (0.960 / 0.828)

card <- function(ybot, ytop, xl, xr, cc) {
  yb <- CY(ybot); yt <- CY(ytop)
  grid.roundrect(unit((xl + xr) / 2, "npc"), unit((yb + yt) / 2, "npc"),
                 width = unit(xr - xl, "npc"), height = unit(yt - yb, "npc"),
                 r = unit(1.2, "mm"), gp = gpar(fill = cc$bg, col = BRD, lwd = 0.7))
  grid.rect(unit(xl + 0.004, "npc"), unit((yb + yt) / 2, "npc"),
            width = unit(0.008, "npc"), height = unit(yt - yb - 0.004, "npc"),
            gp = gpar(fill = cc$bar, col = NA))
}

layer <- function(ybot, ytop, num, title, detail, stats, summary, cc) {
  card(ybot, ytop, 0.018, 0.982, cc)
  hy <- ytop - 0.024
  grid.text(paste("LAYER", num), unit(0.032, "npc"), unit(CY(hy), "npc"),
            just = c("left", "centre"), gp = gpar(fontsize = 8.81, col = INK, fontface = "bold"))
  grid.text(title, unit(0.120, "npc"), unit(CY(hy), "npc"),
            just = c("left", "centre"), gp = gpar(fontsize = 9.16, col = INK, fontface = "bold"))
  for (i in seq_along(detail))
    grid.text(detail[i], unit(0.032, "npc"), unit(CY(hy - 0.021 - (i - 1) * 0.017), "npc"),
              just = c("left", "centre"), gp = gpar(fontsize = 6.69, col = INK2))
  for (i in seq_along(stats$n)) {
    yy <- hy - 0.019 * (i - 1)
    grid.text(stats$n[i], unit(0.690, "npc"), unit(CY(yy), "npc"), just = c("right", "centre"),
              gp = gpar(fontsize = 10.74, col = INK, fontface = "bold"))
    grid.text(stats$lab[i], unit(0.720, "npc"), unit(CY(yy), "npc"), just = c("left", "centre"),
              gp = gpar(fontsize = 7.31, col = INK3))
  }
  sy <- ybot + 0.020
  grid.rect(unit(0.5, "npc"), unit(CY(sy), "npc"), width = unit(0.948, "npc"),
            height = unit(CY(sy + 0.014) - CY(sy - 0.014), "npc"),
            gp = gpar(fill = cc$sum, col = NA))
  grid.text(summary, unit(0.030, "npc"), unit(CY(sy), "npc"), just = c("left", "centre"),
            gp = gpar(fontsize = 7.05, col = INK, fontface = "bold"))
}

down_arrow <- function(y1, y2) {
  grid.lines(unit(c(0.5, 0.5), "npc"), unit(CY(c(y1, y2)), "npc"),
             arrow = grid::arrow(length = unit(1.5, "mm"), type = "open"),
             gp = gpar(col = INK3, lwd = 0.8))
}

cairo_pdf(OUT, width = W, height = Hh, family = "Arial")
grid.newpage()

grid.text("Reproductive disorders and breast disease converge on XBP1 in immune cells",
          unit(0.5, "npc"), unit(0.985, "npc"), gp = gpar(fontsize = 11.10, col = INK, fontface = "bold"))
grid.text("Forty female reproductive traits across fourteen breast outcomes; five nested evidence layers",
          unit(0.5, "npc"), unit(0.964, "npc"), gp = gpar(fontsize = 7.40, col = INK3))

# ---------- top two boxes ----------
# input boxes: purple = exposure, teal = outcome
for (spec in list(list(xl = 0.021, bg = "#F7F3FA", bd = "#CBB9DB", bar = POS),
                  list(xl = 0.509, bg = "#F1F7F6", bd = "#AECCC8", bar = NEG))) {
  yb <- CY(0.884); yt <- CY(0.948)
  grid.roundrect(unit(spec$xl + 0.235, "npc"), unit((yb + yt) / 2, "npc"),
                 width = unit(0.470, "npc"), height = unit(yt - yb, "npc"),
                 r = unit(1.2, "mm"), gp = gpar(fill = spec$bg, col = spec$bd, lwd = 0.8))
  grid.rect(unit(spec$xl + 0.004, "npc"), unit((yb + yt) / 2, "npc"),
            width = unit(0.008, "npc"), height = unit(yt - yb - 0.004, "npc"),
            gp = gpar(fill = spec$bar, col = NA))
}
grid.text("EXPOSURE", unit(0.037, "npc"), unit(CY(0.930), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 8.63, col = INK, fontface = "bold"))
grid.text("FinnGen reproductive traits", unit(0.146, "npc"), unit(CY(0.930), "npc"),
          just = c("left", "centre"), gp = gpar(fontsize = 7.93, col = INK2))
grid.text("40 traits; up to 232,058 women; 36 with estimable rg",
          unit(0.037, "npc"), unit(CY(0.910), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 6.69, col = INK2))
grid.text("18 instrumentable, 162 independent instruments",
          unit(0.037, "npc"), unit(CY(0.893), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 6.69, col = INK3))

grid.text("OUTCOMES", unit(0.525, "npc"), unit(CY(0.930), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 8.63, col = INK, fontface = "bold"))
grid.text("Breast phenotypes", unit(0.664, "npc"), unit(CY(0.930), "npc"),
          just = c("left", "centre"), gp = gpar(fontsize = 7.93, col = INK2))
grid.text("14 outcomes: BCAC 2020 subtypes, FinnGen benign, density",
          unit(0.525, "npc"), unit(CY(0.910), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 6.69, col = INK2))
grid.text("5 malignant, 3 density, 6 benign or other",
          unit(0.525, "npc"), unit(CY(0.893), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 6.69, col = INK3))

# ---------- evidence layers ----------
down_arrow(0.878, 0.870)
layer(0.760, 0.868, "1", "Global genetic architecture",
      c("LD score regression across every admissible exposure-outcome pair;",
        "Benjamini-Hochberg correction within the pair family"),
      list(n = c("560", "483", "7"), lab = c("pairs crossed", "admissible", "survive FDR")),
      "Broad but weak  |  7 pairs survive correction, 5 of them benign rather than malignant", C1)

down_arrow(0.752, 0.744)
layer(0.636, 0.744, "2", "Local convergence",
      c("LAVA across 737 linkage-disequilibrium blocks, 1,528 block-pair tests;",
        "two exposures meet in one block and a third reaches nominal support"),
      list(n = c("1,528", "2", "1"), lab = c("block-pair tests", "exposures converge", "shared block")),
      "Ovarian cyst (rho = 0.753) and menorrhagia (rho = 0.870) at chr22:27.95-29.46 Mb", C2)

down_arrow(0.628, 0.620)
layer(0.512, 0.620, "3", "Molecular localisation",
      c("Two-sample MR; colocalisation against eQTL from four independent resources;",
        "the same variant interrogated in immune cells and in breast epithelium"),
      list(n = c("182", "12/27", "1.000"), lab = c("MR combinations", "immune cell types", "PP.H4 immune")),
      "rs2881496 acts on XBP1 in lymphocytes  |  breast PP.H4 = 0.061, no colocalisation", C3)

down_arrow(0.504, 0.496)
layer(0.388, 0.496, "4", "Exclusion of alternative routes",
      c("Four competing explanations tested at the same variant; circulating mediators",
        "screened in seven panels with clumped instruments, corrected within panel"),
      list(n = c("324", "11", "0"), lab = c("mediators screened", "breast eQTL genes", "survivors")),
      "Epithelium, density, chromatin and 324 circulating mediators all excluded", C4)

# ---------- layer five and open questions ----------
down_arrow(0.380, 0.372)

card(0.150, 0.372, 0.018, 0.640, C5)
grid.text("LAYER 5", unit(0.032, "npc"), unit(CY(0.348), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 8.81, col = INK, fontface = "bold"))
grid.text("Retained model", unit(0.120, "npc"), unit(CY(0.348), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 9.16, col = INK, fontface = "bold"))
for (i in seq_len(5))
  grid.text(c("One shared variant, rs2881496 (chr22:28,921,265), acting",
              "on XBP1 expression in lymphocytes; effect confined to",
              "oestrogen-receptor-positive disease and concordant",
              "across four independent expression resources",
              "(eQTLGen, GTEx breast, TenK10K single-cell and DICE)")[i],
            unit(0.032, "npc"), unit(CY(0.324 - (i - 1) * 0.026), "npc"),
            just = c("left", "centre"), gp = gpar(fontsize = 7.57, col = INK2))
sy <- 0.170
grid.rect(unit((0.018 + 0.640) / 2, "npc"), unit(CY(sy), "npc"), width = unit(0.606, "npc"),
          height = unit(CY(sy + 0.014) - CY(sy - 0.014), "npc"), gp = gpar(fill = C5$sum, col = NA))
grid.text("One variant, one gene, one compartment", unit(0.030, "npc"), unit(CY(sy), "npc"),
          just = c("left", "centre"), gp = gpar(fontsize = 7.05, col = INK, fontface = "bold"))

card(0.150, 0.372, 0.662, 0.982, C6)
grid.text("OPEN", unit(0.676, "npc"), unit(CY(0.348), "npc"), just = c("left", "centre"),
          gp = gpar(fontsize = 8.81, col = INK, fontface = "bold"))
grid.text("What would settle it", unit(0.760, "npc"), unit(CY(0.348), "npc"),
          just = c("left", "centre"), gp = gpar(fontsize = 9.16, col = INK, fontface = "bold"))
for (i in seq_len(5))
  grid.text(c("Perturb XBP1 in primary lymphocytes,",
              "measure spliced XBP1s, not total XBP1,",
              "and use larger female-specific",
              "reproductive GWAS, where instrument",
              "count is the constraint on precision")[i],
            unit(0.676, "npc"), unit(CY(0.324 - (i - 1) * 0.026), "npc"),
            just = c("left", "centre"), gp = gpar(fontsize = 7.57, col = INK2))
grid.rect(unit((0.662 + 0.982) / 2, "npc"), unit(CY(sy), "npc"), width = unit(0.312, "npc"),
          height = unit(CY(sy + 0.014) - CY(sy - 0.014), "npc"), gp = gpar(fill = C6$sum, col = NA))
grid.text("No functional experiment performed", unit(0.674, "npc"), unit(CY(sy), "npc"),
          just = c("left", "centre"), gp = gpar(fontsize = 7.05, col = INK, fontface = "bold"))

dev.off()
cat("done:", OUT, "\n")
