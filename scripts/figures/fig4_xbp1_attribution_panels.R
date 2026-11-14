#!/usr/bin/env Rscript
# ============================================================
#  Figure 4 — XBP1 attribution: immune colocalisation vs breast exclusion
# ============================================================
suppressMessages({library(data.table); library(grDevices)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

COL <- fread(file.path(ROOT, "_fig4data/t10k_coloc.csv"))
BR  <- fread(file.path(ROOT, "repro_breast_L3/sc_eqtl/coloc_GTExV8_breast.csv"))

INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
PUR <- "#8E6FA8"; TEAL <- "#3E7D5A"

# ---------- XBP1 cell-type coloc ----------
X <- COL[Gene == "XBP1"]
setnames(X, c("V1","SNP","PP.H4","Celltype","Gene"), c("id","SNP","H4","Celltype","Gene"), skip_absent = TRUE)
if (!"H4" %in% names(X)) setnames(X, "PP.H4", "H4")
X <- X[order(-H4)]

lineage <- function(ct) {
  if (grepl("^CD4", ct)) "CD4 T"
  else if (grepl("^CD8", ct)) "CD8 T"
  else if (grepl("^(MAIT|gdT|dnT|Treg|ILC|NK)", ct)) "NK / innate-like T"
  else if (grepl("^B_|Plasmablast", ct)) "B"
  else if (grepl("Mono", ct)) "Monocyte"
  else if (grepl("^(cDC|pDC|ASDC)", ct)) "DC"
  else "Other"
}
X[, lin := sapply(Celltype, lineage)]
LINCOL <- c("CD4 T"="#8E6FA8", "CD8 T"="#B5657E", "NK / innate-like T"="#5B8C8A",
            "B"="#4A7C9B", "Monocyte"="#B08A5E", "DC"="#7A6E9B", "Other"="#9AA0A6")

# ==================== Fig 4a — PP.H4 across 27 immune cell types ====================
cairo_pdf(file.path(ROOT, "Fig4a_immune_coloc.pdf"), width = 8.8, height = 8.6, family = "Arial")

par(fig = c(0.285, 0.905, 0.110, 0.762), mar = c(0,0,0,0))
n <- nrow(X); yy <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(0, 1.06), ylim = c(0.4, n + 0.6))
for (th in c(0.5, 0.8)) abline(v = th, lty = 2, col = "#E3E6EA", lwd = 0.8)
for (i in seq_len(n)) {
  y <- yy[i]
  segments(0, y, X$H4[i], y, col = "#DDE1E6", lwd = 1.3)
  points(X$H4[i], y, pch = 21, cex = 1.20, bg = LINCOL[X$lin[i]], col = "#FFFFFF", lwd = 0.4)
}
axis(2, at = yy, labels = X$Celltype, cex.axis = 0.58, tick = FALSE, col.axis = INK2, line = -1.0, las = 1)
axis(1, cex.axis = 0.66, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6, las = 1)
mtext("Posterior probability of a shared causal variant (PP.H4)", side = 1, line = 2.4,
      cex = 0.68, col = INK2, las = 0)
text(0.515, n + 0.40, "0.5", cex = 0.46, col = INK3)
text(0.815, n + 0.40, "0.8", cex = 0.46, col = INK3)

# inset: PP.H4 distribution
par(fig = c(0.285, 0.640, 0.790, 0.905), mar = c(1.9,2.6,1.8,0.4), new = TRUE, mgp = c(1.2,0.35,0))
hh <- hist(X$H4, breaks = seq(0, 1, by = 0.1), plot = FALSE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0, max(hh$counts)*1.15))
for (k in seq_along(hh$counts)) {
  if (hh$counts[k] <= 0) next
  rect(hh$breaks[k], 0, hh$breaks[k+1], hh$counts[k], col = "#C9BFD8", border = "white", lwd = 0.5)
}
axis(1, cex.axis = 0.44, tck = -0.03, col.axis = INK3, lwd = 0.5, las = 1)
axis(2, cex.axis = 0.44, tck = -0.03, col.axis = INK3, lwd = 0.5, las = 1, line = -0.5)
mtext("Cell types", side = 2, line = 1.0, cex = 0.42, col = INK3, las = 0)
mtext("PP.H4 distribution", side = 3, line = 0.35, cex = 0.46, col = INK2)
text(0.50, max(hh$counts)*0.92, sprintf("median PP.H4 = %.2f", median(X$H4)),
     adj = 0.5, cex = 0.44, col = INK3)

# legend and title
par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
lins <- names(LINCOL)
for (i in seq_along(lins)) {
  x <- 0.030 + ((i-1) %% 4) * 0.155
  y <- 0.062 - ((i-1) %/% 4) * 0.026
  points(x, y, pch = 21, cex = 0.80, bg = LINCOL[lins[i]], col = "#FFFFFF", lwd = 0.4)
  text(x + 0.011, y, lins[i], adj = 0, cex = 0.48, col = INK2)
}
text(0.650, 0.115, sprintf("%d / 27 at PP.H4 = 1.000", sum(X$H4 > 0.9999)),
     adj = 0, cex = 0.54, col = INK2, font = 2)
text(0.650, 0.086, sprintf("%d above 0.8    %d above 0.5", sum(X$H4 > 0.8), sum(X$H4 > 0.5)),
     adj = 0, cex = 0.50, col = INK3)
text(0.5, 0.975, "XBP1 expression colocalises with breast cancer risk in immune cells",
     adj = 0.5, cex = 0.84, col = INK2)
text(0.5, 0.946, "TenK10K single-cell eQTL across 27 immune cell types, colocalised against BCAC luminal A",
     adj = 0.5, cex = 0.52, col = INK3)
dev.off(); cat("Fig4a done
")

# ==================== Fig 4b — MR forest ====================
MR <- data.table(
  res = c("eQTLGen (whole blood)","TenK10K CD4_CTL","TenK10K CD8_TCM","TenK10K gdT",
          "TenK10K NK","TenK10K CD14_Mono","DynCD4 Naive 0h","DynCD4 TN 0h",
          "TenK10K CD4_TEM","TenK10K Treg"),
  b   = c(0.03765, 0.0928, 0.0923, 0.0921, 0.0869, 0.0571, 0.0153, 0.0156, -0.106, -0.0146),
  se  = c(0.00418, 0.0214, 0.0213, 0.0213, 0.0214, 0.0299, 0.00698, 0.00712, 0.149, 0.0737),
  p   = c(2.12e-19, 1.47e-5, 1.49e-5, 1.49e-5, 5.06e-5, 0.0563, 0.0287, 0.0287, 0.476, 0.843),
  niv = c("80 SNP","1 SNP","1 SNP","1 SNP","1 SNP","1 SNP","1 SNP","1 SNP","1 SNP","2 SNP")
)
MR[, sig := p < 0.05]
MR <- MR[order(b)]

cairo_pdf(file.path(ROOT, "Fig4b_mr_forest.pdf"), width = 8.2, height = 5.6, family = "Arial")
par(fig = c(0.300, 0.860, 0.135, 0.830), mar = c(0,0,0,0))
n <- nrow(MR); yy <- seq_len(n)
plot.new(); plot.window(xlim = c(-0.32, 0.30), ylim = c(0.4, n + 0.6))
abline(v = 0, col = "#D8DCE1", lwd = 0.9)
for (i in seq_len(n)) {
  lo <- MR$b[i] - 1.96*MR$se[i]; hi <- MR$b[i] + 1.96*MR$se[i]
  segments(lo, yy[i], hi, yy[i], col = "#C9CFD6", lwd = 1.5)
  points(MR$b[i], yy[i], pch = 21, cex = if (MR$sig[i]) 1.10 else 0.85,
         bg = if (MR$sig[i]) PUR else "#EDEDEF", col = if (MR$sig[i]) INK else "#B8BCC2",
         lwd = if (MR$sig[i]) 0.9 else 0.5)
}
axis(2, at = yy, labels = MR$res, cex.axis = 0.56, tick = FALSE, col.axis = INK2, line = -0.9, las = 1)
axis(1, cex.axis = 0.62, col.axis = INK2, tck = -0.013, col = INK3, lwd = 0.6)
mtext("MR effect on luminal A breast cancer per SD of XBP1 expression (log OR)", side = 1, line = 2.3,
      cex = 0.62, col = INK2, las = 0)
for (i in seq_len(n))
  text(0.315, yy[i], ifelse(MR$p[i] < 1e-4, formatC(MR$p[i], format="e", digits=1),
                            formatC(MR$p[i], format="g", digits=2)),
       adj = 0, cex = 0.50, col = if (MR$sig[i]) INK else INK3, xpd = NA)
text(0.315, n + 0.45, "P", adj = 0, cex = 0.55, font = 2, col = INK2, xpd = NA)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "Higher XBP1 expression increases risk across T, NK and monocyte subsets",
     adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.935, "Two-sample MR; eQTLGen uses 80 instruments, single-cell resources use one lead variant each",
     adj = 0.5, cex = 0.52, col = INK3)
points(0.038, 0.045, pch = 21, cex = 1.00, bg = PUR, col = INK, lwd = 0.9)
text(0.050, 0.045, "P < 0.05", adj = 0, cex = 0.50, col = INK2)
points(0.125, 0.045, pch = 21, cex = 1.00, bg = "#EDEDEF", col = "#B8BCC2", lwd = 0.5)
text(0.137, 0.045, "not significant", adj = 0, cex = 0.50, col = INK2)
dev.off(); cat("Fig4b done\n")

# ==================== Fig 4c — tissue contrast ====================
BRl <- BR[trait == "BCAC_LumA"][order(-PP_H4)]

cairo_pdf(file.path(ROOT, "Fig4c_tissue_contrast.pdf"), width = 8.0, height = 5.8, family = "Arial")
par(fig = c(0.180, 0.945, 0.185, 0.820), mar = c(0,0,0,0))
n <- nrow(BRl) + 1
labs <- c(BRl$gene, "XBP1 (immune)")
vals <- c(BRl$PP_H4, max(X$H4))
cols <- c(rep("#C6B7D8", nrow(BRl)), PUR)
ypos <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(0, 1.05), ylim = c(0.4, n + 0.6))
abline(h = n - nrow(BRl) + 0.5, col = "#E3E6EA", lwd = 1.2)
for (i in seq_len(n)) {
  rect(0, ypos[i]-0.34, vals[i], ypos[i]+0.34, col = cols[i], border = NA)
  text(vals[i] + 0.015, ypos[i], sprintf("%.3f", vals[i]), adj = 0, cex = 0.52, col = INK2)
}
axis(2, at = ypos, labels = labs, cex.axis = 0.56, tick = FALSE, col.axis = INK2, line = -0.9, las = 1)
axis(1, cex.axis = 0.64, col.axis = INK2, tck = -0.013, col = INK3, lwd = 0.6)
mtext("PP.H4 with luminal A breast cancer", side = 1, line = 2.3, cex = 0.68, col = INK2, las = 0)
text(0.52, n + 0.35, "Breast tissue: every candidate gene below 0.2", cex = 0.50, col = INK3, adj = 0)
text(0.62, 3.2, "Immune cells: XBP1
reaches 1.000", cex = 0.50, col = PUR, adj = 0, font = 3)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "The shared variant acts in immune cells, not in breast epithelium",
     adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.935, "GTEx v8 breast eQTL (upper bars) versus TenK10K immune eQTL (lower bar) for the same locus",
     adj = 0.5, cex = 0.52, col = INK3)
dev.off(); cat("Fig4c done\n")

# ==================== Fig 4d — four molecular layers ====================
LY <- data.table(
  layer = c("Chromatin accessibility (caQTL)", "DNA methylation (mQTL)",
            "XBP1 expression (eQTL)", "Breast expression (eQTL)"),
  coloc = c(0.572, 0.823, 1.000, 0.061),
  sig   = c(1.85e-24, 3.35e-216, 5.5e-45, 0.31),
  tissue= c("NK / CD4 T","Whole blood","Whole blood / immune","Breast"),
  note  = c("chr22-28442177","cg18940763","rs2881496","no colocalisation")
)
LY[, lab := sprintf("%.0e", sig)]

cairo_pdf(file.path(ROOT, "Fig4d_molecular_layers.pdf"), width = 8.6, height = 5.4, family = "Arial")
par(fig = c(0.320, 0.760, 0.190, 0.820), mar = c(0,0,0,0))
n <- nrow(LY); yy <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(0, 1.05), ylim = c(0.4, n + 0.6))
for (th in c(0.5, 0.8)) abline(v = th, lty = 2, col = "#E3E6EA", lwd = 0.8)
for (i in seq_len(n)) {
  v <- LY$coloc[i]
  cc <- if (v > 0.8) PUR else if (v > 0.5) "#B9A9CE" else "#DDE1E6"
  rect(0, yy[i]-0.32, v, yy[i]+0.32, col = cc, border = NA)
  text(v + 0.015, yy[i], sprintf("%.3f", v), adj = 0, cex = 0.56, col = INK2, font = 2)
}
axis(2, at = yy, labels = LY$layer, cex.axis = 0.56, tick = FALSE, col.axis = INK2, line = -0.9, las = 1)
axis(1, cex.axis = 0.64, col.axis = INK2, tck = -0.013, col = INK3, lwd = 0.6)
mtext("PP.H4 with luminal A breast cancer", side = 1, line = 2.3, cex = 0.66, col = INK2, las = 0)

# right annotation column
par(fig = c(0.770, 0.995, 0.190, 0.820), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0.4, n + 0.6))
for (i in seq_len(n)) {
  text(0.0, yy[i] + 0.13, LY$tissue[i], adj = 0, cex = 0.48, col = INK2)
  text(0.0, yy[i] - 0.13, LY$note[i], adj = 0, cex = 0.44, col = INK3, font = 3)
}
mtext("signal in", side = 3, line = 0.4, cex = 0.46, col = INK3, at = 0.30)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "Four molecular layers at the same variant, three of them colocalise",
     adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.935, "rs2881496 (chr22:28,921,265) is simultaneously the lead eQTL, mQTL and caQTL for this locus",
     adj = 0.5, cex = 0.52, col = INK3)
text(0.038, 0.128, "Breast tissue is the only layer that fails to colocalise - consistent with an immune-cell mechanism",
     adj = 0, cex = 0.50, col = TEAL, font = 3)
dev.off(); cat("Fig4d done\n")
