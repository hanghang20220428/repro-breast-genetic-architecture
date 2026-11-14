#!/usr/bin/env Rscript
# ============================================================
#  Figure 3 — LAVA local genetic correlation (rich version)
#  Data: lava_new_bivar_full.csv (1,528 block x pair)
#  chr22 convergence block = LOC 2473 (chr22:27,952,442-29,457,601)
# ============================================================
suppressMessages({library(data.table); library(grDevices)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

L <- fread(file.path(ROOT, "repro_breast_L3/lava_run/lava_new_bivar_full.csv"))
L <- L[!is.na(rho)]

INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
EXPCOL <- c(Menorrhagia = "#8E6FA8", Ovarian_cyst = "#B5657E",
            Endometriosis = "#5B8C8A", Pelvic_organ_prolapse = "#B08A5E")
EXPLAB <- c(Menorrhagia = "Menorrhagia", Ovarian_cyst = "Ovarian cyst",
            Endometriosis = "Endometriosis", Pelvic_organ_prolapse = "Pelvic organ prolapse")

# BCAC pairs only, one row per block x exposure
B <- L[phen1 == "BCAC_LumA" | phen2 == "BCAC_LumA"]
B[, expo := ifelse(phen2 == "BCAC_LumA", phen1, phen2)]
B <- B[expo %in% names(EXPCOL)]
B[, logp := -log10(p)]
bord <- c("Menorrhagia","Ovarian_cyst","Endometriosis","Pelvic_organ_prolapse")

# ==================== Fig 3a — Manhattan over block index ====================
cairo_pdf(file.path(ROOT, "Fig3a_lava_manhattan.pdf"), width = 11.0, height = 5.6, family = "Arial")

par(fig = c(0.060, 0.978, 0.115, 0.895), mar = c(0,0,0,0))
yl <- max(B$logp, na.rm = TRUE) * 1.10
xl <- c(min(B$LOC), max(B$LOC) + 40)
plot.new(); plot.window(xlim = xl, ylim = c(0, yl))
abline(h = -log10(0.05), lty = 2, col = "#D8DCE1", lwd = 0.8)
for (e in bord) {
  s <- B[expo == e]
  points(s$LOC, s$logp, pch = 21, cex = 0.62, bg = EXPCOL[e], col = "#FFFFFF", lwd = 0.3)
}
# highlight chr22 block
h <- B[LOC == 2473]
points(h$LOC, h$logp, pch = 21, cex = 1.55, bg = NA, col = INK, lwd = 1.1)
segments(2473, -0.3, 2473, max(h$logp) + 0.4, col = INK, lwd = 0.7, lty = 3)
text(2473 - 25, max(h$logp) - 0.15, "chr22
27.95-29.46 Mb", cex = 0.52, col = INK2, adj = 1)
axis(1, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
axis(2, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
mtext("LAVA block index (genomic order)", side = 1, line = 2.3, cex = 0.70, col = INK2, las = 0)
mtext(expression(-log[10](italic(P))~"for local rg"), side = 2, line = 2.4, cex = 0.70, col = INK2, las = 0)
text(xl[1] + 20, -log10(0.05) + yl*0.02, "P = 0.05", adj = 0, cex = 0.52, col = INK3)

# inset: rho distribution
par(fig = c(0.100, 0.285, 0.645, 0.870), mar = c(2.2,2.6,2.0,0.4), new = TRUE, mgp = c(1.4,0.4,0))
dd <- density(B$rho, from = -1, to = 1, n = 400)
plot.new(); plot.window(xlim = c(-1,1), ylim = c(0, max(dd$y)*1.15))
polygon(c(-1, dd$x, 1), c(0, dd$y, 0), col = "#DCD3E6", border = NA)
abline(v = 0, col = INK3, lwd = 0.7)
axis(1, cex.axis = 0.42, tck = -0.03, col.axis = INK3, lwd = 0.5)
axis(2, cex.axis = 0.42, tck = -0.03, col.axis = INK3, lwd = 0.5)
mtext("rho distribution", side = 3, line = 0.35, cex = 0.46, col = INK2)
mtext("local rg", side = 1, line = 1.1, cex = 0.42, col = INK3, las = 0)

# inset: significant blocks per exposure
par(fig = c(0.330, 0.535, 0.640, 0.870), mar = c(2.2,3.2,2.0,0.4), new = TRUE, mgp = c(1.4,0.4,0))
cnt <- B[, .(n = sum(p < 0.05), nf = sum(p < 0.05/nrow(B))), by = expo][match(bord, expo)]
plot.new(); plot.window(xlim = c(0, max(cnt$n)*1.25), ylim = c(0.4, 4.6))
for (i in 1:4) {
  y <- 5 - i
  rect(0, y-0.34, cnt$n[i], y+0.34, col = EXPCOL[bord[i]], border = NA)
  if (!is.na(cnt$nf[i]) && cnt$nf[i] > 0)
    rect(cnt$n[i], y-0.34, cnt$n[i] + cnt$nf[i]*0.55, y+0.34, col = INK, border = NA)
  text(cnt$n[i] + 0.35, y, cnt$n[i], cex = 0.44, col = INK2, adj = 0)
}
axis(2, at = 4:1, labels = c("Menorrhagia","Ovarian cyst","Endometriosis","Prolapse"),
     cex.axis = 0.40, tick = FALSE, col.axis = INK2, line = -0.8, las = 1)
mtext("Significant blocks", side = 3, line = 0.35, cex = 0.46, col = INK2)

# legend and title
par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
lx <- 0.090
for (i in 1:4) {
  points(lx + (i-1)*0.155, 0.055, pch = 21, cex = 0.75, bg = EXPCOL[bord[i]], col = "#FFFFFF", lwd = 0.3)
  text(lx + (i-1)*0.155 + 0.010, 0.055, EXPLAB[bord[i]], adj = 0, cex = 0.50, col = INK2)
}
points(0.700, 0.055, pch = 21, cex = 1.35, bg = NA, col = INK, lwd = 1.0)
text(0.712, 0.055, "chr22 convergence block", adj = 0, cex = 0.50, col = INK2)
text(0.5, 0.975, "Genome-wide local genetic correlation between reproductive disorders and Luminal A breast cancer",
     adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.945, sprintf("%d LD blocks tested; 4 exposures had informative blocks vs Luminal A", uniqueN(B$LOC)),
     adj = 0.5, cex = 0.52, col = INK3)
dev.off(); cat("Fig3a done\n")

# ==================== Fig 3b — chr22 block detail ====================
CH <- L[LOC %in% c(2473, 2474) & !is.na(rho)]
CH[, pairlab := gsub("_", " ", pair)]
CH[, isBC := grepl("BCAC", pair)]
CH <- CH[order(-isBC, -rho)]

cairo_pdf(file.path(ROOT, "Fig3b_chr22_detail.pdf"), width = 8.6, height = 6.6, family = "Arial")
par(fig = c(0.290, 0.985, 0.105, 0.830), mar = c(0,0,0,0))
n <- nrow(CH); yy <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(-1.05, 1.05), ylim = c(0.4, n + 1.2))
abline(v = 0, col = "#D8DCE1", lwd = 0.9)
for (i in seq_len(n)) {
  y <- yy[i]
  lo <- max(CH$rho.lower[i], -1); hi <- min(CH$rho.upper[i], 1)
  if (is.finite(lo) && is.finite(hi) && hi - lo < 1.9)
    segments(lo, y, hi, y, col = "#C9CFD6", lwd = 1.5)
  sig <- !is.na(CH$p[i]) && CH$p[i] < 0.05
  points(CH$rho[i], y, pch = 21, cex = if (sig) 1.15 else 0.85,
         bg = if (sig) if (CH$isBC[i]) "#8E6FA8" else "#C6B7D8" else "#EDEDEF",
         col = if (sig) INK else "#B8BCC2", lwd = if (sig) 0.9 else 0.5)
}
axis(2, at = yy, labels = sprintf("%s  (LOC %d)", CH$pairlab, CH$LOC),
     cex.axis = 0.52, tick = FALSE, col.axis = INK2, line = -0.9, las = 1)
axis(1, cex.axis = 0.64, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
mtext("Local genetic correlation (rho) with 95% CI", side = 1, line = 2.2, cex = 0.68, col = INK2, las = 0)
# P value column
for (i in seq_len(n))
  text(1.02, yy[i], ifelse(is.na(CH$p[i]), "NA", formatC(CH$p[i], format = "e", digits = 1)),
       adj = 0, cex = 0.50, col = if (!is.na(CH$p[i]) && CH$p[i] < 0.05) INK else INK3, xpd = NA)
text(1.02, n + 0.95, "P", adj = 0, cex = 0.56, font = 2, col = INK2, xpd = NA)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "The chr22 convergence block carries shared signal across phenotypes",
     adj = 0.5, cex = 0.82, col = INK2)
text(0.5, 0.935, "LOC 2473 (chr22:27,952,442-29,457,601) and the adjacent LOC 2474 (chr22:29,457,602-30,962,255)",
     adj = 0.5, cex = 0.53, col = INK3)
points(0.038, 0.032, pch = 21, cex = 1.10, bg = "#8E6FA8", col = INK, lwd = 0.9)
text(0.050, 0.055, "reproductive disorder x Luminal A", adj = 0, cex = 0.50, col = INK2)
points(0.300, 0.032, pch = 21, cex = 1.10, bg = "#C6B7D8", col = INK, lwd = 0.9)
text(0.312, 0.032, "reproductive disorder x reproductive disorder", adj = 0, cex = 0.50, col = INK2)
points(0.640, 0.032, pch = 21, cex = 0.90, bg = "#EDEDEF", col = "#B8BCC2", lwd = 0.5)
text(0.652, 0.032, "not significant", adj = 0, cex = 0.50, col = INK2)
dev.off(); cat("Fig3b done\n")

# ==================== Fig 3c — co-convergence scatter ====================
W <- dcast(L[LOC %in% L$LOC], LOC ~ pair, value.var = "rho")
M1 <- L[pair == "BCAC_LumA x Menorrhagia", .(LOC, r1 = rho)]
M2 <- L[pair == "BCAC_LumA x Ovarian_cyst",  .(LOC, r2 = rho)]
CS <- merge(M1, M2, by = "LOC")
CS <- CS[!is.na(r1) & !is.na(r2)]

cairo_pdf(file.path(ROOT, "Fig3c_coconvergence.pdf"), width = 7.4, height = 6.8, family = "Arial")
par(fig = c(0.135, 0.965, 0.135, 0.845), mar = c(0,0,0,0))
plot.new(); plot.window(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05))
abline(h = 0, v = 0, col = "#E3E6EA", lwd = 0.9)
rect(0, 0, 1.05, 1.05, border = NA, col = "#FBF8FD")
points(CS$r1, CS$r2, pch = 21, cex = 0.60, bg = "#D3D7DC", col = "#FFFFFF", lwd = 0.3)
h <- CS[LOC == 2473]
points(h$r1, h$r2, pch = 21, cex = 2.0, bg = "#8E6FA8", col = INK, lwd = 1.3)
text(h$r1 - 0.05, h$r2 + 0.09, "chr22\nLOC 2473", cex = 0.56, col = INK2, adj = 1)
axis(1, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
axis(2, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
mtext("rho: Menorrhagia x Luminal A", side = 1, line = 2.3, cex = 0.68, col = INK2, las = 0)
mtext("rho: Ovarian cyst x Luminal A", side = 2, line = 2.3, cex = 0.68, col = INK2, las = 0)
mtext("concordant positive", side = 3, line = -1.2, adj = 0.88, cex = 0.50, col = "#A08DB8")
r <- cor(CS$r1, CS$r2, use = "complete.obs")
text(-1.0, 0.95, sprintf("r = %.2f   (n = %d blocks)", r, nrow(CS)), adj = 0.5, cex = 0.54, col = INK3)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "Two independent exposures converge on the same block", adj = 0.5, cex = 0.82, col = INK2)
text(0.5, 0.935, "Each point is one LD block; the chr22 block is the only strong outlier in the concordant-positive quadrant",
     adj = 0.5, cex = 0.53, col = INK3)
dev.off(); cat("Fig3c done\n")

# ==================== Fig 3d — significance burden ====================
ST <- B[, .(nblk = .N, nsig = sum(p < 0.05), nbonf = sum(p < 0.05/.N),
            minp = min(p), maxr = max(rho)), by = expo][match(bord, expo)]

cairo_pdf(file.path(ROOT, "Fig3d_block_burden.pdf"), width = 8.0, height = 6.6, family = "Arial")
par(fig = c(0.130, 0.735, 0.255, 0.810), mar = c(0,0,0,0))
bp <- barplot(ST$nsig, col = EXPCOL[ST$expo], border = NA, axes = FALSE,
              ylim = c(0, max(ST$nsig)*1.45), width = 0.62, space = 0.55, names.arg = rep("", 4))
for (i in 1:4) {
  if (!is.na(ST$nbonf[i]) && ST$nbonf[i] > 0)
    rect(bp[i]-0.31, 0, bp[i]+0.31, ST$nbonf[i], col = INK, border = NA)
  text(bp[i], ST$nsig[i] + max(ST$nsig)*0.06, ST$nsig[i], cex = 0.65, col = INK2)
  text(bp[i], -max(ST$nsig)*0.10, sprintf("of %d\ntested", ST$nblk[i]), cex = 0.50, col = INK3, xpd = NA)
}
axis(2, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
mtext("Blocks with significant local rg (P < 0.05)", side = 2, line = 2.4, cex = 0.68, col = INK2, las = 0)
text(bp, -max(ST$nsig)*0.22, EXPLAB[ST$expo], srt = 22, adj = c(1,1), cex = 0.62, col = INK2, xpd = NA)
abline(h = 0, col = INK3, lwd = 0.7)

# inset: minimum P values
par(fig = c(0.775, 0.975, 0.480, 0.780), mar = c(2.2,2.8,2.0,0.4), new = TRUE, mgp = c(1.3,0.4,0))
lp <- -log10(ST$minp)
plot.new(); plot.window(xlim = c(0.4,4.6), ylim = c(0, max(lp)*1.15))
for (i in 1:4) rect(i-0.32, 0, i+0.32, lp[i], col = EXPCOL[bord[i]], border = NA)
axis(2, cex.axis = 0.42, tck = -0.03, col.axis = INK3, lwd = 0.5)
axis(1, at = 1:4, labels = c("Men","Ova","End","Pro"), cex.axis = 0.42, tck = -0.03, col.axis = INK3, lwd = 0.5)
mtext("strongest signal", side = 3, line = 0.4, cex = 0.46, col = INK2)
mtext(expression(-log[10](P[min])), side = 2, line = 1.5, cex = 0.42, col = INK3, las = 0)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "Prolapse shares genetic basis broadly but never reaches experiment-wide significance",
     adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.935, "Dark segments indicate blocks passing an experiment-wide Bonferroni threshold within each exposure",
     adj = 0.5, cex = 0.52, col = INK3)
rect(0.038, 0.045, 0.052, 0.062, col = INK, border = NA)
text(0.058, 0.054, "passes Bonferroni", adj = 0, cex = 0.50, col = INK2)
dev.off(); cat("Fig3d done\n")
