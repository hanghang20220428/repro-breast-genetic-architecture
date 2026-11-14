#!/usr/bin/env Rscript
# ============================================================
#  Figure 2b / 2c / 2d  (rich versions, with insets)
# ============================================================
suppressMessages({library(data.table); library(grDevices)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

M <- fread(file.path(ROOT, "repro_breast_L3/repro_x_breast_matrix/repro_x_breast_rg.csv"))
TT<- fread(file.path(ROOT, "repro_breast_L3/repro_x_breast_tsmr/tsmr_full.csv"))
M <- M[!is.na(rg)]; M <- M[abs(rg) <= 1]
M[, sig := p < 0.05][, fsig := !is.na(fdr) & fdr < 0.05]

NEG <- "#5B8C8A"; POS <- "#8E6FA8"
INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
CATCOL <- c(Cancer = "#8E6FA8", Benign = "#5B8C8A", Density = "#B08A5E", Other = "#9AA0A6")

bord <- c("BCAC_overall","LumA","LumB","HER2","TNBC","Breast_insitu",
          "Breast_disease_NOS","Breast_hypertrophy","Benign_breast","Breast_lump",
          "Density_PMD","Density_DA","Density_NDA","Lactation_disorder")
blab <- c("BC overall","LumA","LumB","HER2","TNBC","In situ",
          "Disease NOS","Hypertrophy","Benign","Lump","PMD","DA","NDA","Lactation")
bcat <- c(rep("Cancer",6), rep("Benign",4), rep("Density",3), "Other"); names(bcat) <- bord

# ==================== Fig 2b ====================
cn <- M[, .(n = sum(sig), nf = sum(fsig), mn = min(rg), mx = max(rg), md = median(rg)),
         by = breast][match(bord, breast)]

cairo_pdf(file.path(ROOT, "Fig2b_sig_counts.pdf"), width = 9.0, height = 7.4, family = "Arial")

# --- top: significant counts ---
par(fig = c(0.075, 0.985, 0.505, 0.845), mar = c(0, 5.2, 0, 0.6))
bp <- barplot(cn$n, col = CATCOL[bcat[cn$breast]], border = NA, axes = FALSE,
              ylim = c(0, max(cn$n)*1.30), width = 0.72, space = 0.42, names.arg = rep("", 14))
for (i in 1:14) if (!is.na(cn$nf[i]) && cn$nf[i] > 0)
  rect(bp[i]-0.36, 0, bp[i]+0.36, cn$nf[i], col = INK, border = NA)
axis(2, at = seq(0, 12, 3), cex.axis = 0.66, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6)
mtext("No. of reproductive\ndisorders (P < 0.05)", side = 2, line = 2.7, cex = 0.62, col = INK2, las = 0)
for (i in 1:14) text(bp[i], cn$n[i] + max(cn$n)*0.05, cn$n[i], cex = 0.58, col = INK2)
for (i in 1:14) if (!is.na(cn$nf[i]) && cn$nf[i] > 0)
  text(bp[i], cn$nf[i]*0.5, cn$nf[i], cex = 0.56, col = "white", font = 2)
abline(h = 0, col = INK3, lwd = 0.7)

# --- bottom: range of rg signs ---
par(fig = c(0.075, 0.985, 0.185, 0.475), mar = c(4.6, 5.2, 0, 0.6), new = TRUE)
plot(NA, xlim = c(0.3, 14.7), ylim = c(-0.75, 0.75), axes = FALSE, xlab = "", ylab = "")
abline(h = 0, col = "#D8DCE1", lwd = 0.9, lty = 1)
for (i in 1:14) {
  segments(i, cn$mn[i], i, cn$mx[i], col = "#C9CFD6", lwd = 1.6)
  points(i, cn$mn[i], pch = 16, cex = 0.55, col = NEG)
  points(i, cn$mx[i], pch = 16, cex = 0.55, col = POS)
  points(i, cn$md[i], pch = 21, cex = 0.85, bg = "#FFFFFF", col = INK, lwd = 0.9)
}
axis(2, at = seq(-0.6, 0.6, 0.3), cex.axis = 0.62, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6)
mtext("rg range", side = 2, line = 2.7, cex = 0.62, col = INK2, las = 0)
text(1:14, -0.71, blab, srt = 45, adj = c(1,1), cex = 0.62, col = INK2, xpd = NA)
legend("bottomleft", bty = "n", cex = 0.54, pch = c(16,16,21), pt.cex = c(0.7,0.7,1.0),
       col = c(NEG, POS, INK), pt.bg = "#FFFFFF",
       legend = c("minimum rg","maximum rg","median rg"), text.col = INK2)

# --- inset: category composition of FDR hits ---
par(fig = c(0.838, 0.978, 0.628, 0.762), mar = c(0.2,0.2,1.5,0.2), new = TRUE)
fd <- M[fsig == TRUE, .N, by = breast]
fd[, cat := bcat[breast]]
cc <- fd[, .N, by = cat]
if (nrow(cc)) {
  pie(cc$N, col = CATCOL[cc$cat], border = "white", radius = 0.92, labels = NA, lwd = 1.1)
  text(0, 0.04, sum(cc$N), cex = 0.92, font = 2, col = INK)
  text(0, -0.34, "FDR hits", cex = 0.48, col = INK3)
}
mtext("FDR hits by class", side = 3, line = 0.4, cex = 0.46, col = INK2)

# --- legend and title ---
par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
yy <- c(0.072, 0.050, 0.028); cc2 <- c("Cancer","Benign","Density")
for (i in 1:3) {
  rect(0.075, yy[i]-0.006, 0.090, yy[i]+0.006, col = CATCOL[cc2[i]], border = NA)
  text(0.095, yy[i], c("Cancer subtypes","Benign / non-malignant","Density")[i], adj = 0, cex = 0.54, col = INK2)
}
rect(0.330, 0.060, 0.345, 0.072, col = INK, border = NA)
text(0.350, 0.066, "FDR < 0.05", adj = 0, cex = 0.54, col = INK2)
text(0.5, 0.985, "Shared signal concentrates on benign proliferative phenotypes",
     adj = 0.5, cex = 0.84, col = INK2)
text(0.5, 0.958, "Upper: number of reproductive disorders correlated with each outcome.  Lower: direction and spread of rg",
     adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig2b done\n")

# ==================== Fig 2c ====================
M[, logp := -log10(p)]
M[, cat := bcat[breast]]; M[is.na(cat), cat := "Other"]

cairo_pdf(file.path(ROOT, "Fig2c_volcano.pdf"), width = 8.4, height = 7.2, family = "Arial")

# main scatter
par(fig = c(0.095, 0.800, 0.145, 0.855), mar = c(0,0,0,0))
yl <- max(M$logp)*1.08
xl <- c(min(M$rg, na.rm=TRUE)-0.03, max(M$rg, na.rm=TRUE)+0.03)
plot.new(); plot.window(xlim = xl, ylim = c(0, yl))
abline(h = -log10(0.05), lty = 2, col = INK3, lwd = 0.8)
abline(v = 0, col = "#E3E6EA", lwd = 0.9)
for (cc in c("Other","Density","Benign","Cancer")) {
  s <- M[cat == cc]
  points(s$rg, s$logp, pch = 21, cex = 0.62, bg = CATCOL[cc], col = "#FFFFFF", lwd = 0.3)
}
F <- M[fsig == TRUE]
points(F$rg, F$logp, pch = 21, cex = 1.35, bg = NA, col = INK, lwd = 0.95)
axis(1, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
axis(2, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
mtext("Genetic correlation (rg)", side = 1, line = 2.2, cex = 0.70, col = INK2, las = 0)
mtext(expression(-log[10](italic(P))), side = 2, line = 2.0, cex = 0.70, col = INK2, las = 0)
text(xl[1]+0.02, -log10(0.05) + yl*0.02, "P = 0.05", adj = 0, cex = 0.52, col = INK3)

# top marginal density
par(fig = c(0.095, 0.800, 0.862, 0.930), mar = c(0,0,0,0), new = TRUE)
dx <- density(M$rg, n = 512, from = xl[1], to = xl[2])
plot.new(); plot.window(xlim = xl, ylim = c(0, max(dx$y)*1.12))
polygon(c(xl[1], dx$x[dx$x <= 0], 0), c(0, dx$y[dx$x <= 0], 0), col = "#9FBBB4", border = NA)
polygon(c(0, dx$x[dx$x >= 0], xl[2]), c(0, dx$y[dx$x >= 0], 0), col = "#C9BFD8", border = NA)
abline(v = 0, col = INK3, lwd = 0.7)
mtext("Density", side = 2, line = 0.3, cex = 0.42, col = INK3, las = 0)

# right marginal density
par(fig = c(0.808, 0.885, 0.145, 0.855), mar = c(0,0,0,0), new = TRUE)
hy <- hist(M$logp, breaks = 30, plot = FALSE)
plot.new(); plot.window(xlim = c(0, max(hy$counts)*1.1), ylim = c(0, yl))
for (b in seq_along(hy$counts)) {
  if (hy$counts[b] <= 0) next
  rect(0, hy$breaks[b], hy$counts[b], hy$breaks[b+1], col = "#C9CFD6", border = "white", lwd = 0.3)
}
mtext("Count", side = 1, line = 0.3, cex = 0.44, col = INK3, las = 0)

# inset: quadrant counts
par(fig = c(0.890, 0.995, 0.640, 0.845), mar = c(1.4,1.6,1.9,0.4), new = TRUE)
q <- matrix(0, 2, 2)
q[1,1] <- M[rg<0 & p<0.05, .N]; q[1,2] <- M[rg>0 & p<0.05, .N]
q[2,1] <- M[rg<0 & p>=0.05, .N]; q[2,2] <- M[rg>0 & p>=0.05, .N]
plot.new(); plot.window(xlim = c(0,2), ylim = c(0,2))
for (i in 1:2) for (j in 1:2) {
  col <- if (i==1 && j==1) "#9FBBB4" else if (i==1 && j==2) "#C9BFD8" else "#F0F1F3"
  rect(j-1, 2-i, j, 3-i, col = col, border = "white", lwd = 1.2)
  text(j-0.5, 2.5-i, q[i,j], cex = 0.62, font = 2, col = INK)
}
mtext("neg rg", side = 1, line = 0.5, at = 0.5, cex = 0.42, col = INK3)
mtext("pos rg", side = 1, line = 0.5, at = 1.5, cex = 0.42, col = INK3)
mtext("sig", side = 2, line = 0.4, at = 1.5, cex = 0.42, col = INK3, las = 0)
mtext("ns", side = 2, line = 0.4, at = 0.5, cex = 0.42, col = INK3, las = 0)
mtext("2x2 breakdown", side = 3, line = 0.5, cex = 0.46, col = INK2)

# inset: median by outcome
par(fig = c(0.890, 0.995, 0.145, 0.615), mar = c(1.6,1.8,1.9,0.3), new = TRUE)
agg2 <- M[, .(n = sum(sig)), by = breast][match(bord, breast)]
plot.new(); plot.window(xlim = c(0, max(agg2$n)*1.15), ylim = c(0.5, 14.5))
for (i in 1:14) {
  y <- 15 - i
  rect(0, y-0.36, agg2$n[i], y+0.36, col = CATCOL[bcat[bord[i]]], border = NA)
}
axis(2, at = 15-1:14, labels = blab, cex.axis = 0.40, tick = FALSE, col.axis = INK2, line = -0.9)
mtext("Pairs per outcome", side = 3, line = 0.5, cex = 0.46, col = INK2)

# labels and legend
par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
lab <- F[order(-logp)][1:min(4, nrow(F))]
lab[, txt := paste0(gsub("_"," ",repro), " - ", blab[match(breast, bord)])]
lab[, ty := yl*(0.93 - 0.072*(seq_len(.N)-1))]
for (i in seq_len(nrow(lab))) {
  segments(lab$rg[i], lab$logp[i], 0.470, lab$ty[i], col = INK3, lwd = 0.32)
  text(0.478, lab$ty[i], lab$txt[i], adj = 0, cex = 0.48, col = INK2)
}
yy <- c(0.740, 0.705, 0.670); cc2 <- c("Cancer","Benign","Density")
for (i in 1:3) {
  points(0.130, yy[i], pch = 21, cex = 0.78, bg = CATCOL[cc2[i]], col = "#FFFFFF", lwd = 0.3)
  text(0.143, yy[i], c("Cancer subtypes","Benign / non-malignant","Density")[i], adj = 0, cex = 0.50, col = INK2)
}
points(0.130, 0.635, pch = 21, cex = 1.25, bg = NA, col = INK, lwd = 0.95)
text(0.143, 0.635, "FDR < 0.05", adj = 0, cex = 0.50, col = INK2)
text(0.5, 0.985, "Effect-size distribution across all testable pairs", adj = 0.5, cex = 0.84, col = INK2)
text(0.5, 0.962, "Marginal densities shown along both axes; inset panels summarise the 2x2 split and per-outcome burden",
     adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig2c done\n")

# ==================== Fig 2d ====================
ivw <- TT[method == "Inverse variance weighted"]
rg_n <- M[, .(n_rg = sum(sig), n_any = .N), by = repro]
mr_n <- ivw[, .(n_mr = sum(p < 0.05, na.rm = TRUE), nIV = max(nsnp, na.rm = TRUE)), by = repro]
D <- merge(rg_n, mr_n, by = "repro", all.x = TRUE)
D[is.na(n_mr), n_mr := 0]
D <- D[order(-n_rg, -n_mr)][1:20]
D[, lab := gsub("_"," ",repro)]

cairo_pdf(file.path(ROOT, "Fig2d_rg_vs_MR.pdf"), width = 8.6, height = 7.4, family = "Arial")

par(fig = c(0.300, 0.800, 0.185, 0.855), mar = c(0,0,0,0))
n <- nrow(D); yy <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(0, max(D$n_rg)+1.6), ylim = c(0.4, n+0.9))
for (i in seq_len(n)) {
  segments(D$n_mr[i], yy[i], D$n_rg[i], yy[i], col = "#D8DCE1", lwd = 2.2)
  points(D$n_rg[i], yy[i], pch = 21, cex = 1.10, bg = POS, col = "#FFFFFF", lwd = 0.4)
  points(D$n_mr[i], yy[i], pch = 21, cex = 1.10, bg = NEG, col = "#FFFFFF", lwd = 0.4)
}
axis(2, at = yy, labels = D$lab, cex.axis = 0.56, tick = FALSE, col.axis = INK2, line = -0.9, las = 1)
axis(1, cex.axis = 0.66, col.axis = INK2, tck = -0.014, col = INK3, lwd = 0.6)
mtext("No. of breast outcomes significant", side = 1, line = 2.2, cex = 0.68, col = INK2, las = 0)
abline(v = 0, col = INK3, lwd = 0.7)

# inset 1: 2x2 concordance
par(fig = c(0.815, 0.995, 0.590, 0.840), mar = c(1.8,2.2,2.0,0.4), new = TRUE)
DA <- merge(rg_n, mr_n, by = "repro", all.x = TRUE); DA[is.na(n_mr), n_mr := 0]
q <- c(nrow(DA[n_rg>0 & n_mr>0]), nrow(DA[n_rg>0 & n_mr==0]),
       nrow(DA[n_rg==0 & n_mr>0]), nrow(DA[n_rg==0 & n_mr==0]))
plot.new(); plot.window(xlim = c(0,2), ylim = c(0,2))
vv <- matrix(q, 2, 2, byrow = TRUE)
ccs <- c("#8E6FA8","#C9BFD8","#E8E8EA","#F2F3F4")
for (i in 1:2) for (j in 1:2) {
  rect(j-1, 2-i, j, 3-i, col = ccs[(i-1)*2+j], border = "white", lwd = 1.2)
  text(j-0.5, 2.5-i, vv[i,j], cex = 0.72, font = 2, col = INK)
}
mtext("MR +", side = 1, line = 0.4, at = 0.5, cex = 0.44, col = INK3)
mtext("MR -", side = 1, line = 0.4, at = 1.5, cex = 0.44, col = INK3)
mtext("rg +", side = 2, line = 0.4, at = 1.5, cex = 0.44, col = INK3, las = 0)
mtext("rg -", side = 2, line = 0.4, at = 0.5, cex = 0.44, col = INK3, las = 0)
mtext("Exposure-level agreement", side = 3, line = 0.5, cex = 0.46, col = INK2)

# inset 2: instrument count distribution
par(fig = c(0.815, 0.995, 0.185, 0.555), mar = c(2.0,2.4,2.0,0.4), new = TRUE)
iv <- ivw[!is.na(nsnp)]
hv <- hist(iv$nsnp, breaks = c(0,2,4,8,16,32,64,128), plot = FALSE)
plot.new(); plot.window(xlim = c(0.6, 7.4), ylim = c(0, max(hv$counts)*1.15))
for (b in seq_along(hv$counts))
  rect(b-0.38, 0, b+0.38, hv$counts[b], col = "#C9BFD8", border = "white", lwd = 0.4)
axis(2, cex.axis = 0.44, tck = -0.03, col.axis = INK3, lwd = 0.5)
axis(1, at = 1:7, labels = c("0-2","3-4","5-8","9-16","17-32","33-64",">64"),
     cex.axis = 0.40, tck = -0.03, col.axis = INK3, lwd = 0.5)
mtext("Instrument count (IVW)", side = 3, line = 0.5, cex = 0.46, col = INK2)
mtext("No. of pairs", side = 2, line = 1.4, cex = 0.42, col = INK3, las = 0)

# legend and title
par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
points(0.240, 0.098, pch = 21, cex = 1.05, bg = POS, col = "#FFFFFF", lwd = 0.4)
text(0.253, 0.098, "LDSC rg  (P < 0.05)", adj = 0, cex = 0.54, col = INK2)
points(0.240, 0.070, pch = 21, cex = 1.05, bg = NEG, col = "#FFFFFF", lwd = 0.4)
text(0.253, 0.070, "TSMR IVW  (P < 0.05)", adj = 0, cex = 0.54, col = INK2)
n_no <- 40 - uniqueN(ivw$repro)
text(0.760, 0.085, sprintf("%d of 40 exposures\nlack an IVW estimate", n_no), cex = 0.50, col = INK3, adj = 0)
text(0.5, 0.985, "Genetic correlation is broad; causal evidence is narrow", adj = 0.5, cex = 0.84, col = INK2)
text(0.5, 0.962, "Each row compares the number of breast outcomes reaching P < 0.05 under LDSC versus TSMR for the same exposure",
     adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig2d done\n")
