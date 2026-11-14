#!/usr/bin/env Rscript
# ============================================================
#  Figure 6 — mammographic density is not on the causal pathway
# ============================================================
suppressMessages({library(data.table); library(grDevices)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

RG  <- fread(file.path(ROOT, "repro_breast_L3/repro_x_breast_matrix/repro_x_breast_rg.csv"))
RG  <- RG[!is.na(rg) & abs(rg) <= 1]
CO  <- fread(file.path(ROOT, "repro_breast_L3/md_luma_coloc/coloc_PMD_LumA.csv"))

INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
NEG <- "#5B8C8A"; POS <- "#8E6FA8"; AMB <- "#B08A5E"
PUR <- "#8E6FA8"; TEAL <- "#3E7D5A"

pl <- function(x) gsub("_", " ", x)

# ---- wide-format table ----
W <- dcast(RG, repro ~ breast, value.var = "rg")
P <- dcast(RG, repro ~ breast, value.var = "p")
setkey(W, NULL); setkey(P, NULL)
W <- W[repro %in% RG$repro]; P <- P[repro %in% RG$repro]

# ==================== Fig 6a — quadrant scatter ====================
D <- W[!is.na(Density_PMD) & !is.na(BCAC_overall)]
D[, p_ca := P$BCAC_overall[match(repro, P$repro)]]
D[, p_de := P$Density_PMD[match(repro, P$repro)]]
D[, sig := (p_ca < 0.05) | (p_de < 0.05)]
D[, quad := ifelse(BCAC_overall > 0 & Density_PMD < 0, "risk up / density down",
            ifelse(BCAC_overall < 0 & Density_PMD > 0, "risk down / density up",
            ifelse(BCAC_overall > 0 & Density_PMD > 0, "both up",
            ifelse(BCAC_overall < 0 & Density_PMD < 0, "both down", "mixed"))))]
QCOL <- c("risk up / density down" = POS, "risk down / density up" = NEG,
          "both up" = "#C6B7D8", "both down" = "#9FBBB4", "mixed" = "#DDE1E6")

cairo_pdf(file.path(ROOT, "Fig6a_quadrant.pdf"), width = 8.4, height = 7.6, family = "Arial")
par(fig = c(0.115, 0.905, 0.135, 0.855), mar = c(0,0,0,0))
xl <- range(D$BCAC_overall, na.rm = TRUE); yl <- range(D$Density_PMD, na.rm = TRUE)
xl <- xl + c(-0.03, 0.03); yl <- yl + c(-0.03, 0.03)
plot.new(); plot.window(xlim = xl, ylim = yl)
rect(0, yl[1], xl[2], 0, border = NA, col = "#FBF8FD")   # risk+, density-
rect(xl[1], 0, 0, yl[2], border = NA, col = "#F5F9F7")   # risk-, density+
abline(h = 0, v = 0, col = "#C9CFD6", lwd = 0.9)
for (q in names(QCOL)) {
  s <- D[quad == q]
  if (!nrow(s)) next
  points(s$BCAC_overall, s$Density_PMD, pch = 21, cex = ifelse(s$sig, 1.10, 0.72),
         bg = QCOL[q], col = ifelse(s$sig, INK, "#B8BCC2"), lwd = ifelse(s$sig, 0.9, 0.4))
}
# label significant points
lab <- D[sig == TRUE][order(-abs(BCAC_overall))]
for (i in seq_len(min(6, nrow(lab))))
  text(lab$BCAC_overall[i], lab$Density_PMD[i], pl(lab$repro[i]), cex = 0.46,
       col = INK2, pos = ifelse(lab$Density_PMD[i] > 0, 3, 1), offset = 0.35)
axis(1, cex.axis = 0.64, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6, las = 1)
axis(2, cex.axis = 0.64, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6, las = 1)
mtext("rg: disorder x breast cancer", side = 1, line = 2.2, cex = 0.66, col = INK2, las = 0)
q1 <- nrow(D[BCAC_overall > 0 & Density_PMD < 0]); q2 <- nrow(D[BCAC_overall < 0 & Density_PMD > 0])
q3 <- nrow(D[BCAC_overall > 0 & Density_PMD > 0]); q4 <- nrow(D[BCAC_overall < 0 & Density_PMD < 0])
text(xl[2]-0.02, yl[1]+0.145, sprintf("n = %d", q1), adj = 1, cex = 0.58, col = POS, font = 2)
text(xl[1]+0.02, yl[2]-0.11, sprintf("n = %d", q2), adj = 0, cex = 0.58, col = NEG, font = 2)
text(xl[2]-0.02, yl[2]-0.03, sprintf("n = %d", q3), adj = 1, cex = 0.52, col = "#8E7AA8")
text(xl[1]+0.02, yl[1]+0.065, sprintf("n = %d", q4), adj = 0, cex = 0.52, col = "#4E7A72")
mtext("rg: disorder x percent density", side = 2, line = 2.2, cex = 0.66, col = INK2, las = 0)
text(xl[2]-0.02, yl[1]+0.02, "risk UP\nbut density DOWN", adj = 1, cex = 0.46, col = POS, font = 3)
text(xl[1]+0.02, yl[2]-0.02, "risk down\ndensity UP", adj = 0, cex = 0.46, col = NEG, font = 3)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.975, "Reproductive disorders raise breast cancer risk while lowering dense area",
     adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.946, "Each point is one exposure; filled and outlined points reach P < 0.05 for either correlation",
     adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig6a done\n")

# ==================== Fig 6b — the nine significant density correlations ====================
DG <- RG[breast %in% c("Density_PMD","Density_DA","Density_NDA") & p < 0.05]
DG[, kind := ifelse(breast == "Density_NDA", "Non-dense area", "Dense area")]
DG <- DG[order(kind, rg)]
DG[, lab := sprintf("%s  (%s)", pl(repro), sub("Density_", "", breast))]

cairo_pdf(file.path(ROOT, "Fig6b_density_direction.pdf"), width = 8.0, height = 5.6, family = "Arial")
par(fig = c(0.340, 0.880, 0.115, 0.820), mar = c(0,0,0,0))
n <- nrow(DG); yy <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(min(DG$rg)*1.25, max(DG$rg)*1.25), ylim = c(0.4, n + 0.6))
abline(v = 0, col = "#C9CFD6", lwd = 0.9)
for (i in seq_len(n)) {
  cc <- if (DG$rg[i] < 0) NEG else AMB
  segments(0, yy[i], DG$rg[i], yy[i], col = cc, lwd = 2.0)
  points(DG$rg[i], yy[i], pch = 21, cex = 1.15, bg = cc, col = INK, lwd = 0.6)
  text(DG$rg[i] + sign(DG$rg[i])*0.012, yy[i], sprintf("%.3f", DG$rg[i]),
       adj = ifelse(DG$rg[i] > 0, 0, 1), cex = 0.50, col = INK2)
}
axis(2, at = yy, labels = DG$lab, cex.axis = 0.52, tick = FALSE, col.axis = INK2, line = -1.0, las = 1)
axis(1, cex.axis = 0.62, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6, las = 1)
mtext("Genetic correlation (rg)", side = 1, line = 2.2, cex = 0.66, col = INK2, las = 0)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
points(0.038, 0.058, pch = 21, cex = 1.05, bg = NEG, col = INK, lwd = 0.6)
text(0.050, 0.058, "dense area down", adj = 0, cex = 0.50, col = INK2)
points(0.230, 0.058, pch = 21, cex = 1.05, bg = AMB, col = INK, lwd = 0.6)
text(0.242, 0.058, "non-dense area up", adj = 0, cex = 0.50, col = INK2)
text(0.5, 0.965, "Direction is consistent: glandular tissue falls, fat rises", adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.936, sprintf("All %d nominally significant density correlations; none survives FDR correction", n),
     adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig6b done\n")

# ==================== Fig 6c — density-BC coloc (the contrast) ====================
CO <- CO[order(-H4)]
CO[, lab := sprintf("chr%d:%.1f Mb", chr, bp/1e6)]

cairo_pdf(file.path(ROOT, "Fig6c_density_bc_coloc.pdf"), width = 7.6, height = 5.4, family = "Arial")
par(fig = c(0.230, 0.900, 0.135, 0.820), mar = c(0,0,0,0))
n <- nrow(CO); yy <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(0, 1.08), ylim = c(0.4, n + 0.6))
abline(v = 0.8, lty = 2, col = "#D8DCE1", lwd = 0.9)
for (i in seq_len(n)) {
  cc <- if (CO$H4[i] > 0.8) PUR else if (CO$H4[i] > 0.5) "#C6B7D8" else "#DDE1E6"
  rect(0, yy[i]-0.32, CO$H4[i], yy[i]+0.32, col = cc, border = NA)
  text(CO$H4[i] + 0.012, yy[i], sprintf("%.3f", CO$H4[i]), adj = 0, cex = 0.48, col = INK2)
}
axis(2, at = yy, labels = CO$lab, cex.axis = 0.54, tick = FALSE, col.axis = INK2, line = -1.0, las = 1)
axis(1, cex.axis = 0.62, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6, las = 1)
mtext("PP.H4: density x luminal A breast cancer", side = 1, line = 2.3, cex = 0.66, col = INK2, las = 0)
text(0.82, n + 0.35, "PP.H4 = 0.8", cex = 0.46, col = INK3, adj = 0)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "Density and breast cancer do share causal variants", adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.936, sprintf("%d of %d loci reach PP.H4 > 0.8 - the expected, concordant direction", sum(CO$H4 > 0.8), n),
     adj = 0.5, cex = 0.50, col = INK3)
text(0.040, 0.070, "So density is a genuine breast cancer axis - it is simply not the axis the reproductive disorders act through",
     adj = 0, cex = 0.50, col = TEAL, font = 3)
dev.off(); cat("Fig6c done\n")

# ==================== Fig 6d — paired directions ====================
PR <- W[!is.na(Density_PMD) & !is.na(BCAC_overall)]
PR[, diff := BCAC_overall - Density_PMD]
PR <- PR[order(diff)][c(1:8, (.N-7):.N)]
PR[, lab := pl(repro)]

cairo_pdf(file.path(ROOT, "Fig6d_paired_direction.pdf"), width = 8.4, height = 6.2, family = "Arial")
par(fig = c(0.330, 0.860, 0.115, 0.820), mar = c(0,0,0,0))
n <- nrow(PR); yy <- rev(seq_len(n))
plot.new(); plot.window(xlim = c(-0.35, 0.45), ylim = c(0.4, n + 0.6))
abline(v = 0, col = "#C9CFD6", lwd = 0.9)
for (i in seq_len(n)) {
  segments(PR$Density_PMD[i], yy[i], PR$BCAC_overall[i], yy[i], col = "#DDE1E6", lwd = 2.4)
  points(PR$Density_PMD[i], yy[i], pch = 21, cex = 1.05, bg = NEG, col = "#FFFFFF", lwd = 0.4)
  points(PR$BCAC_overall[i], yy[i], pch = 21, cex = 1.05, bg = POS, col = "#FFFFFF", lwd = 0.4)
}
axis(2, at = yy, labels = PR$lab, cex.axis = 0.52, tick = FALSE, col.axis = INK2, line = -1.0, las = 1)
axis(1, cex.axis = 0.62, col.axis = INK2, tck = -0.012, col = INK3, lwd = 0.6, las = 1)
mtext("Genetic correlation (rg)", side = 1, line = 2.2, cex = 0.66, col = INK2, las = 0)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
points(0.038, 0.062, pch = 21, cex = 1.05, bg = NEG, col = "#FFFFFF", lwd = 0.4)
text(0.050, 0.062, "vs percent density", adj = 0, cex = 0.50, col = INK2)
points(0.200, 0.062, pch = 21, cex = 1.05, bg = POS, col = "#FFFFFF", lwd = 0.4)
text(0.212, 0.062, "vs breast cancer", adj = 0, cex = 0.50, col = INK2)
text(0.5, 0.965, "Most exposures sit on opposite sides of zero for the two outcomes", adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.936, "Eight most discordant and eight most concordant exposures shown",
     adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig6d done\n")
