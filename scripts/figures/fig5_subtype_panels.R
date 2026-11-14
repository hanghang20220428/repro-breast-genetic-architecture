#!/usr/bin/env Rscript
# ============================================================
#  Figure 5 — subtype specificity (ER+ / Luminal A axis)
# ============================================================
suppressMessages({library(data.table); library(grDevices)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

MB <- fread(file.path(ROOT, "repro_breast_L3/xbp1_metabric.csv"))
MB <- MB[!is.na(XBP1)]
AR <- fread(file.path(ROOT, "repro_breast_L3/androgen/androgen_score_by_celltype.csv"))
RG <- fread(file.path(ROOT, "repro_breast_L3/repro_x_breast_matrix/repro_x_breast_rg.csv"))
RG <- RG[!is.na(rg) & abs(rg) <= 1]

INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
PUR <- "#8E6FA8"; TEAL <- "#3E7D5A"; ROSE <- "#B5657E"
SUBCOL <- c(LumA="#8E6FA8", LumB="#B5657E", Her2="#4A7C9B", Basal="#B08A5E",
            Normal="#5B8C8A", `claudin-low`="#7A6E9B", NC="#C9CFD6")

# ==================== Fig 5a — XBP1 by PAM50 ====================
ord <- c("LumA","LumB","Normal","Her2","claudin-low","Basal","NC")
MB[, sub := factor(CLAUDIN_SUBTYPE, levels = ord)]
gg <- MB[!is.na(sub), .(med = median(XBP1), n = .N), by = sub][order(-med)]
lev <- as.character(gg$sub)
MB[, sub := factor(as.character(sub), levels = lev)]

cairo_pdf(file.path(ROOT, "Fig5a_pam50.pdf"), width = 7.8, height = 5.6, family = "Arial")
par(fig = c(0.115, 0.975, 0.130, 0.830), mar = c(0,0,0,0))
yl <- range(MB$XBP1, na.rm = TRUE)
plot.new(); plot.window(xlim = c(0.4, length(lev)+0.6), ylim = c(yl[1]-0.3, yl[2]+0.9))
for (i in seq_along(lev)) {
  s <- MB[sub == lev[i]]
  boxplot.stats(s$XBP1)$stats
  q <- quantile(s$XBP1, c(.25,.5,.75), na.rm = TRUE)
  lo <- quantile(s$XBP1, .05, na.rm = TRUE); hi <- quantile(s$XBP1, .95, na.rm = TRUE)
  segments(i, lo, i, hi, col = "#C9CFD6", lwd = 1.1)
  rect(i-0.32, q[1], i+0.32, q[3], col = adjustcolor(SUBCOL[lev[i]], 0.45),
       border = SUBCOL[lev[i]], lwd = 1.0)
  segments(i-0.32, q[2], i+0.32, q[2], col = INK, lwd = 1.6)
  # jitter
  set.seed(1); jx <- i + runif(nrow(s), -0.20, 0.20)
  points(jx, s$XBP1, pch = 16, cex = 0.22, col = adjustcolor(SUBCOL[lev[i]], 0.30))
  text(i, yl[2]+0.55, sprintf("n=%d", nrow(s)), cex = 0.48, col = INK3)
  text(i, yl[2]+0.20, sprintf("%.2f", median(s$XBP1)), cex = 0.52, col = INK2, font = 2)
}
axis(1, at = 1:length(lev), labels = lev, cex.axis = 0.62, tick = FALSE, col.axis = INK2, line = -0.8, las = 1)
axis(2, cex.axis = 0.64, col.axis = INK2, tck = -0.013, col = INK3, lwd = 0.6, las = 1)
mtext("XBP1 expression (log2)", side = 2, line = 2.3, cex = 0.68, col = INK2, las = 0)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "XBP1 is highest in luminal tumours and lowest in basal", adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.938, sprintf("METABRIC, n = %d tumours; boxes show median and IQR, whiskers 5th-95th percentile", nrow(MB)),
     adj = 0.5, cex = 0.50, col = INK3)
text(0.042, 0.028, "LumA 13.30  vs  Basal 10.93  (2.37 log2 units)", adj = 0, cex = 0.54, col = INK2, font = 3)
dev.off(); cat("Fig5a done\n")

# ==================== Fig 5b — ER status x HER2 ====================
MB[, grp := paste0(ifelse(ER_STATUS == "Positive", "ER+", "ER-"),
                   " / ", ifelse(HER2_STATUS == "Positive", "HER2+", "HER2-"))]
G2 <- MB[!is.na(ER_STATUS) & !is.na(HER2_STATUS),
         .(med = median(XBP1), n = .N), by = grp][order(-med)]

cairo_pdf(file.path(ROOT, "Fig5b_er_her2.pdf"), width = 7.4, height = 5.6, family = "Arial")
par(fig = c(0.115, 0.965, 0.140, 0.830), mar = c(0,0,0,0))
plot.new(); plot.window(xlim = c(0.4, nrow(G2)+0.6), ylim = c(yl[1]-0.3, yl[2]+0.9))
for (i in seq_len(nrow(G2))) {
  s <- MB[grp == G2$grp[i]]
  q <- quantile(s$XBP1, c(.25,.5,.75), na.rm = TRUE)
  lo <- quantile(s$XBP1, .05, na.rm = TRUE); hi <- quantile(s$XBP1, .95, na.rm = TRUE)
  cc <- if (grepl("^ER\\+", G2$grp[i])) PUR else "#B08A5E"
  segments(i, lo, i, hi, col = "#C9CFD6", lwd = 1.1)
  rect(i-0.28, q[1], i+0.28, q[3], col = adjustcolor(cc, 0.45), border = cc, lwd = 1.0)
  segments(i-0.28, q[2], i+0.28, q[2], col = INK, lwd = 1.6)
  set.seed(2); points(i + runif(nrow(s), -0.21, 0.21), s$XBP1, pch = 16, cex = 0.22,
                      col = adjustcolor(cc, 0.28))
  text(i, yl[2]+0.55, sprintf("n=%d", nrow(s)), cex = 0.46, col = INK3)
  text(i, yl[2]+0.20, sprintf("%.2f", q[2]), cex = 0.52, col = INK2, font = 2)
}
axis(1, at = 1:nrow(G2), labels = G2$grp, cex.axis = 0.60, tick = FALSE, col.axis = INK2, line = -0.8, las = 1)
axis(2, cex.axis = 0.64, col.axis = INK2, tck = -0.013, col = INK3, lwd = 0.6, las = 1)
mtext("XBP1 expression (log2)", side = 2, line = 2.3, cex = 0.68, col = INK2, las = 0)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "The gradient tracks ER status, not HER2 status", adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.938, "ER+ median 13.24 versus ER- 11.70 (Mann-Whitney P = 2.3e-138)", adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig5b done\n")

# ==================== Fig 5c — AR across breast cell types ====================
AR <- AR[order(AR_expr)]
AR[, lab := gsub("_", " ", celltype)]
AR[, hi := celltype %in% c("luminal_cell_hormone_responsive","basal_myoepithelial_cell")]

cairo_pdf(file.path(ROOT, "Fig5c_AR_celltype.pdf"), width = 8.2, height = 6.2, family = "Arial")
par(fig = c(0.360, 0.820, 0.115, 0.830), mar = c(0,0,0,0))
n <- nrow(AR); yy <- seq_len(n)
plot.new(); plot.window(xlim = c(0, max(AR$AR_expr)*1.18), ylim = c(0.4, n + 0.6))
for (i in seq_len(n)) {
  cc <- if (AR$celltype[i] == "luminal_cell_hormone_responsive") PUR
        else if (AR$celltype[i] == "basal_myoepithelial_cell") "#B08A5E" else "#C9CFD6"
  rect(0, yy[i]-0.34, AR$AR_expr[i], yy[i]+0.34, col = cc, border = NA)
  text(AR$AR_expr[i] + 0.03, yy[i], sprintf("%.3f", AR$AR_expr[i]), adj = 0, cex = 0.46,
       col = if (AR$hi[i]) INK2 else INK3)
}
axis(2, at = yy, labels = AR$lab, cex.axis = 0.52, tick = FALSE, col.axis = INK2, line = -1.0, las = 1)
axis(1, cex.axis = 0.62, col.axis = INK2, tck = -0.013, col = INK3, lwd = 0.6, las = 1)
mtext("Mean AR expression (normalised)", side = 1, line = 2.3, cex = 0.66, col = INK2, las = 0)

# cell count, right margin
par(fig = c(0.825, 0.995, 0.115, 0.830), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0.4, n + 0.6))
for (i in seq_len(n))
  text(0, yy[i], format(AR$n_cell[i], big.mark = ","), adj = 0, cex = 0.44,
       col = if (AR$hi[i]) INK2 else INK3)
mtext("cells", side = 3, line = 0.4, cex = 0.46, col = INK3, at = 0.15)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "Androgen receptor is confined to the luminal compartment", adj = 0.5, cex = 0.80, col = INK2)
text(0.5, 0.938, "Hormone Cell Atlas breast atlas; luminal hormone-responsive cells 2.627 versus basal myoepithelial 0.086 (31-fold)",
     adj = 0.5, cex = 0.50, col = INK3)
dev.off(); cat("Fig5c done\n")

# ==================== Fig 5d — subtype specificity of rg ====================
SUB <- c("BCAC_overall","LumA","LumB","HER2","TNBC")
SUBL <- c("Overall","LumA","LumB","HER2","TNBC")
S <- RG[breast %in% SUB]
cnt <- S[, .(n = sum(p < 0.05, na.rm=TRUE)), by = breast][match(SUB, breast)]
top <- S[, .(nsig = sum(p < 0.05)), by = repro][order(-nsig)][1:14]
S2 <- S[repro %in% top$repro]
Rm <- matrix(NA_real_, nrow(top), 5)
Pm <- Rm
for (k in seq_len(nrow(S2))) {
  i <- match(S2$repro[k], top$repro); j <- match(S2$breast[k], SUB)
  Rm[i, j] <- S2$rg[k]; Pm[i, j] <- S2$p[k]
}
lim <- max(abs(Rm), na.rm = TRUE)
pal <- colorRampPalette(c("#5B8C8A", "#F8F9FA", "#8E6FA8"))(201)
colf <- function(v) if (is.na(v)) "#FFFFFF" else pal[round((v+lim)/(2*lim)*200)+1]

cairo_pdf(file.path(ROOT, "Fig5d_subtype_rg.pdf"), width = 8.8, height = 6.6, family = "Arial")
par(fig = c(0.310, 0.800, 0.115, 0.830), mar = c(0,0,0,0))
NR <- nrow(top)
plot.new(); plot.window(xlim = c(0.5, 5.5), ylim = c(0.5, NR + 0.5))
for (i in 1:NR) for (j in 1:5)
  rect(j-0.5, NR+1-i-0.5, j+0.5, NR+1-i+0.5, col = colf(Rm[i,j]), border = "#FFFFFF", lwd = 0.6)
for (i in 1:NR) for (j in 1:5) {
  if (!is.na(Pm[i,j]) && Pm[i,j] < 0.05) points(j, NR+1-i, pch = 16, cex = 0.5, col = INK)
}
text(0.25, NR+1-(1:NR), gsub("_"," ",top$repro), adj = 1, cex = 0.50, col = INK2, xpd = NA)
text(1:5, 0.08, SUBL, adj = c(0.5, 1), cex = 0.56, col = INK2, xpd = NA)

# right: number significant
par(fig = c(0.815, 0.960, 0.415, 0.830), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0, 6.5), ylim = c(0.5, 5.5))
for (j in 1:5) {
  y <- 6 - j
  rect(0, y-0.33, cnt$n[j], y+0.33, col = if (j==2) PUR else "#C9BFD8", border = NA)
  text(cnt$n[j] + 0.25, y, cnt$n[j], adj = 0, cex = 0.52, col = INK2)
}
axis(2, at = 5:1, labels = FALSE, tick = FALSE)
mtext("Significant
exposures", side = 3, line = 0.5, cex = 0.46, col = INK2)

par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
text(0.5, 0.965, "Shared signal is concentrated in luminal A, not in basal-like or HER2 disease",
     adj = 0.5, cex = 0.78, col = INK2)
text(0.5, 0.938, sprintf("Top %d exposures by number of significant subtypes; dots mark P < 0.05", NR),
     adj = 0.5, cex = 0.50, col = INK3)
cbx <- 0.040; cby <- 0.055; cbw <- 0.130; cbh <- 0.011
for (k in 1:200) rect(cbx+(k-1)*cbw/200, cby, cbx+k*cbw/200, cby+cbh, col = pal[k], border = NA)
text(cbx, cby-0.014, sprintf("%.2f", -lim), cex = 0.46, col = INK3)
text(cbx+cbw/2, cby-0.014, "0", cex = 0.46, col = INK3)
text(cbx+cbw, cby-0.014, sprintf("%.2f", lim), cex = 0.46, col = INK3)
text(cbx+cbw/2, cby+cbh+0.012, "rg", cex = 0.46, col = INK2)
points(0.205, cby+0.005, pch = 16, cex = 0.48, col = INK)
text(0.213, cby+0.005, "P < 0.05", adj = 0, cex = 0.46, col = INK3)
dev.off(); cat("Fig5d done\n")
