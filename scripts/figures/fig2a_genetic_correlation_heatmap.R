#!/usr/bin/env Rscript
# ============================================================
#  Figure 2a  (rich)  main heatmap, marginal bars and two insets
# ============================================================
suppressMessages({library(data.table); library(grDevices)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

M <- fread(file.path(ROOT, "repro_breast_L3/repro_x_breast_matrix/repro_x_breast_rg.csv"))
M <- M[!is.na(rg)]
M <- M[abs(rg) <= 1]        # rg must lie in [-1,1]; 21 unusable pairs removed
M[, sig := p < 0.05][, fsig := !is.na(fdr) & fdr < 0.05]

NEG <- "#5B8C8A"; POS <- "#8E6FA8"
INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
CATCOL <- c(Cancer = "#8E6FA8", Benign = "#5B8C8A", Density = "#B08A5E", Other = "#9AA0A6")

bord <- c("BCAC_overall","LumA","LumB","HER2","TNBC","Breast_insitu",
          "Breast_disease_NOS","Breast_hypertrophy","Benign_breast","Breast_lump",
          "Density_PMD","Density_DA","Density_NDA","Lactation_disorder")
blab <- c("BC overall","LumA","LumB","HER2","TNBC","In situ",
          "Disease NOS","Hypertrophy","Benign","Lump","PMD","DA","NDA","Lactation")
bcat <- c(rep("Cancer",6), rep("Benign",4), rep("Density",3), "Other")
names(bcat) <- bord

agg <- M[, .(n_f = sum(fsig), n_s = sum(sig), mr = mean(rg)), by = repro][order(-n_f, -n_s, -mr)]
ro <- agg$repro; NR <- length(ro)

Rm <- matrix(NA_real_, NR, 14); Pm <- Rm; Fm <- Rm
for (k in seq_len(nrow(M))) {
  i <- match(M$repro[k], ro); j <- match(M$breast[k], bord)
  Rm[i,j] <- M$rg[k]; Pm[i,j] <- M$p[k]; Fm[i,j] <- M$fdr[k]
}
lim <- max(abs(Rm), na.rm = TRUE)
pal <- colorRampPalette(c(NEG, "#F8F9FA", POS))(201)
colf <- function(v) if (is.na(v)) "#FFFFFF" else pal[round((v+lim)/(2*lim)*200)+1]

cairo_pdf(file.path(ROOT, "Fig2a_rg_heatmap.pdf"),
          width = 11.4, height = 10.6, family = "Arial")

# ---- main heatmap ----
par(fig = c(0.255, 0.795, 0.135, 0.852), mar = c(0,0,0,0))
plot.new(); plot.window(xlim = c(0.5, 14.5), ylim = c(0.5, NR + 0.5))
for (i in 1:NR) for (j in 1:14)
  rect(j-0.5, NR+0.5-i, j+0.5, NR+0.5-i+1, col = colf(Rm[i,j]), border = "#FFFFFF", lwd = 0.45)
for (i in 1:NR) for (j in 1:14) {
  y <- NR + 1 - i
  if (!is.na(Fm[i,j]) && Fm[i,j] < 0.05)      points(j, y, pch = 16, cex = 0.44, col = INK)
  else if (!is.na(Pm[i,j]) && Pm[i,j] < 0.05) points(j, y, pch = 1, cex = 0.40, col = "#33383D", lwd = 0.6)
}
for (b in c(6.5, 10.5, 13.5)) segments(b, 0.5, b, NR+0.5, col = "#FFFFFF", lwd = 1.8)
text(0.25, NR+1-(1:NR), gsub("_"," ",ro), adj = 1, cex = 0.50, col = INK2, xpd = NA)
text(1:14, 0.15, blab, srt = 50, adj = c(1,1), cex = 0.52, col = INK2, xpd = NA)

# ---- top panel ----
par(fig = c(0.278, 0.795, 0.862, 0.918), mar = c(0,0,0,0), new = TRUE)
cn <- M[, .(n = sum(sig), nf = sum(fsig)), by = breast][match(bord, breast)]
plot.new(); plot.window(xlim = c(0.5,14.5), ylim = c(0, max(cn$n)*1.40))
for (j in 1:14) {
  rect(j-0.34, 0, j+0.34, cn$n[j], col = CATCOL[bcat[j]], border = NA)
  if (!is.na(cn$nf[j]) && cn$nf[j] > 0)
    rect(j-0.34, cn$n[j], j+0.34, cn$n[j]+0.32, col = INK, border = NA)
  text(j, cn$n[j] + max(cn$n)*0.11, cn$n[j], cex = 0.52, col = INK2)
}
abline(h = 0, col = INK3, lwd = 0.6)
mtext("Exposures", side = 2, line = 1.5, cex = 0.46, col = INK2, las = 0)

# ---- right panel ----
par(fig = c(0.800, 0.888, 0.135, 0.852), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0, 14.6), ylim = c(0.5, NR + 0.5))
for (i in 1:NR) {
  y <- NR + 1 - i
  rect(0, y-0.33, agg$n_s[i], y+0.33, col = "#C6B7D8", border = NA)
  if (agg$n_f[i] > 0) rect(agg$n_s[i], y-0.33, agg$n_s[i]+0.42, y+0.33, col = INK, border = NA)
}
axis(1, at = c(0,5,10), cex.axis = 0.46, tck = -0.02, col.axis = INK3, line = -0.85, lwd = 0.5)
mtext("Outcomes", side = 1, line = 1.0, cex = 0.48, col = INK2, las = 0)

# ---- inset 1: rg distribution ----
par(fig = c(0.800, 0.988, 0.600, 0.840), mar = c(2.2,2.6,2.0,0.4), new = TRUE, mgp = c(1.4,0.4,0))
rlo <- min(M$rg, na.rm=TRUE); rhi <- max(M$rg, na.rm=TRUE)
br <- seq(rlo - 0.01, rhi + 0.01, length.out = 32)
h <- hist(M$rg, breaks = br, plot = FALSE)
plot(NA, xlim = c(rlo, rhi), ylim = c(0, max(h$counts)*1.18), axes = FALSE, xlab="", ylab="")
for (b in seq_along(h$counts)) {
  cc <- if (h$mids[b] < 0) "#9FBBB4" else "#C9BFD8"
  rect(h$breaks[b], 0, h$breaks[b+1], h$counts[b], col = cc, border = "white", lwd = 0.35)
}
abline(v = 0, col = INK3, lwd = 0.7)
axis(1, cex.axis = 0.44, tck = -0.03, col.axis = INK3, lwd = 0.5, line = -0.7)
axis(2, cex.axis = 0.44, tck = -0.03, col.axis = INK3, lwd = 0.5, line = -0.6)
mtext("Distribution of rg", side = 3, line = 0.4, cex = 0.48, col = INK2)
text(0.57, max(h$counts)*0.95, sprintf("median %.3f", median(M$rg)), adj = 1, cex = 0.42, col = INK3)
text(0.57, max(h$counts)*0.78, sprintf("range %.2f to %.2f", min(M$rg), max(M$rg)), adj = 1, cex = 0.42, col = INK3)

# ---- inset 2: h2 versus |rg| ----
par(fig = c(0.800, 0.988, 0.135, 0.565), mar = c(2.4,2.8,2.0,0.5), new = TRUE, mgp = c(1.4,0.4,0))
plot(NA, xlim = c(0, max(M$h2_repro, na.rm=TRUE)*1.06), ylim = c(0, 0.85), axes = FALSE, xlab="", ylab="")
abline(h = 0.2, lty = 2, col = "#D8DCE1", lwd = 0.7)
points(M$h2_repro, abs(M$rg), pch = 21, cex = 0.52,
       bg = ifelse(M$sig, CATCOL[bcat[M$breast]], "#E8E8EA"), col = "#FFFFFF", lwd = 0.3)
axis(1, cex.axis = 0.44, tck = -0.03, col.axis = INK3, lwd = 0.5, line = -0.7)
axis(2, cex.axis = 0.44, tck = -0.03, col.axis = INK3, lwd = 0.5, line = -0.6)
mtext("h2 of exposure", side = 1, line = 1.2, cex = 0.44, col = INK2, las = 0)
mtext("|rg|", side = 2, line = 1.3, cex = 0.44, col = INK2, las = 0)
mtext("Power diagnostic", side = 3, line = 0.4, cex = 0.48, col = INK2)
r2 <- cor(M$h2_repro, abs(M$rg), use = "complete.obs")
text(max(M$h2_repro)*0.96, 0.80, sprintf("r = %.2f", r2), adj = 1, cex = 0.42, col = INK3)

# ---- colour scale, legend and title ----
par(fig = c(0,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
cbx <- 0.255; cby <- 0.052; cbw <- 0.125; cbh <- 0.010
for (k in 1:200) rect(cbx+(k-1)*cbw/200, cby, cbx+k*cbw/200, cby+cbh, col = pal[k], border = NA)
text(cbx, cby-0.013, sprintf("%.2f", -lim), cex = 0.46, col = INK3)
text(cbx+cbw/2, cby-0.013, "0", cex = 0.46, col = INK3)
text(cbx+cbw, cby-0.013, sprintf("%.2f", lim), cex = 0.46, col = INK3)
text(cbx+cbw/2, cby+cbh+0.012, "rg", cex = 0.46, col = INK2)

points(0.408, cby+0.005, pch = 16, cex = 0.48, col = INK)
text(0.416, cby+0.005, "FDR < 0.05", adj = 0, cex = 0.46, col = INK3)
points(0.470, cby+0.005, pch = 1, cex = 0.44, col = "#33383D", lwd = 0.6)
text(0.478, cby+0.005, "P < 0.05", adj = 0, cex = 0.46, col = INK3)

yy <- c(0.068, 0.052, 0.036)
cc <- c("Cancer","Benign","Density")
for (i in 1:3) {
  rect(0.545, yy[i]-0.004, 0.556, yy[i]+0.004, col = CATCOL[cc[i]], border = NA)
  text(0.561, yy[i], c("Cancer subtypes","Benign / non-malignant","Density")[i], adj = 0, cex = 0.44, col = INK3)
}

text(0.5, 0.988, "Genetic correlation landscape: reproductive disorders x breast outcomes",
     adj = 0.5, cex = 0.82, col = INK2)
text(0.5, 0.963, sprintf("%d exposures with estimable rg x 14 breast outcomes = %d testable pairs",
     NR, nrow(M)), adj = 0.5, cex = 0.54, col = INK3)


dev.off()
cat("Fig2a done\n")
