ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

#!/usr/bin/env Rscript
# Figure 6e — full-resolution MR of reproductive exposures on benign/density outcomes
suppressMessages({library(data.table); library(ggplot2)})

H <- ROOT
d <- fread(file.path(H, "mediator_analysis", "out", "tsmr_benign_density.csv"))

INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
POS <- "#8E6FA8"; NEG <- "#5B8C8A"; AMB <- "#B08A5E"

# for each outcome, the strongest exposures, with instrument counts annotated
sel <- d[p < 0.20]
sel[, lab := sprintf("%s -> %s  (nIV=%d)", gsub("_", " ", repro),
                     gsub("Density_", "", outcome), nIV)]
sel[, lab := factor(lab, levels = rev(lab[order(b)]))]
sel[, sig := fdr < 0.10]

p <- ggplot(sel, aes(b, lab, colour = sig)) +
  geom_vline(xintercept = 0, colour = INK3, linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(aes(size = nIV), alpha = 0.9) +
  scale_size_continuous(range = c(1.2, 3.4), name = "Instruments") +
  scale_colour_manual(values = c(`FALSE` = AMB, `TRUE` = POS),
                      labels = c("FDR > 0.10", "FDR < 0.10"), name = NULL) +
  labs(x = "Causal effect on benign breast disease or density (95% CI)", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 8.5) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6.8, colour = INK2),
        axis.line = element_line(colour = INK3, linewidth = 0.3),
        axis.ticks = element_line(colour = INK3, linewidth = 0.3))

cairo_pdf(file.path(H, "Fig6e_benign_mr_full.pdf"), width = 7.6, height = 4.6, family = "Arial")
print(p); dev.off()
cat("Fig6e written:", nrow(sel), "pairs
")
print(sel[, .(repro, outcome, nIV, b = round(b, 3), p = signif(p, 2), fdr = round(fdr, 3))])
