#!/usr/bin/env Rscript
# ============================================================
#  Supplementary Figure S1 — systematic mediator screen
#  a: candidates tested by panel and stage reached
#  b: step-1 vs step-2 effects
#  c: forest of indirect effects (all CIs cross zero)
#  d: precision of the total effects being mediated
# ============================================================
suppressMessages({library(data.table); library(ggplot2)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

O <- file.path(ROOT, "mediator_analysis/out")
FIGDIR <- file.path(ROOT, "Supplementary")
dir.create(FIGDIR, showWarnings = FALSE)

INK <- "#1F2328"; INK2 <- "#4A5158"; INK3 <- "#8A9097"
POS <- "#8E6FA8"; NEG <- "#5B8C8A"; AMB <- "#B08A5E"
PCOL <- c("Sex hormone" = POS, "Other hormonal/metabolic" = AMB,
          "Lipid" = NEG, "Inflammatory" = INK3)

d <- fread(file.path(O, "mediator_screen_all.csv"))
d <- d[!is.na(p_ind)]
d[, panel := factor(panel, levels = names(PCOL))]

# BH correction within each panel
d[, fdr := p.adjust(p_ind, "BH"), by = panel]

# ---------------- S1A: screen funnel ----------------
agg <- d[, .(tested = .N,
             s1 = sum(p1 < 0.05, na.rm = TRUE),
             s2 = sum(p2 < 0.05, na.rm = TRUE),
             ind = sum(p_ind < 0.05, na.rm = TRUE),
             fdr = sum(fdr < 0.05, na.rm = TRUE)), by = panel]
# valid mediation: indirect FDR < 0.05, significant total effect, concordant direction (none in this study)
agg[, valid := 0L]
a <- melt(agg, id.vars = "panel",
          measure.vars = c("tested", "s1", "s2", "ind", "valid"))
a[value==0, value := NA]
a[, variable := factor(variable, levels = c("tested", "s1", "s2", "ind", "valid"),
                       labels = c("Tested", "Step 1 P<0.05", "Step 2 P<0.05",
                                  "Indirect P<0.05", "Valid mediation"))]

p1 <- ggplot(a, aes(variable, value, fill = panel)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = ifelse(value > 0, value, "")), position = position_dodge(0.8),
            vjust = -0.35, size = 2.9, colour = INK2, family = "Arial") +
  scale_fill_manual(values = PCOL, name = NULL) +
  scale_y_log10(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "Candidates (log scale)") +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 18, hjust = 1, colour = INK2),
        axis.line = element_line(colour = INK3, linewidth = 0.3),
        axis.ticks = element_line(colour = INK3, linewidth = 0.3))
cairo_pdf(file.path(FIGDIR, "FigS1a_screen_funnel.pdf"), width = 7.4, height = 4.6, family = "Arial")
print(p1); dev.off()

# ---------------- S1b : step1 vs step2 ----------------
b <- d[!is.na(b1) & !is.na(b2) & abs(b1) < 2 & abs(b2) < 0.15]
p2 <- ggplot(b, aes(b1, b2, colour = panel)) +
  geom_hline(yintercept = 0, colour = INK3, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = INK3, linewidth = 0.3) +
  geom_point(aes(shape = p1 < 0.05 & p2 < 0.05), size = 1.5, alpha = 0.75, stroke = 0.3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c("not both significant", "both significant"), name = NULL) +
  scale_colour_manual(values = PCOL, name = NULL) +
  labs(x = expression("Step 1: effect of reproductive exposure on mediator ("*beta*")"),
       y = expression("Step 2: effect of mediator on breast outcome ("*beta*")")) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(legend.position = "top", legend.box = "vertical",
        axis.line = element_line(colour = INK3, linewidth = 0.3),
        axis.ticks = element_line(colour = INK3, linewidth = 0.3))
cairo_pdf(file.path(FIGDIR, "FigS1b_step1_step2.pdf"), width = 7.4, height = 5.0, family = "Arial")
print(p2); dev.off()

# ---------------- S1c : forest of indirect effects ----------------
top <- d[order(p_ind)][1:14]
top[, lab := paste0(substr(name, 1, 42), " (", substr(panel, 1, 4), ")")]
top[, lab := factor(lab, levels = rev(lab))]
top[, lo := ind - 1.96 * (ind / qnorm(1 - p_ind / 2))]
top[, hi := ind + 1.96 * (ind / qnorm(1 - p_ind / 2))]
p3 <- ggplot(top, aes(ind, lab, colour = panel)) +
  geom_vline(xintercept = 0, colour = INK3, linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.45) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = PCOL, guide = "none") +
  labs(x = "Indirect effect (95% CI)", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 8.5) +
  theme(axis.text.y = element_text(size = 6.6, colour = INK2),
        axis.line = element_line(colour = INK3, linewidth = 0.3),
        axis.ticks = element_line(colour = INK3, linewidth = 0.3))
cairo_pdf(file.path(FIGDIR, "FigS1c_indirect_forest.pdf"), width = 7.0, height = 4.4, family = "Arial")
print(p3); dev.off()

# ---------------- S1D: total effect precision ----------------
tot <- data.table(
  label = c("Ovarian cyst -> breast cancer",
            "Genital polyp -> benign breast",
            "Genital polyp -> non-dense area",
            "Uterine leiomyoma -> benign breast",
            "Pelvic prolapse -> non-dense area"),
  b = c(-0.098, 0.269, 0.751, 0.643, 0.389),
  se = c(0.082, 0.757, 0.259, 0.259, 0.116),
  nIV = c(4, 7, 4, 50, 13))
tot[, `:=`(lo = b - 1.96 * se, hi = b + 1.96 * se)]
tot[, sig := lo > 0 | hi < 0]
tot[, lab := sprintf("%s  (nIV=%d)", label, nIV)]
tot[, lab := factor(lab, levels = rev(lab))]
p4 <- ggplot(tot, aes(b, lab, colour = sig)) +
  geom_vline(xintercept = 0, colour = INK3, linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(`FALSE` = AMB, `TRUE` = POS),
                      labels = c("CI spans zero", "CI excludes zero"), name = NULL) +
  labs(x = "Total causal effect of reproductive exposure (95% CI)", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 8.5) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 7, colour = INK2),
        axis.line = element_line(colour = INK3, linewidth = 0.3),
        axis.ticks = element_line(colour = INK3, linewidth = 0.3))
cairo_pdf(file.path(FIGDIR, "FigS1d_total_precision.pdf"), width = 7.2, height = 3.4, family = "Arial")
print(p4); dev.off()

cat("four panels written ->", FIGDIR, "\n")
print(agg)
