#!/usr/bin/env Rscript
# Aggregate LAVA chunks: 4 reproductive traits x BCAC Luminal A
suppressMessages(library(data.table))
ROOT <- file.path(ROOT, "repro_breast_L3/lava_run")

fs <- list.files(file.path(ROOT, "chunks_new"), pattern = "^lava_new_c.*csv$", full.names = TRUE)
cat("chunks found:", length(fs), "\n")
if (!length(fs)) quit(status = 0)

A <- rbindlist(lapply(fs, fread), fill = TRUE)
cat("total rows:", nrow(A), "\n")
print(table(A$type))

b <- A[type == "bivar"]
if (!nrow(b)) { cat("no bivar results\n"); quit(status = 0) }

# keep only phenotype pairs involving BCAC_LumA
b <- b[grepl("BCAC_LumA", phen1) | grepl("BCAC_LumA", phen2)]
cat("bivar pairs with BCAC_LumA:", nrow(b), "\n")

b[, pair := paste(pmin(phen1, phen2), pmax(phen1, phen2), sep = " x ")]
agg <- b[, .(n_loci = .N,
             n_sig = sum(p < 0.05 / .N, na.rm = TRUE),
             min_p = min(p, na.rm = TRUE),
             rho_med = median(rho, na.rm = TRUE),
             rho_min = min(rho, na.rm = TRUE),
             rho_max = max(rho, na.rm = TRUE)), by = pair]
setorder(agg, min_p)
fwrite(b, file.path(ROOT, "lava_new_bivar_full.csv"))
fwrite(agg, file.path(ROOT, "lava_new_summary.csv"))

cat("\n=== summary by pair ===\n")
print(agg, digits = 3)

cat("\n=== Bonferroni-significant loci ===\n")
sig <- b[!is.na(p) & p < 0.05 / .N, ][order(p)]
if (nrow(sig)) print(sig[, .(pair, LOC, rho, rho_se, p, nsnps)], digits = 3) else cat("none\n")
