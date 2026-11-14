#!/usr/bin/env Rscript
# Step 5: coloc.abf (duplicates removed; corrected varbeta)
Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")
suppressMessages({library(coloc); library(data.table)})

ROOT <- Sys.getenv("REPRO_BREAST_ROOT", ".")   # project root; set this or export the variable

BASE <- file.path(ROOT, "repro_breast_L3/coloc")
files <- list.files(BASE, pattern = "^chr.*\\.tsv$", full.names = TRUE)
S_PCOS <- 5171 / 288356; N_PCOS <- 288356
S_LUMA <- 45253 / 136730; N_LUMA <- 136730

# derive beta/se from Z, MAF and N: se = 1/sqrt(2*MAF*(1-MAF)*(N + Z^2))
z2bs <- function(z, maf, N) {
  se <- 1 / sqrt(2 * maf * (1 - maf) * (N + z^2))
  list(beta = z * se, varbeta = se^2)
}

res <- list()
for (f in files) {
  reg <- sub("\\.tsv$", "", basename(f))
  d <- fread(f, data.table = FALSE)
  d <- d[!is.na(d$MAF) & d$MAF > 0.01, ]
  d <- d[!duplicated(d$snp), ]
  if (nrow(d) < 50) { cat(sprintf("%-18s too few SNPs (%d)\n", reg, nrow(d))); next }

  a <- z2bs(d$z_pcos, d$MAF, N_PCOS)
  b <- z2bs(d$z_luma, d$MAF, N_LUMA)
  cat(sprintf("%-18s nSNP=%4d  PCOS min P=%.2e  LumA min P=%.2e\n",
              reg, nrow(d), min(d$p_pcos), min(d$p_luma)))

  r <- try(coloc.abf(
    dataset1 = list(snp = d$snp, beta = a$beta, varbeta = a$varbeta,
                    MAF = d$MAF, N = N_PCOS, type = "cc", s = S_PCOS),
    dataset2 = list(snp = d$snp, beta = b$beta, varbeta = b$varbeta,
                    MAF = d$MAF, N = N_LUMA, type = "cc", s = S_LUMA)), silent = TRUE)
  if (inherits(r, "try-error")) { cat("   -> failed\n"); next }

  s <- r$summary
  top <- r$results[which.max(r$results$SNP.PP.H4), ]
  res[[reg]] <- data.frame(
    region = reg, n_snp = as.integer(s[["nsnps"]]),
    PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
    PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
    top_snp = top$snp, top_snp_pp = top$SNP.PP.H4, stringsAsFactors = FALSE)
}

out <- do.call(rbind, res)
out <- out[order(-out$PP.H4), ]
fwrite(out, file.path(BASE, "coloc_results.csv"))

cat("\n==== coloc.abf：PCOS vs BCAC Luminal A ====\n")
cat(sprintf("%-18s %6s %8s %8s %8s %8s   %-13s %7s\n",
            "region", "nSNP", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "top SNP", "SNP.PP"))
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  flag <- if (r$PP.H4 > 0.8) "  <<< shared causal variant" else
          if (r$PP.H3 > 0.8) "  <<< distinct causal variants in the same region" else ""
  cat(sprintf("%-18s %6d %8.3f %8.3f %8.3f %8.3f   %-13s %7.3f%s\n",
              r$region, r$n_snp, r$PP.H0, r$PP.H1, r$PP.H2, r$PP.H3,
              r$top_snp, r$top_snp_pp, flag))
}
cat("\nH4 = shared causal variant; H3 = distinct causal variants in the same region; H0/H1/H2 = signal in one trait only\n")
