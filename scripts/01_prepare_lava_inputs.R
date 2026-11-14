#!/usr/bin/env Rscript
# Convert LDSC-format sumstats to LAVA .sum.txt, then build input.info.txt
suppressMessages(library(data.table))
IN  <- file.path(ROOT, "ldsc_40traits")
ROOT<- file.path(ROOT, "repro_breast_L3/lava_run")

TRAITS <- data.frame(
  acc  = c("GCST90454221","GCST90454215","GCST90454213","GCST90454214"),
  name = c("Menorrhagia","Ovarian_cyst","Endometriosis","Pelvic_organ_prolapse"),
  cases= c(60359, 32208, 19588, 21385),
  ctrl = c(189515, 199850, 213669, 213236),
  stringsAsFactors = FALSE)

for (i in seq_len(nrow(TRAITS))) {
  f <- file.path(IN, paste0(TRAITS$acc[i], ".sumstats.gz"))
  d <- fread(f, data.table = FALSE)
  d <- d[, 1:7]; names(d) <- c("SNP","A1","A2","eaf","Z","P","N")
  d$se   <- 1 / sqrt(2 * pmin(d$eaf,1-d$eaf) * (1 - pmin(d$eaf,1-d$eaf)) * (d$N + d$Z^2))
  d$BETA <- d$Z * d$se
  d <- d[!is.na(d$BETA) & !is.na(d$se) & d$se > 0 & !is.na(d$P), ]
  d$A1 <- toupper(d$A1); d$A2 <- toupper(d$A2)
  out <- d[, c("SNP","A1","A2","N","BETA","se","P","Z")]
  names(out)[6] <- "SE"
  fwrite(out, file.path(ROOT, paste0(TRAITS$name[i], ".sum.txt")), sep = "\t")
  cat(sprintf("%-24s %d SNPs\n", TRAITS$name[i], nrow(out)))
}

# input.info.txt: 4 new traits + BCAC_LumA (existing)
info <- rbind(
  data.frame(phenotype = TRAITS$name, cases = TRAITS$cases, controls = TRAITS$ctrl,
             filename = file.path(ROOT, paste0(TRAITS$name, ".sum.txt")),
             stringsAsFactors = FALSE),
  data.frame(phenotype = "BCAC_LumA", cases = 45253, controls = 91477,
             filename = file.path(ROOT, "BCAC_LumA.sum.txt"), stringsAsFactors = FALSE))
fwrite(info, file.path(ROOT, "input.info2.txt"), sep = "\t")
cat("\ninput.info2.txt written with", nrow(info), "phenotypes\n")
print(info[, c("phenotype","cases","controls")])
