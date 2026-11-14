# Analysis paths
#
# Set ROOT to the location of the project on your machine.
# Every script in this repository sources this file rather than hard-coding paths.

ROOT     <- normalizePath(Sys.getenv("REPRO_BREAST_ROOT", "."), mustWork = FALSE)
DATA_DIR <- file.path(ROOT, "data")
OUT_DIR  <- file.path(ROOT, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures")

for (d in c(OUT_DIR, FIG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- shared plotting parameters -------------------------------------------
INK  <- "#1F2328"   # primary text and axis lines
INK2 <- "#4A5158"   # secondary text
INK3 <- "#8A9097"   # grid and ticks
POS  <- "#8E6FA8"   # positive association
NEG  <- "#5B8C8A"   # negative association
AMB  <- "#B08A5E"   # ambiguous

FONT <- "Arial"
BASE_SIZE <- 7

save_pdf <- function(plot, file, width, height) {
  cairo_pdf(file.path(FIG_DIR, file), width = width, height = height, family = FONT)
  print(plot)
  dev.off()
}
