# Shared genetic architecture between reproductive disorders and breast disease

Analysis code accompanying the manuscript:

> **Shared genetic architecture links reproductive disorders with breast disease**
> Jiahan Dong, Hanghang Chen

This repository contains the scripts used to generate every figure and table in the
manuscript. It is intended as a record of the analysis, not as a turnkey pipeline:
input files are third-party summary statistics that must be obtained from the sources
listed in `data/DATA_SOURCES.md`.

---

## What the study does

Forty female reproductive traits are screened against fourteen breast outcomes using
genome-wide and local genetic correlation, the shared signal is localised to individual
linkage-disequilibrium blocks, and the gene and cell type carrying it are then resolved
by colocalisation against expression, methylation and chromatin-accessibility QTLs.
Three candidate mediating routes are tested and excluded.

The analysis proceeds in five layers, and each layer answers a different question:

| Layer | Question | Method |
|-------|----------|--------|
| 1 | Is there shared architecture at all? | Cross-trait LD score regression |
| 2 | Is the sharing localised or diffuse? | Local genetic correlation (LAVA) |
| 3 | Is the relationship causal? | Two-sample and multivariable Mendelian randomisation |
| 4 | Which gene and which cell type carry it? | Colocalisation against eQTL, mQTL, caQTL |
| 5 | Does a peripheral mediator transmit it? | Mediation screen across 1,068 candidates |

---

## Repository layout

```
scripts/
  00_config.R                      paths, package loading, shared plotting theme
  01_prepare_lava_inputs.R         builds LAVA input files from summary statistics
  02_aggregate_lava.R              aggregates local correlation across LD blocks
  03_colocalisation.R              colocalisation against molecular QTL resources
  figures/
    fig1_study_design.R            study design schematic
    fig2a_genetic_correlation_heatmap.R
    fig2bcd_genetic_correlation_panels.R
    fig3_local_correlation_panels.R
    fig4_xbp1_attribution_panels.R
    fig5_subtype_panels.R
    fig6_density_panels.R
    fig6e_benign_mr.R
    figs1_mediator_screen.R
data/
  DATA_SOURCES.md                  accessions and download locations
docs/
  SOFTWARE.md                      software and package versions
```

---

## Requirements

R 4.5.0 and Python 3.11. Package versions are listed in `docs/SOFTWARE.md`.
Figures are written as vector PDFs through `cairo_pdf`.

`scripts/00_config.R` defines the shared paths, colour palette and PDF helper used by
the figure scripts:

```r
ROOT     <- normalizePath(Sys.getenv("REPRO_BREAST_ROOT", "."), mustWork = FALSE)
OUT_DIR  <- file.path(ROOT, "output")
```

All data and output paths inside the scripts are expressed relative to `ROOT`.
Either export `REPRO_BREAST_ROOT` before running, or change the single `ROOT` line
at the top of each script.

---

## Running the analysis

Scripts are numbered in the order they were run. Layers 1 to 3 depend on
LDSC, LAVA and the `TwoSampleMR` / `MVMR` packages; layers 4 and 5 additionally
require the molecular QTL summary statistics listed in `data/DATA_SOURCES.md`.

Figure scripts are independent of one another once the upstream result tables exist,
and each writes a single PDF to `OUT_DIR`.

---

## Data

No individual-level data were generated or analysed. All inputs are publicly available
summary-level statistics; accessions and access conditions are given in
`data/DATA_SOURCES.md`.

---

## Citation

If you use this code, please cite the manuscript above.

## Contact

Hanghang Chen — 461652608@qq.com
Department of Breast Surgery, The First Affiliated Hospital of Henan University of
Traditional Chinese Medicine, Zhengzhou, China
