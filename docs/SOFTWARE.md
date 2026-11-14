# Software and versions

| Software | Version | Used for |
|----------|---------|----------|
| R | 4.5.0 | All statistical analyses and figures |
| Python | 3.11 | Data handling and figure assembly |
| LDSC | 1.0.1 | Cross-trait LD score regression |
| LAVA | 0.1.1 | Local genetic correlation |
| TwoSampleMR | 0.5.7 | Two-sample Mendelian randomisation |
| MVMR | 0.3 | Multivariable Mendelian randomisation |
| coloc | 5.2.3 | Colocalisation (`coloc.abf`) |
| moloc | — | Multiple-trait colocalisation |
| SMR | 1.4.2 | Summary-data-based Mendelian randomisation |
| data.table | 1.18.4 | Data handling |
| ggplot2 | 4.0.3 | Figures |
| locuszoomr | — | Regional association plots |

## Figure conventions

All figures were produced as vector PDF through `cairo_pdf` with the Arial family,
without gradients, drop shadows or gridlines. Composite figures were assembled at a
fixed width of 170 mm. Panel letters are upper case and placed at the upper left of
each panel.

## Reproducibility

Random seed was fixed before all resampling and clumping steps. Analysis scripts
reproducing every table and figure are archived with this repository.
