# Data sources

All data used in this study are publicly available summary-level statistics.
No individual-level data were generated or analysed.

## Primary genome-wide association summary statistics

| Resource | Used for | Release | Access |
|----------|----------|---------|--------|
| FinnGen | Female reproductive traits (exposures) | Release 11 | https://r11.finngen.fi |
| FinnGen | Benign breast disease, non-dense area, pelvic organ prolapse | Release 11 | https://r11.finngen.fi |
| Breast Cancer Association Consortium (BCAC) | Breast cancer overall, ER-negative, ER-positive, and five intrinsic subtypes | 2020 release | http://bcac.ccge.medschl.cam.ac.uk |
| UK Biobank | Mammographic density traits, selected reproductive traits | — | https://www.ukbiobank.ac.uk (approved application required) |

## Regulatory and single-cell reference data

| Resource | Used for | Access |
|----------|----------|--------|
| GTEx | Breast mammary tissue eQTL, version 8 (primary) and 11 (re-annotation check) | https://gtexportal.org |
| eQTLGen | Whole-blood eQTL, large-sample instrument | https://www.eqtlgen.org |
| GoDMC | Blood methylation QTL meta-analysis | http://mqtldb.godmc.org.uk |
| TenK10K | Single-cell eQTL across immune cell types | https://doi.org/10.1038/s41586-026-10721-2 |
| Hormone Cell Atlas | Breast and ovary cell-type expression | https://hormonecellatlas.cellgeni.sanger.ac.uk |
| FinnGen multiome | Single-nucleus chromatin accessibility QTL across immune cell types | https://doi.org/10.1038/s41586-026-11078-2 |
| METABRIC | Tumour expression by PAM50 subtype | https://www.cbioportal.org |

## Derived tables

All derived results underlying the figures and tables of the manuscript are provided
as Supplementary Tables S1 to S13 in the published article. The full list of data
sources, with the accession identifier used for each trait, is given in
Supplementary Table S1.

## Notes

Summary statistics were harmonised to GRCh37 where required before local correlation
and colocalisation analyses. Allele harmonisation and palindromic variant handling
follow the default behaviour of the `TwoSampleMR` and `coloc` packages; variants with
ambiguous strand and intermediate allele frequency were excluded.
