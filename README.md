# xci-sex-bias

#### Repository for an analysis of XCI and sex-biased expression

This repository contains scripts used in the analysis of XCI and sex-biased expression.

I ran these analyses using [R](https://cran.r-project.org/) (v4.3).

# Inputs

The following data files are expected:

* GTEx nmXCI ASE data from Gylemo et al: gylemo.csv
* GTEx sex-biased expression from DeCasien et al: mashr_sex.rds
* GTEx sex-biased expression from Oliva et al: oliva.csv
* GTEx sex-biased enhancer activity from Hou et al: hou-et-al-h3k27ac.csv
* TADS: A-172_GSE147123_tad.bed
* ASE from single cell datasets: other-ASE-results.csv
  
# Pipeline

### all scripts required for analysis and figures

```
# Create sex-specific splici transcriptomes and references
scripts/Rscripts.R
```
