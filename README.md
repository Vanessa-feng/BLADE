# BLADE

This repository contains the code used for the BLADE project, a Bayesian approach for estimating epithelial layer structure from pathology image data. The examples here are meant to accompany the manuscript and to make the main analyses easier to rerun on the included processed data.

BLADE works with cell-level distance measurements from a hand-drawn or otherwise selected reference curve. The main model is implemented in C++ through `Rcpp`/`RcppArmadillo`, with R scripts for running examples and reproducing figures.

## Repository Layout

- `code/`
  - `blade.cpp`: main BLADE model used for pathology image and nodule examples.
  - `blade_modified.cpp`: simplified version used for the students-rank example.
  - `blade_dpmm.cpp`: Dirichlet process mixture model comparison.
  - `distance.cpp`: helper code for distance calculation.
  - `function.R`: plotting, posterior similarity matrix, and evaluation helpers.
  - `setup.R`: finds the project root and loads required packages.
- `data/`
  - `epoc_pathology_image_data/`: processed EPOC pathology image inputs used by the examples.
  - `nodule_data/`: processed nodule example data.
  - `students_rank_data/`: small example used to illustrate the model on one-dimensional ranked data.
  - `simulated_data/`: simulation scripts and example simulation settings.
- `demo/`
  - `epoc_blade.R`: runs BLADE on one processed EPOC pathology image patch.
  - `nodule_blade.R`: runs BLADE on one nodule example.
  - `students_rank.R`: runs the simplified BLADE model on the students-rank example.
- `reproduce/`
  - scripts for reproducing the figures and table from the manuscript using the included processed outputs.
- `results/`
  - saved example outputs used by the reproduction scripts.

## Requirements

The code was checked with R 4.5.3 on macOS. It should also run on recent R 4.x versions with a working C++ compiler.

Required R packages:

```r
install.packages(c(
  "RColorBrewer", "ggplot2", "imager", "dplyr", "Rcpp", "RcppArmadillo",
  "mcclust", "truncnorm", "clue", "cluster", "mclust", "ggpubr",
  "survival", "survminer", "reshape2", "readxl", "gridExtra", "dbscan"
))
```

If you do not want to install packages into the system R library, you can use a project-local library:

```r
dir.create(".Rlib", showWarnings = FALSE)
install.packages("clue", lib = ".Rlib", repos = "https://cloud.r-project.org")
```

`code/setup.R` automatically uses `.Rlib` when that folder exists.

## Running the Demos

Run commands from the repository root:

```bash
Rscript demo/epoc_blade.R
Rscript demo/nodule_blade.R
Rscript demo/students_rank.R
```

The scripts write `.RData` outputs into `results/`. The MCMC examples use 2000 iterations by default, so runtime depends on the machine and the size of the example.

## Reproducing Figures and Table

The reproduction scripts use the processed data and saved outputs in `results/`:

```bash
Rscript reproduce/figure_1.R
Rscript reproduce/figure_4.R
Rscript reproduce/figure_5.R
Rscript reproduce/figure_6.R
Rscript reproduce/figure_7.R
Rscript reproduce/figure_S2.R
Rscript reproduce/figure_S3.R
Rscript reproduce/table_1.R
```

Most scripts print plots to the active R graphics device. The `ggsave()` lines are left commented in several files so paths and image sizes can be adjusted before exporting final manuscript figures.

## Notes on Data

The repository includes processed example inputs and saved outputs for reproducibility. Raw clinical images and full clinical records may be subject to data-use restrictions and are not included here. The EPOC examples therefore start from precomputed cell/reference-curve distance data.

## Quick Check

After installing dependencies, this is a useful smoke test:

```bash
Rscript -e "source('code/setup.R'); source('code/function.R'); cat('BLADE setup OK\n')"
Rscript demo/epoc_blade.R
```

If `sourceCpp()` fails, check that R can find `Rcpp`, `RcppArmadillo`, and a working C++ compiler.
