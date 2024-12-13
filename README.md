# Bayesian Layer Detection (BLADE)

## Introduction

Bayesian Layer Detection (BLADE) is a framework designed to infer the number of oral epithelial layers from pathology images. It aims to achieve four main objectives:

1. Performing cell segmentation and feature extraction using AI-based methods.
2. Detecting and distinguishing cell structures into layers.
3. Extracting cell shape information.
4. Analyzing the clinical relevance of epithelial layers in relation to oral cancer progression.

The model leverages spatial data of cells relative to a reference curve, capturing both distance and shape information.

## Directory Structure

- **`code`**:
  - **`blade_mfm.cpp`**: Implements the BLADE framework using a mixture of finite mixtures (MFM) approach.
  - **`blade_dpmm.cpp`**: Implements the second BLADE component using a Dirichlet Process Mixture Model (DPMM) approach for comparison.
  - **`distance.cpp`**: Processes pathology images to compute the shortest pixel-to-reference curve distances for cells.
  - **`preprocessing.R`**: Prepares spatial and type information from pathology images and calculates the pixel-to-reference curve distances, serving as input for BLADE.
  - **`function.R`**: Provides basic custom functions used in the analysis.

- **`demo`**: Contains a demo dataset for BayesLASA.

- **`data`**: Includes datasets used in the manuscript:
  - **`simulated_data`**: Contains scripts for generating simulated datasets along with 10 sample datasets.
  - **`epoc_pathology_image_data`**: Provides data related to pathology images:
    - **`slide_info`**: Includes two cropped sample slides with corresponding spatial and type data.
    - **`ref_info`**: Contains hand-drawn reference curves for the slides and their converted spatial data.
    - **`epoc_distance`**: Provides distance data calculated from slide and reference curve data, used as input for BLADE.
    - **`epoc_clinic.csv`**: Contains clinical data for all patients, including the median number of layers for each patient, derived from BLADE-inferred layer counts across all sample slides. This data is used for survival analysis.

- **`results`**: Stores and displays all output results:
  - **`simulation_output`**: Results from simulated datasets.
  - **`epoc_pathology_image_output`**: Results from pathology image datasets.

- **`reproduce`**: Includes files to reproduce the figures and tables presented in the BayesLASA manuscript.

---

This directory structure and accompanying code provide a comprehensive workflow for the application of BLADE to pathology images, facilitating both methodological advancements and clinical insights into oral cancer progression.

## Usage

Below, we demonstrate the usage of BLADE for infering the number of layers on the OPMD-EPOC dataset. The core layer detection functionality is performed by the `runMCMC` function in `code/blade_mfm.cpp`.

### Purpose

The `runMCMC` function implements the Bayesian Layer Detection framework using Markov Chain Monte Carlo (MCMC) methods. It updates parameters iteratively to infer the number of layers and associated shape parameters for cells.

### Required Arguments

- `z_init` (integer vector): Initial sample allocation vector indicating layer membership.
- `alpha_init` (numerical vector): Initial shape parameter  for each cell.
- `beta_init` (numerical vector): Initial shape parameter  for each cell.
- `theta_init` (numerical vector): Initial mean shortest distance  for each cell.
- `mu_init` (numerical vector): Initial layer-specific mean parameters.
- `Sigma_init` (numerical vector): Initial layer-specific variance parameters.
- `lambda` (numerical vector): Size parameters , corresponding to the length of each cell along the axis orthogonal to the tangent of the reference curve.
- `distance` (list): Precomputed shortest distances between cell pixels and reference curve points.
- `G` (matrix): Boolean adjacency matrix (default is a zero matrix of dimension ).
- `f` (numerical scalar): Hyperparameter controlling the prior (default = 0).
- `mu0` (numerical scalar): Prior mean for layer-specific means.
- `tau` (numerical scalar): Precision parameter for the prior on layer-specific means.
- `alpha` (numerical scalar): Shape parameter for the prior on variance.
- `beta` (numerical scalar): Rate parameter for the prior on variance.
- `K_prior` (integer): Prior on the maximum number of layers.
- `GAMMA` (numerical scalar): Concentration parameter for the Dirichlet process (default = 1).
- `max_iters` (integer): Maximum number of iterations for the MCMC algorithm (default = 100).
- `seed` (integer): Random seed for reproducibility (default = 12569).

### Function Output

- `res` (List): A list containing the following components:

  - `K_iter` (vector): Number of layers inferred at each iteration.
  - `group_iter` (matrix): Sample allocation vector  at each iteration.
  - `theta_iter` (matrix): Mean shortest distance  for each cell at each iteration.
  - `alpha_iter` (matrix): Shape parameter  for each cell at each iteration.
  - `beta_iter` (matrix): Shape parameter  for each cell at each iteration.
  - `mu_iter` (matrix): Layer-specific mean parameters  at each iteration.
  - `Sigma_iter` (matrix): Layer-specific variance parameters  at each iteration.

  ### Case study

  我们实现case study主要用real oral cancer pathology imaging data from the Erlotinib Prevention of Oral Cancer (EPOC) trial at the University of Texas MD Anderson Cancer Center, which comprises a cohort of 136 patients with corresponding histology slides. All imaging data include marked spatial cell pattern data, where each cell consists of multiple pixels.  Here, we demonstrate BLADE 用一张被提前预处理好的sample iamge的数据获取的distance

  ```R
  library(dplyr)
  library(mcclust)
  library(Rcpp)
  library(RcppArmadillo)
  
  sourceCpp("Code_bencong/blade.cpp")
  
  # PPM function
  gene.ppm <- function(z_store, iter=2000, burn=1001){
    n <- nrow(z_store)
    ppm <- matrix(0, nrow = n, ncol = n)
    for (ii in burn:iter){
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          if (z_store[i, ii] == z_store[j, ii]) {
            ppm[i, j] <- ppm[i, j] + 1
            ppm[j, i] <- ppm[j, i] + 1
          }
        }
      } 
    } 
    ppm <- ppm/(iter - burn + 1)
    diag(ppm) <- rep(1, n)
    return(ppm)
  }
  
  # load distance data
  file <- "001_part_1_patch_2"
  load(paste0("data/epoc_pathology_image_data/epoc_distance/distance_", file, ".RData")) # dij
  R <- 1 # there are two reference curves, choose the inner one as reference curve.
  dij_R <- lapply(dij, function(df) df[, R])
  
  # Setting
  Max_iteration <- 2000
  Burnin <- Max_iteration/2
  K_prior <- 20
  sigma <- 15
  set.seed(0)
  
  # Initial
  n <-  length(dij_R)
  dij_R_max <- sapply(dij_R, max)
  dij_R_min <- sapply(dij_R, min)
  lambda_est <- dij_R_max - dij_R_min
  lambda <- lambda_est * 1.01
  
  K <- K_prior
  theta_init <- sapply(dij_R, mean)
  a1 <- (theta_init - sapply(dij_R, min)) / lambda
  b1 <- 1- (sapply(dij_R, max) - theta_init) /lambda
  r_init <- (a1 + b1)/2
  total_init <- runif(n, min = 1, max = 20)
  alpha_init <- total_init * r_init
  beta_init <- total_init * (1 - r_init)
  
  z_init = c(sample(1:K, K, replace=F), sample(1:K, n-K, replace=T))
  mu_init <- rep(0, K)
  Sigma_init <- rep(1, K)
  alpha <- n
  beta <- sigma^2*(alpha-1)
  
  # BLADE
  result <- runMCMC(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                    lambda=lambda, dij_R, G = matrix(0, n, 4), f=0,
                    tau=0.1,  mu0=0, alpha = alpha, beta = beta, K_prior = K,
                      max_iters = Max_iteration, seed = 1)
  
  # Posterior Inference: PPM 
  ppm <- gene.ppm(result$group_iter)
  z_ppm <- minbinder(ppm, method = "comp")$cl
  num_layer <- length(unique(z_ppm))
  
  # -------------------------------------------
  # Plots
  set.seed(245)
  my_colors <- sample(c(brewer.pal(8,"Dark2"), brewer.pal(9,"Set1"), brewer.pal(12,"Set3"), brewer.pal(12,"Paired")), 36, replace = FALSE)
                  
   plot.label <- function(label, loc, boundary, main = "", width, height, color_palette) {
    data <- data.frame(expr = label, x = loc[, "x"], y = loc[, "y"]);
    data1 <- data.frame( bx = boundary[, "x"], by = boundary[, "y"])
    ggplot() + geom_point(data=data, mapping = aes(x = x, y = y, color = expr), size = 1) + 
      scale_color_manual(values = color_palette) +
      geom_point(data = data1, mapping = aes(x = bx, y = by), color = "black", size = 1) + 
      coord_fixed(ratio = 1) +  
      theme_classic() + labs(color = "", title = main) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust = 0.5, face = "bold", size = 20),
            panel.border = element_blank(),
            legend.text = element_text(face = "bold", size = 16),
            legend.position="bottom")+
      xlim(0, width) + ylim(0, height)
  }
  
  
  # sample img
  img_file_path <- paste0("data/epoc_pathology_image_data/slide_info/EPOC premalignant trial-imaging AI study-", file, ".png")
  
  if (!file.exists(img_file_path)) {
    stop(paste("File sample iamge does not exist, skipping:", file))
  }
  
  img_data <- load.image(img_file_path)
  width <- dim(img_data)[1]
  height <- dim(img_data)[2]
  
  # sample slides data and reference curve data
  components_file <- paste0("data/epoc_pathology_image_data/slide_info/", file, ".csv")
  ref_file <- paste0("data/epoc_pathology_image_data/ref_info/", file, ".csv")
  
  if (!file.exists(ref_file)) {
    stop(paste0("File ref does not exist, skipping:", file))
  }
  if (!file.exists(components_file)) {
    stop(paste0("File components does not exist, skipping:", file))
  }
  
  ref_points <- read.csv(ref_file)
  components <- read.csv(components_file)
  
  # re-order z_ppm according to the mean of theta
  num_cell <- length(unique(components$cell_id))
  theta_means <- rowMeans(result$theta_iter[, (Burnin + 1):2000])
  theta_group_means <- tapply(theta_means, result$z_ppm, mean)
  sorted_groups <- sort(theta_group_means, index.return = TRUE)
  result$z_ppm <- match(result$z_ppm, sorted_groups$ix)
  
  # draw the cluster result
  for(i in 1:num_cell){
    components[components$cell_id==unique(components$cell_id)[i], "label"] <- result$z_ppm[i]
  }
  
  p_blade <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                   boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                   color_palette=my_colors)
  
  # draw the color-coded layers
  min_dij = sapply(dij_R, min)
  max_dij = sapply(dij_R, max)
  mean_dij = sapply(dij_R, mean)
  range_dij = cbind(id = 1:n, min_dij, mean_dij, max_dij)
  range_dij = as.data.frame(range_dij)
  range_dij = range_dij[order(mean_dij), ]
  range_dij$re_id = rep(1:n)
  
  range_dij$res <- result$z_ppm[order(mean_dij)]
  original_labels <- unique(range_dij$res)
  new_labels <- 1:length(original_labels) 
  range_dij$res <- match(range_dij$res, original_labels)
  
  ## add step func
  mean_dij_per_group <- range_dij %>%
  					group_by(res) %>%
  					summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id), max_re_id = max(re_id))
  p_colorcoded <- ggplot() +
      geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                                   xend = re_id, yend = max_dij, color = factor(res))) +
      scale_color_manual(values = my_colors) +
      labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
      geom_segment(data = mean_dij_per_group, 
                   aes(x = min_re_id, xend = max_re_id, 
                       y = mean_dij, yend = mean_dij, color = factor(res)),
                   linewidth = 1.5)+ 
      ggtitle("BLADE")
  ```

   <img src="/Users/vanessafeng/Downloads/UTD/oral_cancer/OPMD_EPOC/data_result/sigma_15_alpha_n/label_001_part_1_patch_2.png" alt="label_001_part_1_patch_2" style="zoom:20%;" /><img src="/Users/vanessafeng/Downloads/UTD/oral_cancer/OPMD_EPOC/data_result/sigma_15_alpha_n/vp_001_part_1_patch_2.png" style="zoom:20%;" />
