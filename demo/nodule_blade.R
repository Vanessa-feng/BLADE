setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(dplyr)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("code/blade_mfm.cpp")
source("code/function.R")

# ------------------------------------------
file <- 1
load(paste0(getwd(), "/data/nodule_data/distance_nodule_", file, ".RData"))
components <- read.csv(paste0(getwd(), "/data/nodule_data/nodule_cell_", file, ".csv"))
ref_points <- read.csv(paste0(getwd(), "/data/nodule_data/nodule_ref_", file, ".csv"))

R <- 1 # Choose the bottom reference curve, select_ref_id=1
dij <- lapply(dij, function(df) df[, R])  
n <- length(dij)

# sample 5% pixels
set.seed(1234)
dij_sample <- lapply(dij, function(x){
  sample_size <- ceiling(length(x) * 0.05)
  sample(x, size=sample_size)
})

# initial setting
set.seed(0)
Max_iteration <- 2000
Burnin <- Max_iteration/2+1
sigma <- 14
K_prior <- 5

dij_max <- sapply(dij_sample, max)
dij_min <- sapply(dij_sample, min)
lambda_est <- dij_max-dij_min
lambda = lambda_est * 1.01

theta_init = sapply(dij_sample, mean)
a1 = (theta_init - sapply(dij_sample, min)) / lambda
b1 = 1- (sapply(dij_sample, max) - theta_init) /lambda
r_init = (a1 + b1)/2
total_init = runif(n, min = 20, max = 40)
alpha_init = total_init * r_init
beta_init = total_init * (1 - r_init)

K <- K_prior
mu_init = rep(0, K)
Sigma_init = rep(1, K)
z_init = c(sample(1:K, K, replace=F), sample(1:K, n-K, replace=T))  # make sure there are K clusters

alpha <- n
beta <- sigma^2*(alpha-1)

### run MCMC
start.time <- Sys.time()
result = runMCMC(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                 lambda=lambda, dij_sample, G = matrix(0, n, 4), f=0,
                 tau=0.1,  mu0=0, alpha = alpha, beta = beta, K_prior=K,
                 max_iters = Max_iteration, seed = 1)
end.time <- Sys.time()
blade_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

ppm <- gene.ppm(result$group_iter)
blade_ppm <- minbinder(ppm, method = "comp")$cl
blade_numlayer <- length(unique(blade_ppm))

# store results
result$blade_ppm <- blade_ppm
save(result, file = paste0(getwd(), "/results/nodule_output/result_nodule_", file , ".RData"))
