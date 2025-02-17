setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(imager)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("code/blade_mfm.cpp")
source("code/function.R")

# ------------------------------------------------------------------------------
# load distance file
directory_distance <- paste0(getwd(), "/data/epoc_pathology_image_data/epoc_distance/")
file <- "001_part_1_patch_2"
load(file = paste0(directory_distance, "distance_", file , ".RData"))

# Obtain the selected reference curve id (1 or 2)
# # There are two boundaries, and we choose the inner boundary as the reference curve. 
# # The `select_ref_id` are prepared in the table before the model
df_output <- read.csv(paste0(getwd(), "/results/epoc_pathology_image_output/epoc_output.csv"))
fill_index <- which(df_output$file == file)
if(is.na(df_output[fill_index, "select_ref_id"])){
  cat("select_ref_id=None.\n")
  next
}

R <- df_output[fill_index, "select_ref_id"]
dij_R <- lapply(dij, function(df) df[, R])

# ------------------------------------------------------------------------------
# BLADE Model
Max_iteration <- 2000
K_prior <- 20
sigma <- 15

set.seed(0)
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

### run MCMC
start.time <- Sys.time()
result <- runMCMC(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                  lambda=lambda, dij_R, G = matrix(0, n, 4), f=0,
                  tau=0.1,  mu0=0, alpha = alpha, beta = beta,  K_prior=K, 
                  max_iters = Max_iteration, seed = 1)
end.time <- Sys.time()

time_store <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)
ppm <- gene.ppm(result$group_iter)
z_ppm <- minbinder(ppm, method = "comp")$cl
num_layer <- length(unique(z_ppm))
result$z_ppm <- z_ppm

df_output$num_cell[fill_index] <- n
df_output$numlayer_output[fill_index] <- num_layer

write.csv(df_output, file = paste0(getwd(), "/results/epoc_pathology_image_output/epoc_output.csv"))
save(result, file = paste0(getwd(), "/results/epoc_pathology_image_output/result_", file , ".RData"))
