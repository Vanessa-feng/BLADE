setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(dplyr)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("code/blade_mfm_modified.cpp")
source("code/function.R")

# load students_rank
load(file = paste0(getwd(), "/data/students_rank_data/students_rank.RData"))
attach(students_rank)

# Peking University
Peking_select <- peking[which(peking<=300)]
n <- length(Peking_select)
dij <- list()
for(id in 1:n){
  dij[[id]] <- Peking_select[id]
}
var_name <- "Peking"

# # Tsinghua University
# Tsinghua_select <- Tsinghua[which(Tsinghua<=300)]
# n <- length(Tsinghua_select)
# dij <- list()
# for(id in 1:n){
#   dij[[id]] <- Tsinghua_select[id]
# }
# var_name <- "Tsinghua"


# ----------------------------------------------------
# init setting
Max_iteration <- 2000
sigma <- 20
K_prior <- 2

# Model
set.seed(0)
K <- K_prior
theta_init = sapply(dij, mean)
mu_init = rep(0, K)
Sigma_init = rep(1, K)
z_init = c(sample(1:K, K, replace=F), sample(1:K, n-K, replace=T))  # make sure there are K clusters

alpha <- n
beta <- sigma^2*(alpha-1)

start.time <- Sys.time()
result = runMCMC(z_init-1, theta_init, mu_init, Sigma_init, dij, G = matrix(0, n, 4), f=0,
                 tau=0.1,  mu0=0, alpha = alpha, beta = beta, K_prior=K,
                 max_iters = Max_iteration, seed = 1)
end.time <- Sys.time()

blade_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

# ppm method
ppm <- gene.ppm(result$group_iter)
blade_ppm <- minbinder(ppm, method = "comp")$cl
blade_numlayer <- length(unique(blade_ppm))

result$blade_ppm <- blade_ppm
save(result, file= paste0(getwd(), "/results/students_rank_output/result_", var_name, ".RData"))