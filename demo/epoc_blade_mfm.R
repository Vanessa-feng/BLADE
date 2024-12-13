library(Rcpp)
library(RcppArmadillo)
sourceCpp("Code_bencong/blade.cpp")

library(mclust)
library(mcclust)
library(dplyr)
library(kernlab) 
library(gridExtra)
library(truncnorm)

# ------------------------------------------
# BLADE Model
file <- "001_part_1_patch_2"

img_file_path <- paste0(directory_patches, "EPOC premalignant trial-imaging AI study-", file, ".png")
ref_file_path <- paste0(directory_ref, "EPOC premalignant trial-imaging AI study-", file, ".png")

if (!file.exists(img_file_path) || !file.exists(ref_file_path)) {
  cat("File does not exist, skipping:", file, "\n")
  next
}

# nuclei img
img_data <- load.image(img_file_path)
width <- dim(img_data)[1]
height <- dim(img_data)[2]

ref_file <- paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_ref\\", file, ".csv")
components_file <- paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_patches\\", file, ".csv")
if (!file.exists(ref_file)) {
  cat("File ref does not exist, skipping:", file, "\n")
  next
}
if (!file.exists(components_file)) {
  cat("File components does not exist, skipping:", file, "\n")
  next
}

ref_points <- read.csv(ref_file)
components <- read.csv(components_file)

load(file = paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_distance\\distance_", file , ".RData"))

num_ref <- length(unique(ref_points$boundary_id))
num_cell <- length(unique(components$cell_id))

# ------------------------------------------------------------

fill_index <- which(file==df_output$file)

# BLADE
if(is.na(df_output[fill_index, "select_ref_id"])){
  cat("select_ref_id=None.\n")
  next
}

R <- df_output[fill_index, "select_ref_id"]
dij_R <- lapply(dij, function(df) df[, R])

K_prior <- 20
set.seed(0)

n <-  length(dij_R)
dij_R_max <- sapply(dij_R, max)
dij_R_min <- sapply(dij_R, min)
lambda_est <- dij_R_max - dij_R_min
lambda_store[id] <- mean(lambda_est)

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

# Adaptive alpha
alpha <- n

# # Adaptive beta
sigma <- 15
beta <- sigma^2*(alpha-1)

df_output$alpha[fill_index] <- alpha
df_output$beta[fill_index] <- beta
df_output$sigma[fill_index] <- sigma
df_output$num_cell[fill_index] <- num_cell

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

df_output$running_time[fill_index] <- time_store
result$z_ppm <- z_ppm

cat("file:", file, ", R=", R, ", sigma=", sigma, ", alpha=", alpha, ", beta=", beta,
    ", numlayer_ouput=", num_layer, ", time_ouput=", time_store, "\n")


df_output$numlayer_output[fill_index] <- num_layer
save(result, file = paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\sigma_15_alpha_n\\result_", file , ".RData"))
