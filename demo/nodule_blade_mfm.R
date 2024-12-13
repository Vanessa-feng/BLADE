file <- 1
load(paste0("Code_Xin/Distance_info/distance_nodule_", file, ".RData"))
components <- read.csv(paste0("Code_Xin/Cell_info/nodule_", file, ".csv"))
ref_points <- read.csv(paste0("Code_Xin/Boundary_info/B_nodule_", file, ".csv"))
Max_iteration <- 2000
Burnin <- Max_iteration/2+1
n <- nrow(components)
R <- 2 
dij <- lapply(dij, function(df) df[, R])  # one reference curve
sigma <- 15
K_prior <- 20
K <- K_prior
dij_max <- sapply(dij, max)
dij_min <- sapply(dij, min)
lambda_est <- dij_max-dij_min
lambda = lambda_est * 1.01

theta_init = sapply(dij, mean)
a1 = (theta_init - sapply(dij, min)) / lambda
b1 = 1- (sapply(dij, max) - theta_init) /lambda
r_init = (a1 + b1)/2
total_init = runif(n, min = 20, max = 40)
alpha_init = total_init * r_init
beta_init = total_init * (1 - r_init)

mu_init = rep(0, K)
Sigma_init = rep(1, K)
z_init = c(sample(1:K, K, replace=F), sample(1:K, n-K, replace=T))  # make sure there are K clusters

alpha <- n
beta <- sigma^2*(alpha-1)

### run MCMC
start.time <- Sys.time()
result = runMCMC(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                 lambda=lambda, dij, G = matrix(0, n, 4), f=0,
                 tau=0.1,  mu0=0, alpha = alpha, beta = beta, K_prior=K,
                 max_iters = Max_iteration, seed = 1)
end.time <- Sys.time()
blade_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)

ppm <- gene.ppm(result$group_iter)
blade_ppm <- minbinder(ppm, method = "comp")$cl
blade_numlayer <- length(unique(blade_ppm))

# --------------------------
# plot the layer detection result
range_dij$res <- blade_ppm[order(mean_dij)]
original_labels <- unique(range_dij$res)
new_labels <- 1:length(original_labels) 
range_dij$res <- match(range_dij$res, original_labels)

## add step func
mean_dij_per_group <- range_dij %>%
  group_by(res) %>%
  summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
            max_re_id = max(re_id))

