
output_df <- data.frame(matrix(nrow=Repetition, ncol=47))
colnames(output_df) <- c("repetition", "num_layer", "num_nuclei", "sd_within", "sd_between", 
                         "blade_numlayer", "kmeans_numlayer", "gmm_numlayer", "gap_numlayer", "dbscan_numlayer", "dpmm_numlayer",
                         "blade_ari", "kmeans_ari", "gmm_ari", "gap_ari", "dbscan_ari", "dpmm_ari",
                         "blade_nmi", "kmeans_nmi", "gmm_nmi", "gap_nmi", "dbscan_nmi", "dpmm_nmi",
                         "blade_ami", "kmeans_ami", "gmm_ami", "gap_ami", "dbscan_ami", "dpmm_ami",
                         "blade_fmi", "kmeans_fmi", "gmm_fmi", "gap_fmi", "dbscan_fmi", "dpmm_fmi",
                         "blade_cor", "kmeans_cor", "gmm_cor", "gap_cor", "dbscan_cor", "dpmm_cor",
                         "blade_time", "kmeans_time", "gmm_time", "gap_time", "dbscan_time", "dpmm_time")

for(rep in 1:Repetition){
  num_layer <- round(rtruncnorm(1, a=5, b=35, mean=15, sd=10))
  num_nuclei <- round(rtruncnorm(num_layer, a=20, b=200, mean=100, sd=60))
  sd_within <- sapply(num_nuclei, generate_sd_within, num_nuclei)
  sd_between <- rtruncnorm(num_layer, a=35, b=105, mean=60, sd=10)
  n <- sum(num_nuclei)  # total number of cells
  
  ### Generate data
  cat(rep, "\t")
  ni <- sample(90:150, n, replace=TRUE)
  
  # --------------------------------------------
  # Generate sample dij
  # dij: list with length n=300 
  # x[[j]]: cell 1, array with length ni[j]
  # z_true: true layer
  # --------------------------------------------
  dij <- list()
  z_true <- c()
  theta_true <- c()
  
  h_k <- 0
  for(k in 1:num_layer){
    z_true <- c(z_true, rep(k,num_nuclei[k]))
    h_k <- h_k + sd_between[k]
    theta_true <- c(theta_true, rnorm(num_nuclei[k], h_k, sd_within))
  }
  
  lambda_true = runif(n, 20, 40)
  alpha_true = rgamma(n, 1, 0.5)+ 0.4
  beta_true = rgamma(n, 1, 0.5) + 0.4
  
  # generate cells
  for (i in 1:n){
    cell = rbeta(ni[i]*proportion, alpha_true[i], beta_true[i])*lambda_true[i] +
      theta_true[i]-alpha_true[i]/(alpha_true[i]+beta_true[i])*lambda_true[i]
    dij[[i]] = cell
  }
  
  save(dij, file = paste0("Code_Xin\\Simulation\\Simulation_comparison\\store_result\\distance_", rep , ".RData"))
  # -----------------------------------------
  # plot true layers with color palette
  min_dij = sapply(dij, min)
  max_dij = sapply(dij, max)
  mean_dij = sapply(dij, mean)
  range_dij = cbind(min_dij, max_dij, mean_dij)
  range_dij = as.data.frame(range_dij)
  range_dij = range_dij[order(mean_dij), ]
  range_dij$re_id = rep(1:n)
  
  range_dij$res <- z_true[order(mean_dij)]
  original_labels <- unique(range_dij$res)
  new_labels <- 1:length(original_labels) 
  range_dij$res <- match(range_dij$res, original_labels)
  
  # Create the plot with the true layers
  p_true <- ggplot() +
    geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                                 xend = re_id, yend = max_dij, color = factor(res)),
                 linewidth=1.2) +
    scale_color_manual(values = my_colors) +
    labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
    ggtitle("True layers")
  ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\true_layer\\z_true_", rep , ".png"), p_true,
         width=2200, height=1500, units = "px")
  
  # -------------------------------------------
  ## BLADE
  #### init
  K <- K_prior
  lambda = lambda_true 
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
  
  result$num_layer <- num_layer
  result$num_nuclei <- num_nuclei
  result$sd_within <- sd_within
  result$sd_between <- sd_between
  result$lambda_true <- lambda_true
  result$z_true <- z_true
  
  # ppm method
  ppm <- gene.ppm(result$group_iter)
  result$blade_ppm <- minbinder(ppm, method = "comp")$cl
  
  blade_numlayer <- length(unique(result$blade_ppm))
  blade_ari <- adjustedRandIndex(result$blade_ppm, z_true)
  result$blade_ari <- blade_ari
  
  blade_nmi <- NMI(result$z_true, result$blade_ppm)
  blade_ami <- AMI(result$z_true, result$blade_ppm)
  blade_fmi <- fmi(result$z_true, result$blade_ppm)
  blade_cor <- cor.spearman(result$z_true, result$blade_ppm)
  
  # --------------------------
  # plot the layer detection result
  range_dij$res <- result$blade_ppm[order(mean_dij)]
  original_labels <- unique(range_dij$res)
  new_labels <- 1:length(original_labels) 
  range_dij$res <- match(range_dij$res, original_labels)
  
  ## add step func
  mean_dij_per_group <- range_dij %>%
    group_by(res) %>%
    summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
              max_re_id = max(re_id))
  
  # Create the plot with the modified color palette
  if(blade_numlayer <=36){
    p_blade <- ggplot() +
      geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                                   xend = re_id, yend = max_dij, color = factor(res))) +
      scale_color_manual(values = my_colors) +
      labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
      geom_segment(data = mean_dij_per_group, 
                   aes(x = min_re_id, xend = max_re_id, 
                       y = mean_dij, yend = mean_dij, color = factor(res)),
                   linewidth = 1.5)+ 
      ggtitle("BLADE")
  }else{
    my_colors_more <- rep(my_colors, ceiling(blade_numlayer/36))
    p_blade <- ggplot() +
      geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                                   xend = re_id, yend = max_dij, color = factor(res))) +
      scale_color_manual(values = my_colors_more) +
      labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
      geom_segment(data = mean_dij_per_group, 
                   aes(x = min_re_id, xend = max_re_id, 
                       y = mean_dij, yend = mean_dij, color = factor(res)),
                   linewidth = 1.5)+ 
      ggtitle("BLADE")
  }
  ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\blade_result\\z_blade_", rep , ".png"), p_blade,
         width=2200, height=1500, units = "px")
  
  