output_df <- data.frame(matrix(nrow=Repetition, ncol=47))
colnames(output_df) <- c("repetition", "num_layer", "num_nuclei", "sd_within", "sd_between", 
                         "blade_numlayer", "kmeans_numlayer", "gmm_numlayer", "gap_numlayer", "dbscan_numlayer", "dpmm_numlayer",
                         "blade_ari", "kmeans_ari", "gmm_ari", "gap_ari", "dbscan_ari", "dpmm_ari",
                         "blade_nmi", "kmeans_nmi", "gmm_nmi", "gap_nmi", "dbscan_nmi", "dpmm_nmi",
                         "blade_ami", "kmeans_ami", "gmm_ami", "gap_ami", "dbscan_ami", "dpmm_ami",
                         "blade_fmi", "kmeans_fmi", "gmm_fmi", "gap_fmi", "dbscan_fmi", "dpmm_fmi",
                         "blade_cor", "kmeans_cor", "gmm_cor", "gap_cor", "dbscan_cor", "dpmm_cor",
                         "blade_time", "kmeans_time", "gmm_time", "gap_time", "dbscan_time", "dpmm_time")

for(rep in 181:Repetition){
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
  
  
  # -------------------------------------------------------------------------
  input_d <- sapply(dij, mean)
  # --------------------------------------
  # 1. k-means + silhouette method
  k_candidate <- 5:35
  
  # Calculate silhouette scores for each k
  start.time <- Sys.time()
  silhouette_scores <- sapply(k_candidate, silhouette_score)
  kmeans_numlayer <- k_candidate[which.max(silhouette_scores)]
  end.time <- Sys.time()
  kmeans_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)
  
  # ARI 
  kmeans_output <- kmeans(input_d, centers = kmeans_numlayer, nstart = 5)
  kmeans_ari <- adjustedRandIndex(kmeans_output$cluster, z_true)
  
  result$kmeans_cluster <- kmeans_output$cluster
  result$kmeans_ari <- kmeans_ari
  
  kmeans_nmi <- NMI(result$z_true, result$kmeans_cluster)
  kmeans_ami <- AMI(result$z_true, result$kmeans_cluster)
  kmeans_fmi <- fmi(result$z_true, result$kmeans_cluster)
  kmeans_cor <- cor.spearman(result$z_true, result$kmeans_cluster)
  # --------------------------
  # plot the layer detection result
  range_dij$res <- result$kmeans_cluster[order(mean_dij)]
  original_labels <- unique(range_dij$res)
  new_labels <- 1:length(original_labels) 
  range_dij$res <- match(range_dij$res, original_labels)
  
  ## add step func
  mean_dij_per_group <- range_dij %>%
    group_by(res) %>%
    summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
              max_re_id = max(re_id))
  
  # Create the plot with the modified color palette
  p_kmeans <- ggplot() +
    geom_point(data = range_dij, mapping = aes(x = re_id, y = mean_dij, color = factor(res))) +
    scale_color_manual(values = my_colors) +
    labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
    geom_segment(data = mean_dij_per_group, 
                 aes(x = min_re_id, xend = max_re_id, 
                     y = mean_dij, yend = mean_dij, color = factor(res)),
                 linewidth = 1.2) + 
    ggtitle("Kmeans")
  ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\kmeans_result\\z_kmeans_", rep , ".png"), p_kmeans,
         width=2200, height=1500, units = "px")
  
  
  # ---------------------------------------------------------------------
  # 2. GMM + BIC method
  start.time <- Sys.time()
  gmm <- Mclust(input_d, G=5:35, modelNames="V", verbose=FALSE)
  gmm <- Mclust(input_d, G=gmm$G, n_init=5, verbose=FALSE)
  gmm_numlayer <- length(unique(gmm$classification))
  end.time <- Sys.time()
  gmm_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)
  
  # ARI
  gmm_ari <- adjustedRandIndex(gmm$classification, z_true)
  
  sorted_centers <- sort(gmm$parameters$mean)
  cluster_mapping <- match(gmm$parameters$mean, sorted_centers)
  result$gmm_cluster <- as.integer(factor(cluster_mapping[gmm$classification]))
  result$gmm_ari <- gmm_ari
  
  gmm_nmi <- NMI(result$z_true, result$gmm_cluster)
  gmm_ami <- AMI(result$z_true, result$gmm_cluster)
  gmm_fmi <- fmi(result$z_true, result$gmm_cluster)
  gmm_cor <- cor.spearman(result$z_true, result$gmm_cluster)
  # --------------------------
  # plot the layer detection result
  range_dij$res <- result$gmm_cluster[order(mean_dij)]
  original_labels <- unique(range_dij$res)
  new_labels <- 1:length(original_labels) 
  range_dij$res <- match(range_dij$res, original_labels)
  
  ## add step func
  mean_dij_per_group <- range_dij %>%
    group_by(res) %>%
    summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
              max_re_id = max(re_id))
  
  # Create the plot with the modified color palette
  p_gmm <- ggplot() +
    geom_point(data = range_dij, mapping = aes(x = re_id, y = mean_dij, color = factor(res))) +
    scale_color_manual(values = my_colors) +
    labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
    geom_segment(data = mean_dij_per_group, 
                 aes(x = min_re_id, xend = max_re_id, 
                     y = mean_dij, yend = mean_dij, color = factor(res)),
                 linewidth = 1.2) + 
    ggtitle("GMM")
  ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\gmm_result\\z_gmm_", rep , ".png"), p_gmm,
         width=2200, height=1500, units = "px")
  
  # -----------------------------------------------------------------------
  # # 3. Gap statistic -- distance-based
  k_candicate <- 5:35
  start.time <- Sys.time()
  gap_stat <- clusGap(as.matrix(input_d), FUN = kmeans, nstart = 5, K.max = 35, B = 100, verbose=F)
  gap_numlayer <- which.max(gap_stat$Tab[k_candicate, "gap"])+4
  end.time <- Sys.time()
  gap_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)
  
  # ARI
  gap_output <- kmeans(input_d, centers = gap_numlayer, nstart = 5)
  gap_ari <- adjustedRandIndex(gap_output$cluster, z_true)
  
  result$gap_cluster <- gap_output$cluster
  result$gap_ari <- gap_ari
  
  gap_nmi <- NMI(result$z_true, result$gap_cluster)
  gap_ami <- AMI(result$z_true, result$gap_cluster)
  gap_fmi <- fmi(result$z_true, result$gap_cluster)
  gap_cor <- cor.spearman(result$z_true, result$gap_cluster)
  # --------------------------
  # plot the layer detection result
  range_dij$res <- result$gap_cluster[order(mean_dij)]
  original_labels <- unique(range_dij$res)
  new_labels <- 1:length(original_labels) 
  range_dij$res <- match(range_dij$res, original_labels)
  
  ## add step func
  mean_dij_per_group <- range_dij %>%
    group_by(res) %>%
    summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
              max_re_id = max(re_id))
  
  # Create the plot with the modified color palette
  p_gap <- ggplot() +
    geom_point(data = range_dij, mapping = aes(x = re_id, y = mean_dij, color = factor(res))) +
    scale_color_manual(values = my_colors) +
    labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
    geom_segment(data = mean_dij_per_group, 
                 aes(x = min_re_id, xend = max_re_id, 
                     y = mean_dij, yend = mean_dij, color = factor(res)),
                 linewidth = 1.2) + 
    ggtitle("GAP")
  ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\gap_result\\z_gap_", rep , ".png"), p_gap,
         width=2200, height=1500, units = "px")
  
  
  # -----------------------------------------------------------------------
  # # 4. DBSCAN -- graph-based
  start.time <- Sys.time()
  input_d_matrix <- matrix(input_d, ncol = 1)
  dbscan_result <- dbscan(input_d_matrix, eps = 5, minPts = 5)
  dbscan_numlayer <- length(unique(dbscan_result$cluster[dbscan_result$cluster != 0]))
  end.time <- Sys.time()
  dbscan_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)
  
  # ARI
  no_noise_index <- which(dbscan_result$cluster != 0)
  dbscan_no_noise <- dbscan_result$cluster[no_noise_index]
  dbscan_ari <- adjustedRandIndex(dbscan_no_noise, z_true[no_noise_index])
  
  result$dbscan_cluster <- dbscan_result$cluster
  result$dbscan_ari <- dbscan_ari
  
  dbscan_nmi <- NMI(result$z_true[no_noise_index], dbscan_no_noise)
  dbscan_ami <- AMI(result$z_true[no_noise_index], dbscan_no_noise)
  dbscan_fmi <- fmi(result$z_true[no_noise_index], dbscan_no_noise)
  dbscan_cor <- cor.spearman(result$z_true[no_noise_index], dbscan_no_noise)
  # --------------------------
  # plot the layer detection result
  range_dij = cbind(min_dij, max_dij, mean_dij)
  range_dij = as.data.frame(range_dij)
  range_dij_dbscan = range_dij[no_noise_index,]
  range_dij_dbscan = range_dij_dbscan[order(mean_dij[no_noise_index]), ]
  range_dij_dbscan$re_id = rep(1:length(no_noise_index))
  
  range_dij_dbscan$res <- dbscan_no_noise[order(mean_dij[no_noise_index])]
  original_labels <- unique(range_dij_dbscan$res)
  new_labels <- 1:length(original_labels) 
  range_dij_dbscan$res <- match(range_dij_dbscan$res, original_labels)
  
  ## add step func
  mean_dij_per_group <- range_dij_dbscan %>%
    group_by(res) %>%
    summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
              max_re_id = max(re_id))
  
  # Create the plot with the modified color palette
  if(dbscan_numlayer <= 36){
    p_dbscan <- ggplot() +
      geom_point(data = range_dij_dbscan, mapping = aes(x = re_id, y = mean_dij, color = factor(res))) +
      scale_color_manual(values = my_colors) +
      labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
      geom_segment(data = mean_dij_per_group, 
                   aes(x = min_re_id, xend = max_re_id, 
                       y = mean_dij, yend = mean_dij, color = factor(res)),
                   linewidth = 1.2) + 
      ggtitle("DBSCAN")
  }else{
    my_colors_more <- rep(my_colors, ceiling(dbscan_numlayer/36))
    p_dbscan <- ggplot() +
      geom_point(data = range_dij_dbscan, mapping = aes(x = re_id, y = mean_dij, color = factor(res))) +
      scale_color_manual(values = my_colors_more) +
      labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
      geom_segment(data = mean_dij_per_group, 
                   aes(x = min_re_id, xend = max_re_id, 
                       y = mean_dij, yend = mean_dij, color = factor(res)),
                   linewidth = 1.2) + 
      ggtitle("DBSCAN")
  }
  
  ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\dbscan_result\\z_dbscan_", rep , ".png"), p_dbscan,
         width=2200, height=1500, units = "px")
  
  # -------------------------------------------------------------------------
  # 5. run MCMC_dpmm
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
  
  start.time <- Sys.time()
  result_dpmm = runMCMC_dpmm(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                             lambda=lambda, dij, tau=0.1,  mu0=0, alpha = alpha, beta = beta,
                             max_iters = Max_iteration, seed = 1)
  end.time <- Sys.time()
  dpmm_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)
  
  # ppm method
  dpmm_ppm <- gene.ppm(result_dpmm$group_iter)
  result$dpmm_ppm <- minbinder(dpmm_ppm, method = "comp")$cl
  
  dpmm_numlayer <- length(unique(result$dpmm_ppm))
  dpmm_ari <- adjustedRandIndex(result$dpmm_ppm, z_true)
  result$dpmm_ari <- dpmm_ari
  
  dpmm_nmi <- NMI(result$z_true, result$dpmm_ppm)
  dpmm_ami <- AMI(result$z_true, result$dpmm_ppm)
  dpmm_fmi <- fmi(result$z_true, result$dpmm_ppm)
  dpmm_cor <- cor.spearman(result$z_true, result$dpmm_ppm)
  
  # # --------------------------
  # # plot the layer detection result
  range_dij = cbind(min_dij, max_dij, mean_dij)
  range_dij = as.data.frame(range_dij)
  range_dij = range_dij[order(mean_dij), ]
  range_dij$re_id = rep(1:n)
  range_dij$res <- result$dpmm_ppm[order(mean_dij)]
  original_labels <- unique(range_dij$res)
  new_labels <- 1:length(original_labels)
  range_dij$res <- match(range_dij$res, original_labels)
  
  ## add step func
  mean_dij_per_group <- range_dij %>%
    group_by(res) %>%
    summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
              max_re_id = max(re_id))
  
  # Create the plot with the modified color palette
  if(dpmm_numlayer <=36){
    p_dpmm <- ggplot() +
      geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                                   xend = re_id, yend = max_dij, color = factor(res))) +
      scale_color_manual(values = my_colors) +
      labs(x = "cells", y = "distance (d_ij)", color = "layers") +
      geom_segment(data = mean_dij_per_group,
                   aes(x = min_re_id, xend = max_re_id,
                       y = mean_dij, yend = mean_dij, color = factor(res)),
                   linewidth = 1.2)+
      ggtitle("DPMM")
  }else{
    my_colors_more <- rep(my_colors, ceiling(dpmm_numlayer/36))
    p_dpmm <- ggplot() +
      geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                                   xend = re_id, yend = max_dij, color = factor(res))) +
      scale_color_manual(values = my_colors_more) +
      labs(x = "cells", y = "distance (d_ij)", color = "layers") +
      geom_segment(data = mean_dij_per_group,
                   aes(x = min_re_id, xend = max_re_id,
                       y = mean_dij, yend = mean_dij, color = factor(res)),
                   linewidth = 1.2)+
      ggtitle("DPMM")
  }
  ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\dpmm_result\\z_dpmm_", rep , ".png"), p_dpmm,
         width=2200, height=1500, units = "px")
  
  
  #--------------------------------------------------------------------------
  # p_summary <- grid.arrange(p_true, p_blade, p_kmeans, p_gmm, p_gap, p_dbscan, nrow = 2, ncol = 3)
  # ggsave(filename=paste0("Code_Xin\\Simulation\\Simulation_comparison\\true_layer\\z_summary_", rep , ".png"), p_summary,
  #        width=2200, height=1500, units = "px")
  
  save(result, file = paste0("Code_Xin\\Simulation\\Simulation_comparison\\store_result\\result_", rep , ".RData"))
  
  output_df[rep,] <- c(rep, num_layer, n, mean(sd_within), mean(sd_between), 
                       blade_numlayer, kmeans_numlayer, gmm_numlayer, gap_numlayer, dbscan_numlayer, dpmm_numlayer,
                       blade_ari, kmeans_ari, gmm_ari, gap_ari, dbscan_ari, dpmm_ari,
                       blade_nmi, kmeans_nmi, gmm_nmi, gap_nmi, dbscan_nmi, dpmm_nmi,
                       blade_ami, kmeans_ami, gmm_ami, gap_ami, dbscan_ami, dpmm_ami,
                       blade_fmi, kmeans_fmi, gmm_fmi, gap_fmi, dbscan_fmi, dpmm_fmi,
                       blade_cor, kmeans_cor, gmm_cor, gap_cor, dbscan_cor, dpmm_cor,
                       blade_time, kmeans_time, gmm_time, gap_time, dbscan_time, dpmm_time)
  # ------------------------------------------------------------------------
  cat("\n")
  cat("repetition:", rep, "\t num_layer=",num_layer, "\t num_nuclei=", n, ", sd_within=", mean(sd_within), ", sd_between=", mean(sd_between),
      ", blade_ari=", blade_ari, ", kmeans_ari=", kmeans_ari, ", gmm_ari=", gmm_ari, 
      ", gap_ari=",gap_ari, ", dbscan_ari=",dbscan_ari, "\n")
  cat("========================================================\n")
}

write.csv(output_df,
          file = paste0("Code_Xin\\Simulation\\Simulation_comparison\\4_5_Sensitive_analysis_simulation_comparison_random_6_methods.csv"), row.names = FALSE)

