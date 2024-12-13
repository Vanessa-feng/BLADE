library(Rcpp)
library(RcppArmadillo)
library(ggplot2)

proportion <- 0.25
Repetition <- 100

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
  
  save(dij, file = paste0("data\\simulated_data\\Simulated_distance\\distance_", rep , ".RData"))
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
  ggsave(filename=paste0("data\\simulated_data\\Simulated_distance\\simulated_true_layer\\z_true_", rep , ".png"), p_true,
         width=2200, height=1500, units = "px")
  
}
  