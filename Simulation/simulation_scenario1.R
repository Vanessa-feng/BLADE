Max_iteration <- 2000
Burnin <- Max_iteration/2
Repetition <- 10
num_layer_list <- c(10, 20, 30)
m_list <- c(20, 70, 120)
prop_list <- c(0.05, 0.1, 0.25, 0.5, 1)
sd_between <- 50
sd_within_list <- c(10, 15, 20)

# alpha <- 300
# beta <- 50000

output_df <- data.frame(matrix(nrow=length(num_layer_list)*length(m_list)*
                                 length(prop_list)*length(sd_within_list)*Repetition, ncol=9))
colnames(output_df) <- c("scenario_id", "num_layer","num_cells_each_layer", "proportion", "sd_within", "repetition",
                         "num_layer_output", "ari", "running_time")
index <- 1
scenario_id <- 1
for(num_layer in num_layer_list){
  for(mi in m_list){
    m <- rep(mi, num_layer)  # number of cell for each layer 
    n <- sum(m)  # total number of cells
    
    for(prop in prop_list){
      for(sd_within in sd_within_list){
        for(rep in 1:Repetition){
          ### Generate data
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
          
          for(k in 1:num_layer){
            z_true <- c(z_true, rep(k,m[k]))
            theta_true <- c(theta_true, rnorm(m[k], sd_between*k, sd_within))
          }
          
          lambda_true = runif(n, 20, 40)
          alpha_true = rgamma(n, 1, 0.5)+ 0.4
          beta_true = rgamma(n, 1, 0.5) + 0.4
          
          # generate cells
          for (i in 1:n){
            cell = rbeta(ni[i]*prop, alpha_true[i], beta_true[i])*lambda_true[i] +
              theta_true[i]-alpha_true[i]/(alpha_true[i]+beta_true[i])*lambda_true[i]
            dij[[i]] = cell
          }
          #### init
          K = 20
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
          sigma <- 15
          beta <- sigma^2*(alpha-1)
          
          ### run MCMC
          start.time <- Sys.time()
          result = runMCMC(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                           lambda=lambda, dij, G = matrix(0, n, 4), f=0,
                           tau=0.1,  mu0=0, alpha = alpha, beta = beta,  K_prior=K,
                           max_iters = Max_iteration, seed = 1)
          end.time <- Sys.time()
          
          # ppm method
          ppm <- gene.ppm(result$group_iter)
          result$z_ppm <- minbinder(ppm, method = "comp")$cl
          save(result, file = paste0("C:\\oral_cancer\\Code_Xin\\Simulation\\Simulation1\\result_layer_", 
                                     num_layer, "_mi_", mi, "_prop_", prop, "_sdw_", sd_within, "_rep_", rep , ".RData"))
          # check
          # result$group_iter[,Max_iteration]
          # adjustedRandIndex(result$group_iter[,Max_iteration], z_true)
          
          # # use ARI: adjust rand index
          # numlayer_store_rep[rep] <- as.numeric(names(which.max(table(result$K_iter[Burnin:Max_iteration, 1]))))
          numlayer_output <- length(unique(result$z_ppm))
          
          # find most frequency class in Burnin to Max_iteration
          # burnin_iterations <- result$group_iter[, Burnin:Max_iteration]
          # class_freq <- apply(burnin_iterations, 1, table)
          # most_freq_class <- sapply(class_freq, function(x) as.numeric(names(which.max(x))))
          
          ari_output <- adjustedRandIndex(result$z_ppm, z_true)
          
          # running time
          time_output <- round(as.numeric(difftime(end.time, start.time, units = "secs")),3)
          
          
          output_df[index,] <- c(scenario_id, num_layer, mi, prop, sd_within, rep, 
                                 numlayer_output, ari_output, time_output)
          
          cat("\n")
          cat("scenario_id:", scenario_id, "\t repetition:", rep, "\t num_layer=",num_layer, ", mi=", mi, ", proportion=", prop, 
              ", sd_within=", sd_within, ", numlayer_output=", numlayer_output, 
              ", ari_ouput=", ari_output, ", time_ouput=", time_output, "\n")
          index <- index + 1
        }
        scenario_id <- scenario_id + 1
        cat("========================================================\n")
      }
    }
  }
}

write.csv(output_df,
          file = paste0("C:\\oral_cancer\\Code_Xin\\Simulation\\4_1_Sensitive_analysis_simulation_samplesize_numlayer_proportion_sigma_2.csv"), 
          row.names = FALSE)

# draw the simulation plots
library(ggplot2)

# Plotting of ari with error bar
summary_df_ari <- output_df %>%
  group_by(sd_within, proportion, num_cells_each_layer, num_layer) %>%
  summarise(
    ari_mean = mean(ari),
    ari_min = min(ari), 
    ari_max = max(ari), 
    .groups = 'drop'
  )
p <- ggplot(summary_df_ari, aes(x = as.factor(sd_within), fill = as.factor(proportion), y = ari_mean)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ari_min, ymax = ari_max), 
                position = position_dodge(width = 0.9), width = 0.25) +  
  facet_grid(num_cells_each_layer ~ num_layer) +
  labs(x = "sd_within", y = "ARI", fill = "proportion",  
       title = "ARI by proportion(0.05,0.1,0.25,0.5,1) and sd_within(10,15,20)", 
       subtitle = "Grouped by number of layers(10,20,30) and number of cells in each layer(20,70,120)") +
  scale_fill_manual(values = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) +
  theme_minimal()

ggsave(filename="C:\\oral_cancer\\Code_Xin\\Simulation\\4_1_Sensitive_analysis_simulation_samplesize_numlayer_proportion_sigma_ari_2_col_large.png", p)

# Plotting of time
summary_df_time <- output_df %>%
  group_by(sd_within, proportion, num_cells_each_layer, num_layer) %>%
  summarise(
    running_time_mean = mean(running_time),
    running_time_min = min(running_time), 
    running_time_max = max(running_time), 
    .groups = 'drop'
  )
p2 <- ggplot(summary_df_time, aes(x = as.factor(sd_within), fill = as.factor(proportion), y = running_time_mean)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = running_time_min, ymax = running_time_max), 
                position = position_dodge(width = 0.9), width = 0.25) + 
  facet_grid(num_cells_each_layer ~ num_layer) +
  labs(x = "sd_within", y = "running_time", fill = "proportion",  
       title = "Running time by proportion(0.05,0.1,0.25,0.5,1) and sd_within(10,15,20)", 
       subtitle = "Grouped by number of layers(10,20,30) and number of cells in each layer(20,70,120)") +
  scale_fill_manual(values = c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4", "#ffffcc")) +
  theme_minimal()

ggsave(filename="C:\\oral_cancer\\Code_Xin\\Simulation\\4_1_Sensitive_analysis_simulation_samplesize_numlayer_proportion_sigma_time_2.png", p2)

# plotting of number of layer output
summary_df_numlayer <- output_df %>%
  group_by(sd_within, proportion, num_cells_each_layer, num_layer) %>%
  summarise(
    num_layer_output_mean = mean(num_layer_output),
    num_layer_output_min = min(num_layer_output), 
    num_layer_output_max = max(num_layer_output), 
    .groups = 'drop'
  )
p3 <- ggplot(summary_df_numlayer, aes(x = as.factor(sd_within), fill = as.factor(proportion), y = num_layer_output_mean)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = num_layer_output_min, ymax = num_layer_output_max), 
                position = position_dodge(width = 0.9), width = 0.25) + 
  facet_grid(num_cells_each_layer ~ num_layer) +
  labs(x = "sd_within", y = "num_layer_output", fill = "proportion",  
       title = "# of layer output by proportion(0.05,0.1,0.25,0.5,1) and sd_within(10,15,20)", 
       subtitle = "Grouped by number of layers(10,20,30) and number of cells in each layer(20,70,120)") +
  scale_fill_manual(values = c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4", "#ffffcc")) +
  theme_minimal()
ggsave(filename="C:\\oral_cancer\\Code_Xin\\Simulation\\4_1_Sensitive_analysis_simulation_samplesize_numlayer_proportion_sigma_numlayeroutput_2.png", p3)
