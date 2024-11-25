Max_iteration <- 2000
Burnin <- Max_iteration/2
Repetition <- 10
proportion <- 0.05 # proportion of sampling pixels in each cell
sd_between <- 50 # initial theta~N(sd_between, sd_within^2)

num_layer <- 15
num_nuclei_list <- c(20, 70, 120)
sd_within_list <- c(10, 15, 20)
sigma_list <- c(10, 15, 20)
K_list <- c(10, 20, 30)

output_df <- data.frame(matrix(nrow=length(num_nuclei_list)*length(sd_within_list)*
                                 length(sigma_list)*length(K_list)*Repetition, ncol=9))
colnames(output_df) <- c("scenario_id", "num_nuclei", "sd_within", "sigma", "K_prior", "repetition",
                         "num_layer_output", "ari", "running_time")
index <- 1
scenario_id <- 1
for(num_nuclei in num_nuclei_list){
  m <- rep(num_nuclei, num_layer)  # number of cell for each layer
  n <- sum(m)  # total number of cells
  for(sd_within in sd_within_list){
    for(sigma in sigma_list){
      for(K_prior in K_list){
        for(rep in 1:Repetition){
          if(!is.na(output_df$scenario_id[(scenario_id-1)*10+rep])){
            index <- index + 1
            next
          }
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
          
          for(k in 1:num_layer){
            z_true <- c(z_true, rep(k,m[k]))
            theta_true <- c(theta_true, rnorm(m[k], sd_between*k, sd_within))
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
          
          # # Adaptive alpha
          # alpha <- K*n   
          # 
          # # Adaptive beta
          # st2 <- var(unlist(dij))
          # beta <- n*K*st2/(16*q)
          
          ### run MCMC
          start.time <- Sys.time()
          result = runMCMC(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                           lambda=lambda, dij, G = matrix(0, n, 4), f=0,
                           tau=0.1,  mu0=0, alpha = alpha, beta = beta, K_prior=K,
                           max_iters = Max_iteration, seed = 1)
          end.time <- Sys.time()
          
          # ppm method
          ppm <- gene.ppm(result$group_iter)
          result$z_ppm <- minbinder(ppm, method = "comp")$cl
          
          numlayer_output <- length(unique(result$z_ppm))
          ari_output <- adjustedRandIndex(result$z_ppm, z_true)
          
          # running time
          time_output <- round(as.numeric(difftime(end.time, start.time, units = "secs")),3)
          output_df[index,] <- c(scenario_id, num_nuclei, sd_within, sigma, K_prior, rep, 
                                 numlayer_output, ari_output, time_output)
          cat("\n")
          cat("scenario_id:", scenario_id, "\t repetition:", rep, "\t num_nuclei=",num_nuclei, ", sd_within=", sd_within, 
              ", sigma=", sigma, ", K_prior=", K_prior, ", numlayer_ouput=", numlayer_output, 
              ", ari_ouput=",ari_output, 
              ", time_ouput=",time_output, "\n")
          index <- index + 1
        }
        scenario_id <- scenario_id + 1
        cat("========================================================\n")
        
        
      }
    }
  }
}

write.csv(output_df,
          file = paste0("C:\\oral_cancer\\Code_Xin\\Simulation\\4_2_Sensitive_analysis_simulation_numlayer_sd_sigma_K_2.csv"), row.names = FALSE)

# ===============================================================================
# draw the simulation data plot
library(ggplot2)

output_df <- read.table(paste0("D:\\UTD\\oral_cancer\\Code_Xin\\Simulation\\4_2_Sensitive_analysis_simulation_numlayer_sd_sigma_K.txt"), 
                        header=T)
# output_df$mean_val <- output_df$beta/(output_df$alpha-1)
# Plotting for ari
p <- ggplot(output_df, aes(x = as.factor(sigma), fill = as.factor(K_prior), y = ari)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_text(aes(label = as.character(round(mean_val))), 
  # vjust = -0.2, position = position_dodge(width = 1), color = "gray40", size=3) +
  facet_grid(sd_within ~ num_nuclei) +
  labs(x = "sigma", y = "K_prior", fill="K_prior",  
       title="ARI by sigma(10, 15, 20, 25) and K_prior(10, 20)", 
       subtitle="Grouped by number of nuclei per cluster(10, 100) and sd_within(15,20,25)") +
  scale_fill_manual(values = c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4","#ffffcc"))+
  theme_minimal()

ggsave(filename="D:\\UTD\\oral_cancer\\Code_Xin\\Simulation\\4_2_Sensitive_analysis_simulation_numlayer_sigma_q_K_ari.png", p)



# Plotting for time
p2 <- ggplot(output_df, aes(x = as.factor(sigma), fill = as.factor(K_prior), y = running_time)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(sd_within ~ num_nuclei) +
  labs(x = "q", y = "running time", fill="K_prior",  
       title="Running time by sigma(10, 15, 20, 25) and K_prior(10, 20)", 
       subtitle="Grouped by number of nuclei per cluster(10, 100) and sd_within(15,20,25)") +
  scale_fill_manual(values = c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4","#ffffcc"))+
  theme_minimal()

ggsave(filename="D:\\UTD\\oral_cancer\\Code_Xin\\Simulation\\4_2_Sensitive_analysis_simulation_numlayer_sigma_q_K_time.png", p2)


# Plotting for number of layer output
p3 <- ggplot(output_df, aes(x = as.factor(sigma), fill = as.factor(K_prior), y = num_layer_output)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(sd_within ~ num_nuclei) +
  labs(x = "q", y = "# of layer output", fill="K_prior",  
       title="# of layer output by sigma(10, 15, 20, 25) and K_prior(10, 20)", 
       subtitle="Grouped by number of nuclei per cluster(10, 100) and sd_within(15,20,25)") +
  scale_fill_manual(values = c("#253494", "#2c7fb8", "#41b6c4", "#a1dab4","#ffffcc"))+
  theme_minimal()

ggsave(filename="D:\\UTD\\oral_cancer\\Code_Xin\\Simulation\\4_2_Sensitive_analysis_simulation_numlayer_sigma_q_K_numlayeroutput.png", p3)


# ------------------------------------------
# matrix plot
output_df$diff <- abs(num_layer - output_df$num_layer_output)

# output_df$diff <- num_layer - output_df$num_layer_output


mp1 <- ggplot(output_df[which(output_df$num_nuclei==20),], 
              aes(as.factor(sigma), as.factor(K_prior), col = ari)) + 
  geom_point(size=10) + 
  scale_colour_gradient(low = "blue", high = "red") +
  xlab("sigma") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "number of nuclei per layer=20, ARI")

mp2 <- ggplot(output_df[which(output_df$num_nuclei==70),], 
              aes(as.factor(sigma), as.factor(K_prior), col = ari)) + 
  geom_point(size=10) + 
  scale_colour_gradient(low = "blue", high = "red") +
  xlab("sigma") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "number of nuclei per layer=100, ARI")

mp3 <- ggplot(output_df[which(output_df$num_nuclei==10),], 
              aes(as.factor(sigma), as.factor(K_prior), col = diff)) + 
  geom_point(size=10) + 
  scale_colour_gradient(low = "red", high = "blue") +
  xlab("sigma") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "number of nuclei per layer=10, diff_num_layer")

mp4 <- ggplot(output_df[which(output_df$num_nuclei==100),], 
              aes(as.factor(sigma), as.factor(K_prior), col = diff)) + 
  geom_point(size=10) + 
  scale_colour_gradient(low = "red", high = "blue") +
  xlab("sigma") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "number of nuclei per layer=100, diff_num_layer")

library(grid)
library(gridExtra)
mp_output <- grid.arrange(mp1, mp2, mp3, mp4, ncol = 2)

# -------------------------------------------------------------------

output_df_2 <- read.csv(file = paste0("C:\\oral_cancer\\Code_Xin\\Simulation\\4_2_Sensitive_analysis_simulation_numlayer_sd_sigma_K_2.csv"))
output_summary_df <- output_df_2 %>% group_by(scenario_id, num_nuclei, sd_within, sigma, K_prior) %>%
  summarise(mean_mean_num_layer_output = mean(num_layer_output),
            mean_ari = mean(ari),
            mean_running_time = mean(running_time)
  )

ari_min <- min(output_summary_df$mean_ari)
ari_max <- max(output_summary_df$mean_ari)

# Create the plots with the same color scale
mp5 <- ggplot(output_summary_df[which(output_summary_df$num_nuclei %in% c(20,70,120) & output_summary_df$sd_within==10),], 
              aes(as.factor(sigma), as.factor(K_prior), col = mean_ari)) + 
  geom_point(size=14) + 
  geom_text(aes(label = round(mean_ari, 2)), size = 4, vjust = 0.3, col="white") +
  scale_colour_gradient(low = "blue", high = "red", limits = c(ari_min, ari_max)) +
  xlab("sigma_prior") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "sigma_true=10") +
  facet_wrap(~num_nuclei, ncol=1)

mp6 <- ggplot(output_summary_df[which(output_summary_df$num_nuclei %in% c(20,70,120) & output_summary_df$sd_within==15),], 
              aes(as.factor(sigma), as.factor(K_prior), col = mean_ari)) + 
  geom_point(size=14) + 
  geom_text(aes(label = round(mean_ari, 2)), size = 4, vjust = 0.3, col="white") +
  scale_colour_gradient(low = "blue", high = "red", limits = c(ari_min, ari_max)) +
  xlab("sigma_prior") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "sigma_true=15") +
  facet_wrap(~num_nuclei, ncol=1)

mp7 <- ggplot(output_summary_df[which(output_summary_df$num_nuclei %in% c(20,70,120) & output_summary_df$sd_within==20),], 
              aes(as.factor(sigma), as.factor(K_prior), col = mean_ari)) + 
  geom_point(size=14) + 
  geom_text(aes(label = round(mean_ari, 2)), size = 4, vjust = 0.3, col="white") +
  scale_colour_gradient(low = "blue", high = "red", limits = c(ari_min, ari_max)) +
  xlab("sigma_prior") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "sigma_true=20") +
  facet_wrap(~num_nuclei, ncol=1)


title <- "ARI by prior setting K_prior(10,20,30) and sigma_prior(10,15,20)"
subtitle <- "Grouped by number of cells in each layer (20, 70, 120) and sigma_true (10, 15, 20)"
my_output <- grid.arrange(
  arrangeGrob(mp5, mp6, mp7, ncol = 3),  # Arrange the plots in a grid
  top = textGrob(title, gp = gpar(fontsize = 20, fontface = "bold")),  # Add title
  bottom = textGrob(subtitle, gp = gpar(fontsize = 15))  # Add subtitle
)

# ----------------------------------------------------------

# combine by ggplot facet_wrap.
my_output <- ggplot(output_summary_df, 
                    aes(as.factor(sigma), as.factor(K_prior), col = mean_ari)) + 
  geom_point(size=14) + 
  geom_text(aes(label = round(mean_ari, 2)), size = 4, vjust = 0.3, col="white") +
  scale_colour_gradient(low = "blue", high = "red", limits = c(ari_min, ari_max)) +
  xlab("sigma_prior") +
  ylab("K_prior") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "ARI by prior setting K_prior(10,20,30) and sigma_prior(10,15,20)",
       subtitle = "Grouped by number of cells in each layer (20, 70, 120) and sigma_true (10, 15, 20)") +
  facet_grid(num_nuclei~sd_within) 

ggsave(filename="C:\\oral_cancer\\Code_Xin\\Simulation\\4_2_Sensitive_analysis_simulation_numlayer_sigmatrue_sigmaprior_K_2.png", my_output)
