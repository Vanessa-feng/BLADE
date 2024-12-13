library(png)
library(EBImage)
library(dplyr)
library(ggplot2)
library(MASS)
library(imager)
library(lubridate)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

--------------------------------------------
  # Compare with true layers
  output_df <- read.csv(file = paste0("Code_Xin\\Simulation\\Simulation_comparison\\4_5_Sensitive_analysis_simulation_comparison_random_6_methods.csv"))
  
  # Define a function to create the plots
  create_layer_plot <- function(predicted, method_name) {
    ggplot(output_df, aes_string(x = predicted, y = "num_layer")) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      coord_fixed() +  # Ensure equal scaling of both axes
      scale_x_continuous(limits = c(5, 45)) +
      scale_y_continuous(limits = c(5, 45)) +
      labs(x = "Predicted # Layers",
           y = "Actual # Layers",
           title = method_name) +
      theme_minimal() +
      theme(aspect.ratio = 1)  # Optional: to make the plot square-shaped
  }
  
  # Create plots for each method
  p_blade <- create_layer_plot("blade_numlayer", "BLADE(MFM)")
  p_dpmm <- create_layer_plot("dpmm_numlayer", "BLADE(DPMM)")
  p_kmeans <- create_layer_plot("kmeans_numlayer", "K-means")
  p_gmm <- create_layer_plot("gmm_numlayer", "GMM")
  p_gap <- create_layer_plot("gap_numlayer", "GAP")
  p_dbscan <- create_layer_plot("dbscan_numlayer", "DBSCAN")
  
  # Display the plots
  p_comparison3 <- grid.arrange(p_blade, p_dpmm, p_kmeans,p_gmm, p_gap, p_dbscan, nrow=2)
  
  t_test_result <- t.test(output_df$blade_ari, output_df$gmm_ari, paired = F)
  print(t_test_result)
  
  # Compare BLADE(MFM) and BLADE(DPMM)
  correlation <- cor(output_df$blade_numlayer, output_df$dpmm_numlayer)
  ggplot(output_df, aes(x = blade_numlayer, y = dpmm_numlayer)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_fixed() +  # Ensure equal scaling of both axes
    scale_x_continuous(limits = c(5, 35)) +
    scale_y_continuous(limits = c(5, 35)) +
    labs(x = "BLADE (MFM)",
         y = "BLADE (DPMM)",
         title = "BLADE(MFM) versus BLADE(DPMM)") +
    theme_minimal() +
    theme(aspect.ratio = 1) + 
    annotate("text", x = 30, y = 10, label = paste0("cor = ", round(correlation, 2)), color="blue", size = 5)
  