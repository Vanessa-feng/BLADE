setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

# randomly sample scenario -> simulated oral sample images

library(ggplot2)
library(gridExtra)

output_df <- read.csv(file = "results/simulation_output/simulation_comparison.csv")

# Compare with true layers
# Define a function to create the plots
create_layer_plot <- function(predicted, method_name) {
  correlation <- cor(output_df[[predicted]], output_df$num_layer, use = "complete.obs")
  ggplot(output_df, aes_string(x = predicted, y = "num_layer")) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_fixed() +  # Ensure equal scaling of both axes
    scale_x_continuous(limits = c(5, 45)) +
    scale_y_continuous(limits = c(5, 45)) +
    labs(x = "Est. of # of layer",
         y = "Manual counts",
         title = method_name) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12)
    ) +
    annotate("text", x = 30, y = 40, label = paste0("\u03C1 = ", round(correlation, 2)),
             color = "blue", size = 4, hjust = 1)
}

# Create plots for each method
p_blade <- create_layer_plot("blade_numlayer", "BLADE")
p_dpmm <- create_layer_plot("dpmm_numlayer", "BLADE (DPMM)")
p_kmeans <- create_layer_plot("kmeans_numlayer", "k-means")
p_gmm <- create_layer_plot("gmm_numlayer", "GMM")
p_gap <- create_layer_plot("gap_numlayer", "GAP")
p_dbscan <- create_layer_plot("dbscan_numlayer", "DBSCAN")

# Display the plots
p_comparison3 <- grid.arrange(p_blade, p_dpmm, p_dbscan, p_gap, p_gmm, p_kmeans,  nrow=2)
print(p_comparison3)
# ggsave(filename="results/simulation_output/simu_comparison3.png", p_comparison3, width=2400, height=1800, units = "px")
