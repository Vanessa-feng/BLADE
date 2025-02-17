setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(ggplot2)
library(gridExtra)

--------------------------------------------
# Case study: compare with true layers for 6 methods
output_df <- read.csv(file = paste0(getwd(), "/results/epoc_pathology_image_output/output_comparison.csv"))
  
# Define a function to create the scatter plots
create_layer_plot <- function(predicted, method_name) {
  correlation <- cor(output_df[[predicted]], output_df$mean_manual_count, use = "complete.obs")
  ggplot(output_df, aes_string(x = predicted, y = "mean_manual_count")) +
    geom_point(size=2) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_fixed() +  # Ensure equal scaling of both axes
    scale_x_continuous(limits = c(5, 45)) +
    scale_y_continuous(limits = c(5, 45)) +
    labs(x = paste(method_name ,"output"),
         y = "Manual counts",
         title = method_name) +
    theme_minimal() +
    theme(aspect.ratio = 1) + 
    annotate("text", x = 40, y = 7, label = paste0("Correlation = ", round(correlation, 2)), 
             color = "blue", size = 3, hjust = 1)
}

# Create plots for each method
p_blade <- create_layer_plot("blade_numlayer", "BLADE(MFM)")
p_kmeans <- create_layer_plot("kmeans_numlayer", "K-means")
p_gmm <- create_layer_plot("gmm_numlayer", "GMM")
p_gap <- create_layer_plot("gap_numlayer", "GAP")
p_dbscan <- create_layer_plot("dbscan_numlayer", "DBSCAN")
p_op <- create_layer_plot("op_numlayer", "Onion Peeling")

# Display the plots
p_comparison3 <- grid.arrange(p_blade, p_kmeans,p_gmm, p_gap, p_dbscan, p_op,  nrow=2)
print(p_comparison3)