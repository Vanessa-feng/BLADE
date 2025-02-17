setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(png)
library(EBImage)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

file <- 1
Max_iteration <- 2000
Burnin <- Max_iteration/2

# Draw plot for each pixel
plot.label <- function(label, loc, boundary, main = "", color_palette) {
  data <- data.frame(expr = label, x = loc[, 1], y = loc[, 2]);
  data1 <- data.frame( bx = boundary[,1], by = boundary[, 2])
  ggplot() + geom_point(data=data, mapping = aes(x = x, y = y, color = expr), size = 1) + 
    scale_color_manual(values = color_palette) +
    geom_point(data = data1, mapping = aes(x = bx, y = by), color = "black", size = 2) + 
    coord_fixed(ratio = 1) +  
    theme_classic() + labs(color = "", title = main) + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust = 0.5, face = "bold", size = 20),
          panel.border = element_blank(),
          legend.text = element_text(face = "bold", size = 16),
          legend.position="bottom")
}

# plot precision
plot.preci <- function(label, loc, boundary, main = "") {
  
  data <- data.frame(expr = label, x = loc[, 1], y = loc[, 2])
  data1 <- data.frame(bx = boundary[, 1], by = boundary[, 2])
  
  ggplot() +
    geom_point(data = data, mapping = aes(x = x, y = y, color = expr), size = 1) +
    geom_point(data = data1, mapping = aes(x = bx, y = by), color = "black", size = 1) +
    scale_color_gradient(
      low = "black", high = "LightSkyBlue", limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2), 
      labels = seq(0, 1, by = 0.2)  
    ) +
    coord_fixed(ratio = 1) +
    theme_classic() +
    labs(color = "", title = main) +
    guides(color = guide_colorbar(barwidth = 1, barheight = 8)) +  # Adjust bar dimensions
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      panel.border = element_blank(),
      legend.text = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
}

components <- read.csv(paste0("../data/nodule_data/nodule_cell_", file, ".csv"))
ref_points <- read.csv(paste0("../data/nodule_data/nodule_ref_", file, ".csv"))
load(file = paste0("../results/nodule_output/result_nodule_", file, ".RData"))  # result

# ----------------------------------------------------------------------------
# Use most frequency cluster of group_iter after burnin as result
# Analysis the group_iter is converge or not-----use precision

# re-order z_ppm according to the mean of theta
theta_means <- rowMeans(result$theta_iter[, (Burnin + 1):2000])
theta_group_means <- tapply(theta_means, result$blade_ppm, mean)
sorted_groups <- sort(theta_group_means, index.return = TRUE)
result$blade_ppm <- match(result$blade_ppm, sorted_groups$ix)

cluster_result <- apply(result$group_iter, 1, function(x) {
  clusters <- table(x[Burnin:Max_iteration])
  max_cluster <- names(clusters)[which.max(clusters)]
  precision <- clusters[max_cluster] / (Max_iteration-Burnin+1)
  return(list(max_cluster = max_cluster, precision = precision))
})

cluster_df <- data.frame(Cell = seq_along(cluster_result), 
                         Max_Cluster = sapply(cluster_result, function(x) x$max_cluster),
                         Precision = sapply(cluster_result, function(x) x$precision))

components$cell_id <- as.numeric(as.factor(components$cell_id))
for(i in 1:n){
  components[which(components$cell_id==i), "label"] <- result$blade_ppm[i] # or cluster_df$Max_Cluster[i]
  components[which(components$cell_id==i), "precision"] <- cluster_df$Precision[i]
}

# draw uncertainty plot and clustering results
uncertainty_plot <- plot.preci(label=1-components$precision, loc=components[,c("x","y")], 
                         boundary=ref_points[ref_points$boundary_id==R,])

cluster_plot <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")], 
                           boundary=ref_points[ref_points$boundary_id==R,],
                           color_palette=my_colors)
print(cluster_plot)
print(uncertainty_plot)