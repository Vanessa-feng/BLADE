setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

source("code/function.R")

library(imager)
library(dplyr)
library(ggplot2)

file <- 1
Max_iteration <- 2000
Burnin <- Max_iteration/2

directory_patches <- "data/nodule_data/"
img_file_path <- paste0(directory_patches, "nodule0", file, ".png")

if (!file.exists(img_file_path)) {
  cat("File does not exist, skipping:", file, "\n")
  next
}

# nuclei img
img_data <- load.image(img_file_path)
width <- dim(img_data)[1]
height <- dim(img_data)[2]


components <- read.csv(paste0("data/nodule_data/nodule_cell_", file, ".csv"))
ref_points <- read.csv(paste0("data/nodule_data/nodule_ref_", file, ".csv"))
load(file = paste0("results/nodule_output/result_nodule_", file, ".RData"))  # result

# ----------------------------------------------------------------------------
# Use most frequency cluster of group_iter after burnin as result
# Analysis the group_iter is converge or not-----use precision

# re-order z_ppm according to the mean of theta
n <- length(result$blade_ppm)
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
cluster_plot <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")], 
                           boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                           color_palette=my_colors)

uncertainty_plot <- plot.preci(label=1-components$precision, loc=components[,c("x","y")], 
                               boundary=ref_points[ref_points$boundary_id==R,])

print(cluster_plot)
print(uncertainty_plot)
