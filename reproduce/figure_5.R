setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(Rcpp)
library(RcppArmadillo)
sourceCpp("Code_bencong/blade.cpp")
sourceCpp("Code_Xin/Simulation/dpmm.cpp")

# randomly sample scenario -> simulated oral sample images
library(mclust)
library(mcclust)
library(cluster)
library(kernlab) 
library(dbscan)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(kernlab) 
library(gridExtra)
library(truncnorm)
library(reshape2)

Max_iteration <- 2000
Burnin <- Max_iteration/2
# ---------------------------------------
# # 1. BLADE
file <- "001_part_1_patch_2"

img_file_path <- paste0(directory_patches, "EPOC premalignant trial-imaging AI study-", file, ".png")
ref_file_path <- paste0(directory_ref, "EPOC premalignant trial-imaging AI study-", file, ".png")

if (!file.exists(img_file_path) || !file.exists(ref_file_path)) {
  cat("File does not exist, skipping:", file, "\n")
  next
}

# nuclei img
img_data <- load.image(img_file_path)
width <- dim(img_data)[1]
height <- dim(img_data)[2]

ref_file <- paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_ref\\", file, ".csv")
components_file <- paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_patches\\", file, ".csv")
if (!file.exists(ref_file)) {
  cat("File ref does not exist, skipping:", file, "\n")
  next
}
if (!file.exists(components_file)) {
  cat("File components does not exist, skipping:", file, "\n")
  next
}

ref_points <- read.csv(ref_file)
components <- read.csv(components_file)

load(paste0("results/epoc_pathology_image_output/result_", file, ".RData"))
     
# re-order z_ppm according to the mean of theta
theta_means <- rowMeans(result$theta_iter[, (Burnin + 1):2000])
theta_group_means <- tapply(theta_means, result$z_ppm, mean)
sorted_groups <- sort(theta_group_means, index.return = TRUE)
result$z_ppm <- match(result$z_ppm, sorted_groups$ix)

# draw the cluster result
num_cell <- nrow(result$group_iter)
for(i in 1:num_cell){
  components[components$cell_id==unique(components$cell_id)[i], "label"] <- result$z_ppm[i]
}

p1 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                 boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                 color_palette=my_colors)

# -----------------------------------------------------------------------
set.seed(0618)

num_ref <- length(unique(ref_points$boundary_id))
num_cell <- length(unique(components$cell_id))
comparison_output$num_cell[id] <- num_cell

# select reference
fill_index <- which(file==df_output$file)
if(is.na(df_output[id, "select_ref_id"])){
  cat("select_ref_id=None.\n")
  next
}

R <- df_output[fill_index, "select_ref_id"]
dij_R <- lapply(dij, function(df) df[, R])
input_d <- sapply(dij_R, mean)

# # 2. kmeans + silhouette method
# choose k
k_candidate <- 5:35
silhouette_score <- function(k) {
  km <- kmeans(input_d, centers = k, nstart = 10)
  sil <- silhouette(km$cluster, dist(input_d))
  mean(sil[, 3]) # Average silhouette width
}

# Calculate silhouette scores for each k
silhouette_scores <- sapply(k_candidate, silhouette_score)

# Find the optimal k based on the maximum silhouette score
optimal_k_silhouette <- k_candidate[which.max(silhouette_scores)]
# plot(k_candidate, silhouette_scores, type = "b", pch = 19,
#      xlab = "Number of Clusters (k)", ylab = "Average Silhouette Score",
#      main = "Silhouette Method for Optimal k")
# abline(v = optimal_k_silhouette, col = "red", lty = 2)
# text(optimal_k_silhouette, max(silhouette_scores), labels = optimal_k_silhouette, pos = 3)


# model
kmeans_output <- kmeans(input_d, centers = optimal_k_silhouette, nstart = 20)
sorted_centers <- sort(kmeans_output$centers)
cluster_mapping <- match(kmeans_output$centers, sorted_centers)
z_kmeans <- cluster_mapping[kmeans_output$cluster]

# draw the cluster result
for(i in 1:num_cell){
  components[components$cell_id==unique(components$cell_id)[i], "label"] <- z_kmeans[i]
}

p2 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                 boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                 color_palette=my_colors)

# -----------------------------------------------------------------------
# # 3. GMM+BIC
# choose k
gmm <- Mclust(input_d, G=10:35, modelNames="E", verbose=FALSE)
optimal_k_gmm <- gmm$G
comparison_output$K_GMM_BIC[id] <- optimal_k_gmm

# plot(gmm, what="BIC", xlab="Number of Clusters (k)")
# abline(v = optimal_k_gmm, col = "red", lty = 2)
# title(main = paste0("BIC Values for Different Number of Clusters (", file, ")\nOptimal k =", optimal_k_gmm)) 

# model
gmm <- Mclust(input_d, G=10, n_init=10)
sorted_centers <- sort(gmm$parameters$mean)
cluster_mapping <- match(gmm$parameters$mean, sorted_centers)
z_gmm <- as.integer(factor(cluster_mapping[gmm$classification]))

# draw the cluster result
for(i in 1:num_cell){
  components[components$cell_id==unique(components$cell_id)[i], "label"] <- z_gmm[i]
}

p3 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                 boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                 color_palette=my_colors)

# -----------------------------------------------------------------------
# # 4. Gap
# choose k
k_candicate <- 5:35
gap_stat <- clusGap(as.matrix(input_d), FUN = kmeans, nstart = 25, K.max = 35, B = 100, verbose=F)
optimal_k_gap <- which.max(gap_stat$Tab[k_candicate, "gap"])+4
comparison_output$K_Gap[id] <- optimal_k_gap

gap_values <- gap_stat$Tab[k_candicate, "gap"]
# plot(k_candicate, gap_values, type = "b", pch = 19,
#      xlab = "Number of Clusters (K)", ylab = "Gap Statistic",
#      main = paste0("Gap Statistics for K-Means Clustering (", file, ")"))
# points(optimal_k_gap, gap_stat$Tab[optimal_k_gap - 4, "gap"], col = "red", pch = 19)
# text(optimal_k_gap, gap_stat$Tab[optimal_k_gap - 4, "gap"], labels = paste("K =", optimal_k_gap), pos = 4, col = "red")

# model
kmeans_output <- kmeans(input_d, centers = optimal_k_gap, nstart = 20)
sorted_centers <- sort(kmeans_output$centers)
cluster_mapping <- match(kmeans_output$centers, sorted_centers)
z_gap <- cluster_mapping[kmeans_output$cluster]
# draw the cluster result
for(i in 1:num_cell){
  components[components$cell_id==unique(components$cell_id)[i], "label"] <- z_gap[i]
}

p4 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                 boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                 color_palette=my_colors)

# -----------------------------------------------------------------------
# # 5. DBSCAN
# choose eps
input_d_matrix <- matrix(input_d, ncol = 1)
dbscan_result <- dbscan(input_d_matrix, eps = 5, minPts = 5)
optimal_k_dbscan <- length(unique(dbscan_result$cluster[dbscan_result$cluster != 0]))
comparison_output$K_dbscan[id] <- optimal_k_dbscan

eps_list <- 5:25
optimal_k_dbscan_list <- rep(0,length(eps_list))
for(i in 1:length(eps_list)){
  dbscan_result <- dbscan(input_d_matrix, eps = eps_list[i], minPts = 5)
  optimal_k_dbscan <- length(unique(dbscan_result$cluster[dbscan_result$cluster != 0]))
  optimal_k_dbscan_list[i] <- optimal_k_dbscan
}

# plot layering assignment
cluster_no_0 <- dbscan_result$cluster[dbscan_result$cluster != 0]
cluster_means <- sapply(unique(cluster_no_0), function(i) {
  colMeans(input_d_matrix[dbscan_result$cluster == i, , drop = FALSE])
})
cluster_mapping <- match(cluster_means, sort(cluster_means))
z_dbscan <- sapply(dbscan_result$cluster, function(x){
  if(x==0){
    x
  }else{
    cluster_mapping[x]
  }
})

for (i in 1:num_cell) {
  components[components$cell_id == unique(components$cell_id)[i], "label"] <- z_dbscan[i]
}

# Draw the cluster result
p5 <-  plot.label(label=as.factor(components$label[which(components$label!=0)]), loc=components[components$label!=0,c("x","y")],
                  boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                  color_palette=my_colors)
  

