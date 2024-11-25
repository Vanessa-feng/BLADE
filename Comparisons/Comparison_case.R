set.seed(0618)
n_patch <- nrow(df_output)
comparison_output <- data.frame(file = df_output$file,
                                select_ref_id = df_output$select_ref_id,
                                num_cell = numeric(n_ref), 
                                K_kmeans_silhouette = numeric(n_ref),
                                K_GMM_BIC = numeric(n_ref))
# id_record <- integer(0)
for(id in 1:701){
  cat("=======================================================================\n")
  file <- df_output$file[id]
  cat("id:", id, "\t", file, "\n")
  
  
  img_file_path <- paste0(directory_patches, "EPOC premalignant trial-imaging AI study-", file, ".png")
  ref_file_path <- paste0(directory_ref, "EPOC premalignant trial-imaging AI study-", file, ".png")
  
  if (!file.exists(img_file_path) || !file.exists(ref_file_path)) {
    cat("File does not exist, skipping:", file, "\n")
    next
  }
  #-----------------------
  # # check update
  # if(file %in% df_output$file){
  #   next
  # }
  
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
  
  load(file = paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_distance\\distance_", file , ".RData"))
  
  num_ref <- length(unique(ref_points$boundary_id))
  num_cell <- length(unique(components$cell_id))
  comparison_output$num_cell[id] <- num_cell
  # ---------------------------------------
  # select reference
  fill_index <- which(file==df_output$file)
  if(is.na(df_output[id, "select_ref_id"])){
    cat("select_ref_id=None.\n")
    next
  }
  R <- df_output[fill_index, "select_ref_id"]
  dij_R <- lapply(dij, function(df) df[, R])
  input_d <- sapply(dij_R, mean)
  
  
  # --------------------------------------
  # 1. k-means + silhouette method
  # k_candidate <- 5:35
  # silhouette_score <- function(k) {
  #   km <- kmeans(input_d, centers = k, nstart = 10)
  #   sil <- silhouette(km$cluster, dist(input_d))
  #   mean(sil[, 3]) # Average silhouette width
  # }
  # 
  # # Calculate silhouette scores for each k
  # silhouette_scores <- sapply(k_candidate, silhouette_score)
  # 
  # # Find the optimal k based on the maximum silhouette score
  # optimal_k_silhouette <- k_candidate[which.max(silhouette_scores)]
  # comparison_output$K_kmeans_silhouette[id] <- optimal_k_silhouette
  # 
  # png(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\kmeans_result\\silhouette_", file , ".png"), width=800, height=600)
  # plot(k_candidate, silhouette_scores, type = "b", pch = 19,
  #      xlab = "Number of Clusters (k)", ylab = "Average Silhouette Score",
  #      main = "Silhouette Method for Optimal k")
  # abline(v = optimal_k_silhouette, col = "red", lty = 2)
  # text(optimal_k_silhouette, max(silhouette_scores), labels = optimal_k_silhouette, pos = 3)
  # dev.off()
  
  # --------------------------------------
  # 2. k-means + silhouette method + constrain: diff between clusters < 150
  # k_candidate <- integer(0)
  # silhouette_scores <- integer(0)
  # for (k in 5:35) {
  #   km <- kmeans(input_d, centers = k, nstart = 20)
  #   cluster_diff <- diff(sort(km$centers[,1]))
  # 
  #   # Check if all cluster differences are within the 10-100 range
  #   if (all(cluster_diff <= 150)) {
  #     k_candidate <- c(k_candidate, k)
  #     sil <- silhouette(km$cluster, dist(input_d))
  #     silhouette_scores <- c(silhouette_scores, mean(sil[,3]))  # Average silhouette width
  #   }
  # }
  # if(length(silhouette_scores)!=0){
  #   optimal_k_silhouette <- k_candidate[which.max(silhouette_scores)]
  #   comparison_output$K_kmeans_silhouette[id] <- optimal_k_silhouette
  # }else{
  #   cat("There is no required k for kmeans method\n")
  #   id_record <- c(id_record, id)
  # }
  # 
  # silhouette_df <- data.frame(k = k_candidate, silhouette_score = silhouette_scores)
  # 
  # # Plot silhouette scores
  # ggplot(silhouette_df, aes(x = k, y = silhouette_score)) +
  #   geom_line() +
  #   geom_point() +
  #   ggtitle(paste0("Silhouette Scores for Different k for ", file)) +
  #   xlab("Number of Clusters") +
  #   ylab("Silhouette Score")
  
  
  # ----------------------------------------
  # 3. GMM + BIC method
  gmm <- Mclust(input_d, G=10:35, modelNames="E", verbose=FALSE)
  optimal_k_gmm <- gmm$G
  comparison_output$K_GMM_BIC[id] <- optimal_k_gmm
  
  png(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\gmm_result\\bic_", file , ".png"), width=800, height=600)
  plot(gmm, what="BIC", xlab="Number of Clusters (k)")
  abline(v = optimal_k_gmm, col = "red", lty = 2)
  title(main = "BIC Values for Optimal K")  # Add the title separately
  
  title(main = paste0("BIC Values for Different Number of Clusters (", file, ")\nOptimal k =", optimal_k_gmm))  # Add the title separately
  dev.off()
  
  # ----------------------------------------
  # 4. GMM + BIC + constrain: diff between clusters < 150
  # k_candidate <- integer(0)
  # BIC_record <- integer(0)
  # for(k in 5:35){
  #   gmm <- Mclust(input_d, G=k, verbose=F)
  #   gmm_diff <- diff(sort(gmm$parameters$mean))
  #   if (all(gmm_diff <= 150)) {
  #     k_candidate <- c(k_candidate, k)
  #     BIC_record <- c(BIC_record, gmm$bic)
  #   }
  # }
  # if(length(BIC_record)!=0){
  #   optimal_k_gmm <- k_candidate[which.max(BIC_record)]
  #   comparison_output$K_GMM_BIC[id] <- optimal_k_gmm
  # }else{
  #   cat("There is no required k for GMM method\n")
  #   id_record <- c(id_record, id)
  # }
  # 
  
  # # Plot GMM+BIC values
  # gmm_df <- data.frame(k = k_candidate, BIC_record = BIC_record)
  # 
  # # Plot gmm scores
  # ggplot(gmm_df, aes(x = k, y = BIC_record)) +
  #   geom_line() +
  #   geom_point() +
  #   ggtitle(paste0("BIC by GMM for Different k for ", file)) +
  #   xlab("Number of Clusters") +
  #   ylab("BIC")
  
  # ----------------------------------------
  # # 5. Gap statistic -- distance-based
  k_candicate <- 5:35
  gap_stat <- clusGap(as.matrix(input_d), FUN = kmeans, nstart = 25, K.max = 35, B = 100, verbose=F)
  optimal_k_gap <- which.max(gap_stat$Tab[k_candicate, "gap"])+4
  comparison_output$K_Gap[id] <- optimal_k_gap
  
  gap_values <- gap_stat$Tab[k_candicate, "gap"]
  png(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\gap_result\\gap_", file , ".png"), width=800, height=600)
  plot(k_candicate, gap_values, type = "b", pch = 19,
       xlab = "Number of Clusters (K)", ylab = "Gap Statistic",
       main = paste0("Gap Statistics for K-Means Clustering (", file, ")"))
  points(optimal_k_gap, gap_stat$Tab[optimal_k_gap - 4, "gap"], col = "red", pch = 19)
  text(optimal_k_gap, gap_stat$Tab[optimal_k_gap - 4, "gap"], labels = paste("K =", optimal_k_gap), pos = 4, col = "red")
  dev.off()
  
  plot(k_candicate, gap_values, type = "b", pch = 19,
       xlab = "Number of Clusters (K)", ylab = "Gap Statistic",
       main = "Gap Statistics for K-Means Clustering")
  points(optimal_k_gap, gap_stat$Tab[optimal_k_gap - 4, "gap"], col = "red", pch = 19)
  abline(v = optimal_k_gap, col = "red", lty = 2)
  # 
  # ----------------------------------------
  # # 6. DBSCAN -- graph-based
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
  plot(eps_list, optimal_k_dbscan_list, type="b", pch = 19, xlab="eps", ylab="Optimal k", 
       main="Optimal K with varying eps")
  # ----------------------------------------
  # 7. FANNY (Fuzzy clustering) -- distance-based
  ## Fanny aims to minimize the objective function
  ## Slow
  
  # fanny_clusters <- function(data, k) {
  #   fanny_result <- fanny(data, k = k)
  #   return(fanny_result$objective)
  # }
  # 
  # k_range <- 5:35
  # fanny_results <- sapply(k_range, function(k) fanny(input_d, k)$objective[1])
  # plot(k_range, fanny_results, type = "b", xlab = "Number of Clusters (K)", ylab = "Objective Function Value", main = "Fanny Method")
  # k_optimal <- k_range[which.min(fanny_results)]
  
  
  # ----------------------------------------
  
  
  
  # -----------------------------------------------------------------------------------
  # plot labelled nuclei
  # ---------------------------------------
  # # 1. kmeans
  kmeans_output <- kmeans(input_d, centers = optimal_k_silhouette, nstart = 20)
  sorted_centers <- sort(kmeans_output$centers)
  cluster_mapping <- match(kmeans_output$centers, sorted_centers)
  z_kmeans <- cluster_mapping[kmeans_output$cluster]
  # draw the cluster result
  for(i in 1:num_cell){
    components[components$cell_id==unique(components$cell_id)[i], "label"] <- z_kmeans[i]
  }
  
  p1 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                   boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                   color_palette=my_colors)
  ggsave(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\kmeans_result\\label_", file , ".png"), p1,
         width=width, height=height, units = "px")
  # 
  # 
  # # 3. GMM+BIC
  gmm <- Mclust(input_d, G=10, n_init=10)
  sorted_centers <- sort(gmm$parameters$mean)
  cluster_mapping <- match(gmm$parameters$mean, sorted_centers)
  z_gmm <- as.integer(factor(cluster_mapping[gmm$classification]))
  
  # draw the cluster result
  for(i in 1:num_cell){
    components[components$cell_id==unique(components$cell_id)[i], "label"] <- z_gmm[i]
  }
  
  p2 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                   boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                   color_palette=my_colors)
  ggsave(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\gmm_result\\label_", file , ".png"), p2,
         width=width, height=height, units = "px")
  
  
  # ------------------------------------------------
  # # 5. Gap
  kmeans_output <- kmeans(input_d, centers = optimal_k_gap, nstart = 20)
  sorted_centers <- sort(kmeans_output$centers)
  cluster_mapping <- match(kmeans_output$centers, sorted_centers)
  z_gap <- cluster_mapping[kmeans_output$cluster]
  # draw the cluster result
  for(i in 1:num_cell){
    components[components$cell_id==unique(components$cell_id)[i], "label"] <- z_gap[i]
  }
  
  p5 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                   boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                   color_palette=my_colors)
  ggsave(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\gap_result\\label_", file , ".png"), p5,
         width=width, height=height, units = "px")
  
  
  
  # ------------------------------------------------
  # # 6. DBSCAN
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
  p6 <-  plot.label(label=as.factor(components$label[which(components$label!=0)]), loc=components[components$label!=0,c("x","y")],
                    boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                    color_palette=my_colors)
  # Save the plot
  ggsave(filename = paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\dbscan_result\\label_", file, ".png"),
         plot = p6, width = width, height = height, units = "px")
  
}
write.csv(comparison_output, 
          paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\output_comparison.csv"), 
          row.names = F)
