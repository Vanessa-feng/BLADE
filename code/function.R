library(RColorBrewer)
library(ggplot2)


set.seed(245)
my_colors <- sample(c(brewer.pal(8,"Dark2"), brewer.pal(9,"Set1"), brewer.pal(12,"Set3"), brewer.pal(12,"Paired")), 36, replace = FALSE)

# ------------------------------------------------------------------------------
# Functions for BLADE
# PPM method to obtain the assignment
gene.ppm <- function(z_store, iter=2000, burn=1001){
  n <- nrow(z_store)
  ppm <- matrix(0, nrow = n, ncol = n)
  for (ii in burn:iter){
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (z_store[i, ii] == z_store[j, ii]) {
          ppm[i, j] <- ppm[i, j] + 1
          ppm[j, i] <- ppm[j, i] + 1
        }
      }
    } 
  } 
  ppm <- ppm/(iter - burn + 1)
  diag(ppm) <- rep(1, n)
  return(ppm)
}

# Draw plot for each pixel
plot.label <- function(label, loc, boundary, main = "", width, height, color_palette) {
  data <- data.frame(expr = label, x = loc[, "x"], y = loc[, "y"]);
  data1 <- data.frame( bx = boundary[, "x"], by = boundary[, "y"])
  ggplot() + geom_point(data=data, mapping = aes(x = x, y = y, color = expr), size = 1) + 
    scale_color_manual(values = color_palette) +
    geom_point(data = data1, mapping = aes(x = bx, y = by), color = "black", size = 1) + 
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
          legend.position="bottom")+
    xlim(0, width) + ylim(0, height)
}
# ------------------------------------------------------------------------------
# Simulation: generate samples
# generate sd within
generate_sd_within <- function(x, num_nuclei) {
  # more nuclei matches larger sd_within
  norm_value <- (x - min(num_nuclei)) / (max(num_nuclei) - min(num_nuclei))
  sd_within <- rtruncnorm(1, a = 8, b = 22, mean = 15 + 7 * norm_value, sd = 2)
  return(sd_within)
}

# ------------------------------------------------------------------------------
# Comparisons: 
# k-means + silhouette method
silhouette_score <- function(k) {
  km <- kmeans(input_d, centers = k, nstart = 5)
  sil <- silhouette(km$cluster, dist(input_d))
  mean(sil[, 3]) # Average silhouette width
}

# ------------------------------------------------------------------------------
# Evaluations
# cor spearman
cor.spearman <- function(z_true, z_pred){
  conf_matrix <- table(z_true, z_pred)
  
  # Pad the confusion matrix to make it square
  n <- max(nrow(conf_matrix), ncol(conf_matrix))
  padded_matrix <- matrix(0, n, n)
  padded_matrix[1:nrow(conf_matrix), 1:ncol(conf_matrix)] <- conf_matrix
  
  assignment <- solve_LSAP(padded_matrix, maximum = TRUE)
  spearman_cor <- cor(assignment[z_true], z_pred)
  return(spearman_cor)
}

# Evaluation: Fowlkes-Mallows Index
fmi <- function(z_true, z_pred) {
  n <- length(z_true)
  TP <- sum(outer(z_true, z_true, "==") & outer(z_pred, z_pred, "==")) / 2
  FP <- sum(outer(z_true, z_true, "!=") & outer(z_pred, z_pred, "==")) / 2
  FN <- sum(outer(z_true, z_true, "==") & outer(z_pred, z_pred, "!=")) / 2
  fmi_value <- TP / sqrt((TP + FP) * (TP + FN))
  return(fmi_value)
}