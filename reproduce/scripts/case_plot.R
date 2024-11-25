
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
fill_index <- which(file==df_output$file)

# BLADE
if(is.na(df_output[fill_index, "select_ref_id"])){
  cat("select_ref_id=None.\n")
  next
}

R <- df_output[fill_index, "select_ref_id"]
dij_R <- lapply(dij, function(df) df[, R])


load(file = paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\sigma_15_alpha_n\\result_", file , ".RData"))


hist(result$K_iter[,1], nclass=50, xlab="Number of Clusters (k)", main="Histogram of K")
# ---------------------------------------
# re-order z_ppm according to the mean of theta
theta_means <- rowMeans(result$theta_iter[, (Burnin + 1):2000])
theta_group_means <- tapply(theta_means, result$z_ppm, mean)
sorted_groups <- sort(theta_group_means, index.return = TRUE)
result$z_ppm <- match(result$z_ppm, sorted_groups$ix)

# draw the cluster result
for(i in 1:num_cell){
  components[components$cell_id==unique(components$cell_id)[i], "label"] <- result$z_ppm[i]
}

p2 <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                 boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                 color_palette=my_colors)
ggsave(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\sigma_15_alpha_n\\label_", file , ".png"), p2,
       width=width, height=height, units = "px")

# # ---------------------------------------
# # draw vertical plot with label
load(paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\sigma_15_alpha_n\\result_", file , ".RData"))
n <- length(dij_R)
num_layer <- length(unique(result$z_ppm))
min_dij = sapply(dij_R, min)
max_dij = sapply(dij_R, max)
mean_dij = sapply(dij_R, mean)
range_dij = cbind(min_dij, max_dij, mean_dij)
range_dij = as.data.frame(range_dij)
range_dij = range_dij[order(mean_dij), ]
range_dij$re_id = rep(1:n)

range_dij$res <- result$z_ppm[order(mean_dij)]
original_labels <- unique(range_dij$res)
new_labels <- 1:length(original_labels) 
range_dij$res <- match(range_dij$res, original_labels)

## add step func
mean_dij_per_group <- range_dij %>%
  group_by(res) %>%
  summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
            max_re_id = max(re_id))

# Create the plot with the modified color palette
p3 <- ggplot() +
  geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                               xend = re_id, yend = max_dij, color = factor(res))) +
  scale_color_manual(values = my_colors) +
  labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
  geom_segment(data = mean_dij_per_group, 
               aes(x = min_re_id, xend = max_re_id, 
                   y = mean_dij, yend = mean_dij, color = factor(res)),
               linewidth = 1.5)
# without labelled
ggplot() +
  geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                               xend = re_id, yend = max_dij), color="gray50") +
  labs(x = "cells", y = "distance (d_ij)") 

ggsave(filename=paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\sigma_15_alpha_n\\vp_", file , ".png"), p3,
       width=2200, height=1500, units = "px")