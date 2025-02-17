setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

source("code/function.R")

# components and reference curve
file <- "001_part_1_patch_2"
directory_patches <- paste0(getwd(), "/data/epoc_pathology_image_data/slide_info/")
directory_ref <- paste0(getwd(), "/data/epoc_pathology_image_data/ref_info/")
directory_distance <- paste0(getwd(), "/data/epoc_pathology_image_data/epoc_distance/")

img_file_path <- paste0(directory_patches, "EPOC premalignant trial-imaging AI study-", file, ".png")
ref_file_path <- paste0(directory_ref, "EPOC premalignant trial-imaging AI study-", file, ".png")

if (!file.exists(img_file_path) || !file.exists(ref_file_path)) {
  stop(paste("Image patch or reference curve does not exist, skipping:", file))
}

# nuclei img
img_data <- load.image(img_file_path)
width <- dim(img_data)[1]
height <- dim(img_data)[2]

components_file <- paste0(directory_patches, file, ".csv")
ref_file <- paste0(directory_ref, file, ".csv")
if (!file.exists(ref_file)) {
  stop(paste0("File ref does not exist, skipping:", file))
}
if (!file.exists(components_file)) {
  stop(paste0("File components does not exist, skipping:", file))
}

ref_points <- read.csv(ref_file)
components <- read.csv(components_file)
load(file = paste0(directory_distance, "distance_", file , ".RData")) # dij
load(paste0(getwd(), "/results/epoc_pathology_image_output/result_", file , ".RData")) # result

# load select_ref_id
df_output <- read.csv(paste0(getwd(), "/results/epoc_pathology_image_output/epoc_output.csv"))
fill_index <- which(df_output$file == file)
if(is.na(df_output[fill_index, "select_ref_id"])){
  cat("select_ref_id=None.\n")
  next
}
R <- df_output[fill_index, "select_ref_id"]
dij_R <- lapply(dij, function(df) df[, R])

# ------------------------------------------------------------------------------
# re-order z_ppm according to the mean of theta
Max_iteration <- 2000
Burnin <- Max_iteration/2
num_cell <- length(unique(components$cell_id))
theta_means <- rowMeans(result$theta_iter[, (Burnin + 1):2000])
theta_group_means <- tapply(theta_means, result$z_ppm, mean)
sorted_groups <- sort(theta_group_means, index.return = TRUE)
result$z_ppm <- match(result$z_ppm, sorted_groups$ix)

# draw the cluster result
for(i in 1:num_cell){
  components[components$cell_id==unique(components$cell_id)[i], "label"] <- result$z_ppm[i]
}

p_blade <- plot.label(label=as.factor(components$label), loc=components[,c("x","y")],
                      boundary=ref_points[ref_points$boundary_id==R,], width=width, height=height,
                      color_palette=my_colors)
ggsave(filename=paste0(paste0(getwd(), "/results/epoc_pathology_image_output/blade_", file , "png")), p_blade,
       width=2200, height=1500, units = "px")

# ------------------------------------------------------------------------------
# draw vertical plot with label
min_dij = sapply(dij_R, min)
max_dij = sapply(dij_R, max)
mean_dij = sapply(dij_R, mean)
range_dij = cbind(min_dij, max_dij, mean_dij)
range_dij = as.data.frame(range_dij)
range_dij = range_dij[order(mean_dij), ]
range_dij$re_id = rep(1:num_cell)

range_dij$res <- result$z_ppm[order(mean_dij)]
original_labels <- unique(range_dij$res)
new_labels <- 1:length(original_labels) 
range_dij$res <- match(range_dij$res, original_labels)

# add step func
mean_dij_per_group <- range_dij %>%
  group_by(res) %>%
  summarise(mean_dij = mean(mean_dij), min_re_id = min(re_id),
            max_re_id = max(re_id))

# input data
p_input <- ggplot() +
  geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                               xend = re_id, yend = max_dij), color="gray50") +
  labs(x = "cells", y = "distance (d_ij)") 

# Create the color-coded layer detection plot
p_colorcoded <- ggplot() +
  geom_segment(data = range_dij, mapping = aes(x = re_id, y = min_dij,
                                               xend = re_id, yend = max_dij, color = factor(res))) +
  scale_color_manual(values = my_colors) +
  labs(x = "cells", y = "distance (d_ij)", color = "layers") + 
  geom_segment(data = mean_dij_per_group, 
               aes(x = min_re_id, xend = max_re_id, 
                   y = mean_dij, yend = mean_dij, color = factor(res)),
               linewidth = 1.5)

ggsave(filename=paste0(paste0(getwd(), "/results/epoc_pathology_image_output/colorcoded_", file , "png")), p_colorcoded,
       width=2200, height=1500, units = "px")

