# ==========================================================================================
# patches with manual reference curves
# ==========================================================================================
# Find interest cells (R&Y cells between two manual ref curves)
# ------------------------------------------------------------------------------------------
directory_patches <- "C:\\oral_cancer\\OPMD_EPOC\\Data\\img_patches\\"
patches_files <- list.files(directory_patches, full.names = F)

# Extract file names
extract_parts <- function(file_name){
  match <- regmatches(file_name, regexpr("study-([0-9]+(-[0-9]+)?)_part_[0-9]_patch_[0-9]+", file_name))
  gsub("study-", "", match)
}

files_names <- sapply(patches_files, extract_parts, USE.NAMES = FALSE)


df_output <- data.frame(file=character(), numlayer_output=numeric(), time_output=numeric(), stringsAsFactors = FALSE)

# # update df_output
# df_output <- read.csv("C:\\oral_cancer\\OPMD_EPOC\\Data\\output.csv")

for(id in 1:length(overlap_names)){  
  cat("=======================================================================\n")
  file <- overlap_names[id]
  cat("id:", id, "\t", file, "\n")
  
  img_file_path <- paste0(directory_patches, "EPOC premalignant trial-imaging AI study-", file, ".png")
  ref_file_path <- paste0(directory_ref, "EPOC premalignant trial-imaging AI study-", file, ".png")
  
  if (!file.exists(img_file_path) || !file.exists(ref_file_path)) {
    cat("File does not exist, skipping:", file, "\n")
    next
  }
  #-----------------------
  # # check update
  if(file %in% df_output$file){
    next
  }else{
    df_output <- rbind(df_output, c(file, NA, NA, NA, NA, NA, NA, 0))
  }
  
  # # run update reference_curve, otherwise, next
  # ref_file_info <- file.info(ref_file_path)
  # if (ref_file_info$mtime < lubridate::ymd("2024-07-29")) {
  #   next
  # }
  
  # nuclei img
  img_data <- load.image(img_file_path)
  width <- dim(img_data)[1]
  height <- dim(img_data)[2]
  
  interest_pixels <- R(img_data) > 0.8 & G(img_data) < 0.5 & B(img_data) < 0.5 |
    R(img_data) > 0.8 & G(img_data) > 0.8 & B(img_data) < 0.5
  interest_img <- as.cimg(interest_pixels)
  img_labels <- label(interest_img)
  components <- as.data.frame(img_labels) %>%
    filter(value>0) %>%
    mutate(id=row_number(), cell_id=value) 
  components <- dplyr::select(components, id, cell_id, x, y)
  
  
  # filter # pixels >50
  cell_count <- components %>% group_by(cell_id) %>% summarise(pixel_count = n())
  cell_ids_large <- cell_count %>% filter(pixel_count > 30, pixel_count < 5000) %>% pull(cell_id)
  components <- components %>% filter(cell_id %in% cell_ids_large)
  if(length(unique(components$cell_id)) < 100){
    cat("Error: (1) few cells, skip. # cells = ",length(unique(components$cell_id)), "\n")
    next
  }
  
  # ----------------------------------------------------
  # ref points
  ref_img <- load.image(ref_file_path)
  interest_ref <- R(ref_img) > 0.8 & G(ref_img) < 0.5 & B(ref_img) < 0.5
  ref_bw_img <- as.cimg(interest_ref)
  
  # Label the connected components
  ref_labels <- label(ref_bw_img)
  
  # save the reference lines and regions
  ref_regions <- as.data.frame(ref_labels)
  
  if(length(unique(ref_regions$value)) < 5){
    cat("Error: check reference curves, skip. # regions = ", length(unique(ref_regions$value)), "\n")
    next
  }
  plot(ref_regions[sample(which(ref_regions$value==1),800), c("x", "y")], xlim=c(1,width), ylim=c(1,height))
  
  value_counts <- ref_regions %>%
    group_by(value) %>%
    summarise(count=n()) %>%
    arrange(desc(count))
  
  # ref curves are the least two regions after labelling
  ref_points <- ref_regions %>% filter(value==1 | value==3) %>%
    mutate(id = row_number(), boundary_id = ifelse(value == 1, 1, 2)) %>%
    sample_frac(0.25)
  ref_points <- dplyr::select(ref_points, id, boundary_id, x, y)
  ref_points$y <- height-ref_points$y
  
  # save coordination of ref lines
  write.csv(ref_points, 
            paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_ref\\", file, ".csv"), 
            row.names = F)
  
  # -----------------------------------------------------
  # create interest img mat within the two ref curves. (the closed region constains the most number of cells)
  max_nrow <- -1
  max_cell_ids_in <- NULL
  select_region <- -1
  for(region_id in 1:3){
    closed_region <- ref_regions %>% filter(value==value_counts$value[region_id])
    cell_ids_in <- components %>% inner_join(closed_region, by=c("x", "y")) %>%
      distinct(cell_id)
    current_nrow <- nrow(cell_ids_in)
    if(current_nrow > max_nrow){
      select_region <- value_counts$value[region_id]
      max_nrow <- current_nrow
      max_cell_ids_in <- cell_ids_in
    }
  }
  # closed_region <- ref_regions %>% filter(value==0)
  # cell_ids_in <- components %>% inner_join(closed_region, by=c("x", "y")) %>%
  #   distinct(cell_id)
  # components <- components %>% filter(cell_id %in% cell_ids_in$cell_id)
  
  components <- components %>% filter(cell_id %in% max_cell_ids_in$cell_id)
  if(length(unique(components$cell_id)) < 100){
    cat("Error: (2) few cells, skip. # cells = ",length(unique(components$cell_id)), "\n")
    next
  }
  components$y <- height-components$y
  # save coordination of interest pixels
  write.csv(components, 
            paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_patches\\", file, ".csv"), 
            row.names = F)
  
  p1 <- ggplot() +
    geom_point(data = ref_points, aes(x = x, y = y), color = "red", size = 0.5) +
    geom_point(data = components, aes(x = x, y = y), color = "black", size = 0.5) +
    labs(title = "Selected Pixels and Reference Lines",
         x = "X Coordinate",
         y = "Y Coordinate") +
    theme_minimal() +
    xlim(0, width) + ylim(0, height)
  ggsave(paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_patches\\img_patches_interest_",file, ".png"), 
         p1, width=width, height=height, units = "px")
  
  p2 <- ggplot() +
    geom_point(data = as.data.frame(ref_points[which(ref_points$boundary_id==1),]), aes(x = x, y = y), color = "red", size = 0.5) +
    geom_point(data = ref_points[which(ref_points$boundary_id==2),], aes(x = x, y = y), color = "green", size = 0.5) +
    geom_point(data = components, aes(x = x, y = y), color = "black", size = 0.5) +
    labs(title = "Selected Pixels and Reference Lines",
         x = "X Coordinate",
         y = "Y Coordinate") +
    theme_minimal() +
    xlim(0, width) + ylim(0, height)
  ggsave(paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_ref\\img_patches_ref_",file, ".png"),
         p2, width=width, height=height, units = "px")
  
  
  # ---------------------------------------------------------
  # Get dij
  num_ref <- length(unique(ref_points$boundary_id))
  num_cell <- length(unique(components$cell_id))
  proportion <- 0.25  # speed up
  
  dij <- list()
  for(i in 1:num_cell){
    pixel <- components[components$cell_id == unique(components$cell_id)[i], c("x", "y")]
    num_pixel <- nrow(pixel)
    pixel_sample <- pixel[sample(1:num_pixel, proportion*num_pixel),]
    dist_mat <- matrix(0, dim(pixel_sample)[1], num_ref)
    for(j in 1:num_ref){
      dist_mat[, j] = distance(as.matrix(pixel_sample), 
                               as.matrix(ref_points[ref_points$boundary_id == j, c("x", "y")]))
    }
    dij = append(dij, list(dist_mat))
  }
  
  save(dij, file = paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_distance\\distance_", file , ".RData"))
}

# ------------------------------------------
# BLADE Model
df_output <- read.csv("C:\\oral_cancer\\OPMD_EPOC\\Data\\output_sigma_15_alpha_n.csv")
for(id in 1:length(overlap_names)){
  cat("=======================================================================\n")
  file <- overlap_names[id]
  # cat("id:", id, "\t", file, "\n")
  # if(file %in% df_output$file){
  #   next
  # }
  fill_index <- 
    select_ref_id <- df_ref_output[which(df_ref_output$file==file), "select_ref_id"]
  df_output[which(df_output$file==file),"select_ref_id"] <- select_ref_id
}

df_output$select_ref_id <- as.numeric(df_output$select_ref_id)

lambda_store <- numeric(nrow(df_output))
for(id in 1:nrow(df_output)){
  # for(id in sampled_id_list){
  # if(!is.na(df_output$numlayer_output[id]) & df_output$numlayer_output[id]!=0){
  #   next
  # }
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
  
  # ------------------------------------------------------------
  
  fill_index <- which(file==df_output$file)
  
  # BLADE
  if(is.na(df_output[fill_index, "select_ref_id"])){
    cat("select_ref_id=None.\n")
    next
  }
  
  R <- df_output[fill_index, "select_ref_id"]
  dij_R <- lapply(dij, function(df) df[, R])
  
  K_prior <- 20
  set.seed(0)
  
  n <-  length(dij_R)
  dij_R_max <- sapply(dij_R, max)
  dij_R_min <- sapply(dij_R, min)
  lambda_est <- dij_R_max - dij_R_min
  lambda_store[id] <- mean(lambda_est)
  
  lambda <- lambda_est * 1.01
  
  K <- K_prior
  theta_init <- sapply(dij_R, mean)
  a1 <- (theta_init - sapply(dij_R, min)) / lambda
  b1 <- 1- (sapply(dij_R, max) - theta_init) /lambda
  r_init <- (a1 + b1)/2
  total_init <- runif(n, min = 1, max = 20)
  alpha_init <- total_init * r_init
  beta_init <- total_init * (1 - r_init)
  
  z_init = c(sample(1:K, K, replace=F), sample(1:K, n-K, replace=T))
  mu_init <- rep(0, K)
  Sigma_init <- rep(1, K)
  
  # Adaptive alpha
  alpha <- n
  
  # # Adaptive beta
  sigma <- 15
  beta <- sigma^2*(alpha-1)
  
  df_output$alpha[fill_index] <- alpha
  df_output$beta[fill_index] <- beta
  df_output$sigma[fill_index] <- sigma
  df_output$num_cell[fill_index] <- num_cell
  
  ### run MCMC
  start.time <- Sys.time()
  result <- runMCMC(z_init-1, alpha_init, beta_init, theta_init, mu_init, Sigma_init,
                    lambda=lambda, dij_R, G = matrix(0, n, 4), f=0,
                    tau=0.1,  mu0=0, alpha = alpha, beta = beta,  K_prior=K, 
                    max_iters = Max_iteration, seed = 1)
  end.time <- Sys.time()
  
  time_store <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 3)
  ppm <- gene.ppm(result$group_iter)
  z_ppm <- minbinder(ppm, method = "comp")$cl
  num_layer <- length(unique(z_ppm))
  
  df_output$running_time[fill_index] <- time_store
  result$z_ppm <- z_ppm
  
  cat("file:", file, ", R=", R, ", sigma=", sigma, ", alpha=", alpha, ", beta=", beta,
      ", numlayer_ouput=", num_layer, ", time_ouput=", time_store, "\n")
  
  
  df_output$numlayer_output[fill_index] <- num_layer
  save(result, file = paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\data_result\\sigma_15_alpha_n\\result_", file , ".RData"))
  
  
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
}

write.csv(df_output, 
          paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\output_sigma_15_alpha_n.csv"), 
          row.names = F)

