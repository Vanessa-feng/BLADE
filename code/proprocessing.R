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

