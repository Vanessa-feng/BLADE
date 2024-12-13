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

write.csv(df_output, 
          paste0("C:\\oral_cancer\\OPMD_EPOC\\Data\\output_sigma_15_alpha_n.csv"), 
          row.names = F)

