setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(dplyr)
library(reshape2)
library(ggpubr)
library(gridExtra)

# Load the results of 200 simulations for 6 different methods
output_df <- read.csv(file = paste0(getwd(), "/results/simulation_output/simulation_comparison.csv"))

# ----------------------------------------------------
# Add significance lines and labels to the plot
my_comparisons <- list( c("BLADE", "BLADE (DPMM)"), 
                        c("BLADE", "DBSCAN"), c("BLADE", "GAP"),
                        c("BLADE", "GMM"), c("BLADE", "k-means"))

# Compare ARI
ari_df <- melt(output_df, 
               measure.vars = c("blade_ari", "dpmm_ari", "kmeans_ari", "gmm_ari", "gap_ari", "dbscan_ari"),
               variable.name = "Method", 
               value.name = "ARI")
ari_df$Method <- recode(ari_df$Method, 
                        "blade_ari" = "BLADE", 
                        "dpmm_ari" = "BLADE (DPMM)",
                        "kmeans_ari" = "k-means", 
                        "gmm_ari" = "GMM", 
                        "gap_ari" = "GAP", 
                        "dbscan_ari" = "DBSCAN")

# ----------------------------------------------------
# Difference of output (the number of layers)
output_df_long <- melt(output_df, 
                       measure.vars = c("blade_numlayer", "dpmm_numlayer", "kmeans_numlayer", 
                                        "gmm_numlayer", "gap_numlayer", "dbscan_numlayer"),
                       variable.name = "Method", 
                       value.name = "Predicted")
# output_df_long$AbsDiff <- abs(output_df_long$Predicted - output_df_long$num_layer)
output_df_long$Diff <- output_df_long$Predicted - output_df_long$num_layer
output_df_long$Method <- recode(output_df_long$Method, 
                                "blade_numlayer" = "BLADE", 
                                "dpmm_numlayer" = "BLADE (DPMM)",
                                "kmeans_numlayer" = "k-means", 
                                "gmm_numlayer" = "GMM", 
                                "gap_numlayer" = "GAP", 
                                "dbscan_numlayer" = "DBSCAN")

# t.test(output_df_long[output_df_long$Method=="BLADE", "Diff"], 
#        output_df_long[output_df_long$Method=="BLADE (DPMM)", "Diff"], paired=T, alternative="less")

# ---------------------------------------------------------
cor_df <- melt(output_df, 
               measure.vars = c("blade_cor", "dpmm_cor", "kmeans_cor", "gmm_cor", "gap_cor", "dbscan_cor"),
               variable.name = "Method", value.name = "COR")
cor_df$Method <- recode(cor_df$Method, 
                        "blade_cor" = "BLADE", 
                        "dpmm_cor" = "BLADE (DPMM)",
                        "kmeans_cor" = "k-means", 
                        "gmm_cor" = "GMM", 
                        "gap_cor" = "GAP", 
                        "dbscan_cor" = "DBSCAN")

# Change the order of methods (according to initial of methods)
ari_df$Method <- factor(ari_df$Method, levels = c("BLADE", "BLADE (DPMM)", "DBSCAN", "GAP", "GMM", "k-means"))
output_df_long$Method <- factor(output_df_long$Method, levels = c("BLADE", "BLADE (DPMM)", "DBSCAN", "GAP", "GMM", "k-means"))
cor_df$Method <- factor(cor_df$Method, levels = c("BLADE", "BLADE (DPMM)", "DBSCAN", "GAP", "GMM", "k-means"))

# --------------------------------------------------------
# Plots

p_ari <- ggplot(ari_df, aes(x = Method, y = ARI, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() + 
  labs(title = "ARI", x = NULL, y = "Value") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # title: larger, bold, centered
    axis.title.y = element_text(size = 14, face = "plain"),             # y-axis label: large & bold
    axis.text.x = element_text(size = 12, face = "plain", angle = 45, hjust = 1),      # x-axis text: larger
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))


p_ari <- p_ari +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     label.y = c(1.1,1.2,1.3,1.4,1.5))


p_diff <- ggplot(output_df_long, aes(x = Method, y = Diff, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +  
  labs(title = "Difference", 
       x = "", 
       y = "Value") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # title: larger, bold, centered
    axis.title.y = element_text(size = 14, face = "plain"),             # y-axis label: large & bold
    axis.text.x = element_text(size = 12, face = "plain", angle = 45, hjust = 1),      # x-axis text: larger
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))

p_diff <- p_diff +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     label.y = c(10,14,18,22,26))


p_cor <- ggplot(cor_df, aes(x = Method, y = COR, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +
  labs(title = "COR", 
       x = "", y = "Value") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # title: larger, bold, centered
    axis.title.y = element_text(size = 14, face = "plain"),             # y-axis label: large & bold
    axis.text.x = element_text(size = 12, face = "plain", angle = 45, hjust = 1),      # x-axis text: larger
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))
p_cor <- p_cor +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     method.args = list(alternative="greater"),
                     label.y = c(1.1,1.3,1.5,1.7,1.9))


p_comparison1 <- grid.arrange(p_ari, p_cor, p_diff, nrow=1)
print(p_comparison1)

# ggsave(filename="results/simulation_output/simu_comparison1.png", p_comparison1, width=2800, height=1400, units = "px")

