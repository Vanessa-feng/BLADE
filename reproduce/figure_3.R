setwd("C:/Xin/oral_cancer/Code_Xin/BLADE_git")

library(dplyr)
library(reshape2)
library(ggpubr)
library(gridExtra)

# Load the results of 200 simulations for 6 different methods
output_df <- read.csv(file = paste0(getwd(), "/results/simulation_output/simulation_comparison.csv"))

# ----------------------------------------------------
# Compare ARI
ari_df <- melt(output_df, 
               measure.vars = c("blade_ari", "dpmm_ari", "kmeans_ari", "gmm_ari", "gap_ari", "dbscan_ari"),
               variable.name = "Method", 
               value.name = "ARI")
ari_df$Method <- recode(ari_df$Method, 
                        "blade_ari" = "BLADE(MFM)", 
                        "dpmm_ari" = "BLADE(DPMM)",
                        "kmeans_ari" = "K-means", 
                        "gmm_ari" = "GMM", 
                        "gap_ari" = "GAP", 
                        "dbscan_ari" = "DBSCAN")

p_ari <- ggplot(ari_df, aes(x = Method, y = ARI, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() + 
  labs(title = "ARI", 
       x = "", 
       y = "Adjusted Rand Index (ARI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))

# Add significance lines and labels to the plot
my_comparisons <- list( c("BLADE(MFM)", "BLADE(DPMM)"), c("BLADE(MFM)", "K-means"), c("BLADE(MFM)", "GMM"), c("BLADE(MFM)", "GAP"), 
                        c("BLADE(MFM)", "DBSCAN"))
p_ari <- p_ari +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     label.y = c(1.1,1.2,1.3,1.4,1.5))

# ----------------------------------------------------
# running time
time_df <- melt(output_df, 
                measure.vars = c("blade_time", "dpmm_time", "kmeans_time", "gmm_time", "gap_time", "dbscan_time"),
                variable.name = "Method", 
                value.name = "Running_time")
time_df$Method <- recode(time_df$Method, 
                         "blade_time" = "BLADE(MFM)", 
                         "dpmm_time" = "BLADE(DPMM)",
                         "kmeans_time" = "K-means", 
                         "gmm_time" = "GMM", 
                         "gap_time" = "GAP", 
                         "dbscan_time" = "DBSCAN")

mean_time_df <- time_df %>%
  group_by(Method) %>%
  summarize(Mean_Time = mean(Running_time, na.rm = TRUE))


p_time <- ggplot(time_df, aes(x = Method, y = Running_time, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() + 
  labs(title = "Running time", 
       x = "", 
       y = "Running time") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") + 
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))

p_time <- p_time +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     method.args = list(alternative="greater"),
                     label.y = c(1500,1600,1700,1800,1900))

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
                                "blade_numlayer" = "BLADE(MFM)", 
                                "dpmm_numlayer" = "BLADE(DPMM)",
                                "kmeans_numlayer" = "K-means", 
                                "gmm_numlayer" = "GMM", 
                                "gap_numlayer" = "GAP", 
                                "dbscan_numlayer" = "DBSCAN")

p_diff <- ggplot(output_df_long, aes(x = Method, y = Diff, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +  
  labs(title = "Difference", 
       x = "", 
       y = "Difference") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))
p_diff <- p_diff +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     label.y = c(10,14,18,22,26))

t.test(output_df_long[output_df_long$Method=="BLADE(MFM)", "Diff"], 
       output_df_long[output_df_long$Method=="BLADE(DPMM)", "Diff"], paired=T, alternative="less")

# --------------------------------------------------------------------------------
# Compare nmi, ami, fmi and COR
nmi_df <- melt(output_df, 
               measure.vars = c("blade_nmi", "dpmm_nmi", "kmeans_nmi", "gmm_nmi", "gap_nmi", "dbscan_nmi"),
               variable.name = "Method", value.name = "NMI")
nmi_df$Method <- recode(nmi_df$Method, 
                        "blade_nmi" = "BLADE(MFM)", 
                        "dpmm_nmi" = "BLADE(DPMM)",
                        "kmeans_nmi" = "K-means", 
                        "gmm_nmi" = "GMM", 
                        "gap_nmi" = "GAP", 
                        "dbscan_nmi" = "DBSCAN")

ami_df <- melt(output_df, 
               measure.vars = c("blade_ami", "dpmm_ami", "kmeans_ami", "gmm_ami", "gap_ami", "dbscan_ami"),
               variable.name = "Method", value.name = "AMI")
ami_df$Method <- recode(ami_df$Method, 
                        "blade_ami" = "BLADE(MFM)", 
                        "dpmm_ami" = "BLADE(DPMM)",
                        "kmeans_ami" = "K-means", 
                        "gmm_ami" = "GMM", 
                        "gap_ami" = "GAP", 
                        "dbscan_ami" = "DBSCAN")

fmi_df <- melt(output_df, 
               measure.vars = c("blade_fmi", "dpmm_fmi", "kmeans_fmi", "gmm_fmi", "gap_fmi", "dbscan_fmi"),
               variable.name = "Method", value.name = "FMI")
fmi_df$Method <- recode(fmi_df$Method, 
                        "blade_fmi" = "BLADE(MFM)", 
                        "dpmm_fmi" = "BLADE(DPMM)",
                        "kmeans_fmi" = "K-means", 
                        "gmm_fmi" = "GMM", 
                        "gap_fmi" = "GAP", 
                        "dbscan_fmi" = "DBSCAN")

cor_df <- melt(output_df, 
               measure.vars = c("blade_cor", "dpmm_cor", "kmeans_cor", "gmm_cor", "gap_cor", "dbscan_cor"),
               variable.name = "Method", value.name = "COR")
cor_df$Method <- recode(cor_df$Method, 
                        "blade_cor" = "BLADE(MFM)", 
                        "dpmm_cor" = "BLADE(DPMM)",
                        "kmeans_cor" = "K-means", 
                        "gmm_cor" = "GMM", 
                        "gap_cor" = "GAP", 
                        "dbscan_cor" = "DBSCAN")

p_nmi <- ggplot(nmi_df, aes(x = Method, y = NMI, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +
  labs(title = "NMI", 
       x = "", y = "Normalized Mutual Information (NMI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))
p_nmi <- p_nmi +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     method.args = list(alternative="greater"),
                     label.y = c(1.1,1.2,1.3,1.4,1.5))

p_ami <- ggplot(ami_df, aes(x = Method, y = AMI, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +
  labs(title = "AMI", 
       x = "", y = "Adjusted Mutual Information (AMI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))
p_ami <- p_ami +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     method.args = list(alternative="greater"),
                     label.y = c(1.1,1.2,1.3,1.4,1.5))

p_fmi <- ggplot(fmi_df, aes(x = Method, y = FMI, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +
  labs(title = "FMI", 
       x = "", y = "Fowlkes-Mallows Index (FMI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))
p_fmi <- p_fmi +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     method.args = list(alternative="greater"),
                     label.y = c(1.1,1.2,1.3,1.4,1.5))

p_cor <- ggplot(cor_df, aes(x = Method, y = COR, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +
  labs(title = "COR", 
       x = "", y = "Correlation (COR)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")+ 
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))
p_cor <- p_cor +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     method.args = list(alternative="greater"),
                     label.y = c(1.1,1.3,1.5,1.7,1.9))


p_comparison1 <- grid.arrange(p_ari, p_cor, p_diff, nrow=1)
p_comparison2 <- grid.arrange(p_nmi, p_ami, p_fmi, nrow=1)

print(p_comparison1)
print(p_comparison2)