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


# Compare nmi, ami, fmi
nmi_df <- melt(output_df, 
               measure.vars = c("blade_nmi", "dpmm_nmi", "kmeans_nmi", "gmm_nmi", "gap_nmi", "dbscan_nmi"),
               variable.name = "Method", value.name = "NMI")
nmi_df$Method <- recode(nmi_df$Method, 
                        "blade_nmi" = "BLADE", 
                        "dpmm_nmi" = "BLADE (DPMM)",
                        "kmeans_nmi" = "k-means", 
                        "gmm_nmi" = "GMM", 
                        "gap_nmi" = "GAP", 
                        "dbscan_nmi" = "DBSCAN")

ami_df <- melt(output_df, 
               measure.vars = c("blade_ami", "dpmm_ami", "kmeans_ami", "gmm_ami", "gap_ami", "dbscan_ami"),
               variable.name = "Method", value.name = "AMI")
ami_df$Method <- recode(ami_df$Method, 
                        "blade_ami" = "BLADE", 
                        "dpmm_ami" = "BLADE (DPMM)",
                        "kmeans_ami" = "k-means", 
                        "gmm_ami" = "GMM", 
                        "gap_ami" = "GAP", 
                        "dbscan_ami" = "DBSCAN")

fmi_df <- melt(output_df, 
               measure.vars = c("blade_fmi", "dpmm_fmi", "kmeans_fmi", "gmm_fmi", "gap_fmi", "dbscan_fmi"),
               variable.name = "Method", value.name = "FMI")
fmi_df$Method <- recode(fmi_df$Method, 
                        "blade_fmi" = "BLADE", 
                        "dpmm_fmi" = "BLADE (DPMM)",
                        "kmeans_fmi" = "k-means", 
                        "gmm_fmi" = "GMM", 
                        "gap_fmi" = "GAP", 
                        "dbscan_fmi" = "DBSCAN")

# Change the order of methods (according to initial of methods)
nmi_df$Method <- factor(nmi_df$Method, levels = c("BLADE", "BLADE (DPMM)", "DBSCAN", "GAP", "GMM", "k-means"))
ami_df$Method <- factor(ami_df$Method, levels = c("BLADE", "BLADE (DPMM)", "DBSCAN", "GAP", "GMM", "k-means"))
fmi_df$Method <- factor(fmi_df$Method, levels = c("BLADE", "BLADE (DPMM)", "DBSCAN", "GAP", "GMM", "k-means"))

# --------------------------------------------------------
# Plots
p_nmi <- ggplot(nmi_df, aes(x = Method, y = NMI, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) + 
  theme_minimal() +
  labs(title = "NMI", 
       x = "", y = "Value") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # title: larger, bold, centered
    axis.title.y = element_text(size = 14, face = "plain"),             # y-axis label: large & bold
    axis.text.x = element_text(size = 12, face = "plain", angle = 45, hjust = 1),      # x-axis text: larger
    legend.position = "none"
  ) +
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
       x = "", y = "Value") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # title: larger, bold, centered
    axis.title.y = element_text(size = 14, face = "plain"),             # y-axis label: large & bold
    axis.text.x = element_text(size = 12, face = "plain", angle = 45, hjust = 1),      # x-axis text: larger
    legend.position = "none"
  ) +
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
       x = "", y = "Value") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # title: larger, bold, centered
    axis.title.y = element_text(size = 14, face = "plain"),             # y-axis label: large & bold
    axis.text.x = element_text(size = 12, face = "plain", angle = 45, hjust = 1),      # x-axis text: larger
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#F8766D", "#F564E3", "#B79F00", "#00BA38", "#00BFC4", "#619CFF"))
p_fmi <- p_fmi +
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test", 
                     paired = TRUE,         # Specify paired t-test
                     label = "p.signif",    # Label significance
                     method.args = list(alternative="greater"),
                     label.y = c(1.1,1.2,1.3,1.4,1.5))

p_comparison2 <- grid.arrange(p_nmi, p_ami, p_fmi, nrow=1)
print(p_comparison2)

# ggsave(filename="results/simulation_output/simu_comparison2.png", p_comparison2, width=2800, height=1400, units = "px")

