##########################################################
# MutationalSignatures - STATS                           #
# Bar charts                                             #
# compare PDX                                           #
# RESPONDERS V NONRESPONDERS                             #
#                                                        #
# Filtered out HPV/Lymohoma                              #
#                                                        #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 21-11-25                               #
##########################################################

setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(dplyr)
library(tidyr)
#BiocManager::install("maftools")
library(maftools)
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(forcats)
library(rstatix)
library(ggpubr)
library(stringr)
library(scales)


########################################################## READ DATA ##########################################################

resp_cosm <- readRDS("./2_pipelines/MAFtools/PDX_RvNR/Filtered_HPV_Lymphoma_VAF0.1/HN0039A_Filtered_PDX_RvNR__Responders_cosmic_match.rds")
resp_cosm

nonresp_cosm <- readRDS("./2_pipelines/MAFtools/PDX_RvNR/Filtered_HPV_Lymphoma_VAF0.1/HN0039A_Filtered_PDX_RvNR__NonResponders_cosmic_match.rds")
nonresp_cosm

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/NEWMutationalPatterns/PDX_Resp/"
output <- "./3_output/NEWMutationalPatterns/PDX_Resp/"
dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

title <- "HN0039A_Filtered_PDX_NRvR"

resp_SBS <- c("SBS13", "SBS18", "SBS5", "SBS7b")
nonresp_SBS <- c("SBS6", "SBS40", "SBS13")
target_sigs <- unique(c(resp_SBS, nonresp_SBS))

########################################################## SET VARIABLES ##########################################################
## Long data from cosine similarity matrices
resp_mat <- as.data.frame(resp_cosm$cosine_similarities)
nonresp_mat <- as.data.frame(nonresp_cosm$cosine_similarities)

resp_long <- resp_mat %>%
  rownames_to_column("DeNovo") %>%
  pivot_longer(-DeNovo, names_to = "COSMIC", values_to = "cosine_similarities") %>%
  mutate(Cohort = "Responders")

nonresp_long <- nonresp_mat %>%
  rownames_to_column("DeNovo") %>%
  pivot_longer(-DeNovo, names_to = "COSMIC", values_to = "cosine_similarities") %>%
  mutate(Cohort = "NonResponders")

cosmic_combined <- bind_rows(resp_long, nonresp_long) %>%
  mutate(Cohort = factor(Cohort, levels = c("Responders","NonResponders")))

## COSMIC signatures to compare
cosmic_subset <- cosmic_combined %>%
  filter(COSMIC %in% target_sigs) %>%
  mutate(COSMIC = factor(COSMIC, levels = target_sigs))

## Stats per COSMIC
p_values <- cosmic_subset %>%
  group_by(COSMIC) %>%
  filter(n_distinct(Cohort) == 2) %>% 
  wilcox_test(cosine_similarities ~ Cohort, exact = FALSE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "COSMIC", data = cosmic_subset) %>%
  mutate(COSMIC = factor(COSMIC, levels = target_sigs)) %>% 
  mutate(label = paste0("p = ", signif(p.adj, 3)," (", p.adj.signif, ")"))

## y positions just above the max per COSMIC
y_top <- cosmic_subset %>% group_by(COSMIC) %>%
  summarise(ymax = max(cosine_similarities, na.rm = TRUE), .groups = "drop")
p_values <- p_values %>%
  left_join(y_top, by = "COSMIC") %>%
  mutate(y.position = pmax(y.position, ymax * 1.05, na.rm = TRUE)) %>% 
  mutate(label = paste0("p = ", signif(p.adj, 3)," (", p.adj.signif, ")"))

dodge <- position_dodge(width = 0.8)
p <- ggplot(cosmic_subset, aes(x = COSMIC, y = cosine_similarities, fill = Cohort)) +
  geom_boxplot(position = dodge, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
              alpha = 0.6, size = 1.2) +
  scale_fill_manual(values = c("Responders" = "#1B7837FF", "NonResponders" = "#762A83FF")) +
  stat_pvalue_manual(data = p_values,
                     label = "label",
                     xmin = "xmin", xmax = "xmax",
                     y.position = "y.position",
                     tip.length = 0.01, bracket.size = 0.6, size = 4,
                     inherit.aes = FALSE) +
  labs(x = "COSMIC SBS",
       y = "Cosine similarity",
       fill = "Cohort",
       title = "PDX Nonresponders vs Responders: Cosine Similarity to Selected COSMIC signatures") +
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold"))
ggsave(file.path(output, paste0(title, "_Top6_Signatures_PDX_NRvR.pdf")), p, width = 11, height = 6)
print(p)



