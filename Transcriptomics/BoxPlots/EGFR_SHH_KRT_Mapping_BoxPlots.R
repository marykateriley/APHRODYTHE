##########################################################
# MAP GENE EXPRESSION CPM                                #
# PRX GRAFT, PRX HOST, HUMAN                             #
# EGFR, SHH, KRT                                         #
# Avg 3-6 week response                                  #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 21-11-24                               #
##########################################################

setwd("/Users/mary-kateriley/Library/CloudStorage/OneDrive-Queen'sUniversityBelfast/Data Analysis/HNC/Basal PDX/Reanalysis Basal '24")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ggpubr) # add_xy_position()
library(rstatix) # wilcox_test
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(circlize)
library(clusterProfiler)
library(cowplot)
library(patchwork)


########################################################## READ DATA ##########################################################

PRX_graft <- read.csv("./0_data/edgeR/PRX/CPM/All_RemovedSamples__CPMCounts_PRX_GraftHNSCC_201124.csv", row.names = 1)
PRX_host <- read.csv("./0_data/edgeR/PRX/CPM/All_RemovedSamples__CPMCounts_PRX_HostHNSCC_201124.csv", row.names = 1)
Human <- read.csv("./0_data/edgeR/Human/CPM/All_RemovedSamples__CPMCounts_Human_HNSCC_201124.csv", row.names = 1)

colData <- read.csv("./0_data/Complete_Curated_Clinical_Information.csv", header = TRUE)

# Read in SHH GMT and EGFR ligands Files
SHH <- read.gmt("./0_data/GMT/SHH/BIOCARTA_SHH_PATHWAY.v2023.2.Hs.gmt")
SHH_subset <- SHH$gene[SHH$gene %in% c("SHH", "PTCH1", "GLI1", "GLI2", "GLI3", "SMO","SUFU")]
krt <- c("KRT8", "KRT18", "KRT5", "KRT14", "KRT1", "KRT10")
EGFR <- read.csv("./0_data/GMT/EGFR/EGFR_ligands.csv", header = T)
EGFR <- EGFR$EGF_Ligands
EGFR_extra <- c(EGFR, "CAV1", "SOX2")

SHH_mouse <- read.gmt("./0_data/GMT/SHH_Mouse/BIOCARTA_SHH_PATHWAY.v2025.1.Mm.gmt")
SHH_subset_mouse <- SHH_mouse$gene[SHH_mouse$gene %in% c("Shh", "Ptch1", "Gli1", "Gli2", "Gli3", "Smo","Sufu")]
krt_mouse <- c("Krt8", "Krt18", "Krt5", "Krt14", "Krt1", "Krt10")
EGFR_mouse <- c("Egfr","Areg","Ereg","Egf","Btc","Tgfa","Epgn","Hbegf")
EGFR_mouse_extra <- c(EGFR_mouse, "Cav1", "Sox2")

########################################################## SET VARIABLES ##########################################################

output <- "./3_output/Gene_Mapping/edgeR/BoxPlots/"
title <- "HNSCC: "
subtitle <- "Basal HNSCC PDX: Filtered Technical; Lymphoma; HPV"

# Don't use scientific notation unless the numbers are 999+ 
options(scipen = 999)


########################################################## DATA CURATION ##########################################################

# Remove any sample that does not have Definitive_Response_3wk (NAs) from colData
colData_filtered <- colData  %>%
  filter(!is.na(Definitive_Response_3_6wk))

# Ensure all samples in PRX_graft match those samples that have response
PRX_graft <- PRX_graft[, substr(colnames(PRX_graft), 1, 7) %in% colData_filtered$Case_ID]


# Ensure all samples in PRX_host match those samples that have response
PRX_host <- PRX_host[, substr(colnames(PRX_host), 1, 7) %in% colData_filtered$Case_ID]

# Ensure all samples in Human match those samples that have response
Human <- Human[, substr(colnames(Human), 1, 7) %in% colData_filtered$Case_ID]




## RESPONDERS + NONRESPONDERS
## GRAFT ##
# Extract the Case_IDs from PRX_graft column names
prx_graft_case_ids <- substr(colnames(PRX_graft), 1, 7)

# Match the Case_IDs with the Definitive_Response_3_6wk in colData_filtered
response_data <- colData_filtered %>%
  filter(Case_ID %in% prx_graft_case_ids) %>%
  select(Case_ID, Definitive_Response_3_6wk)

# Repeat the response data for each sample in PRX_graft
response_data_repeated <- response_data[rep(1:nrow(response_data), times = table(factor(prx_graft_case_ids, levels = response_data$Case_ID))),]

# Add the response category for each sample in PRX_graft
response_data_repeated$Sample <- prx_graft_case_ids

# Count the number of Responders and Non-Responders in PRX_graft
graft_resp <- sum(response_data_repeated$Definitive_Response_3_6wk == "Responder")
graft_nonresp <- sum(response_data_repeated$Definitive_Response_3_6wk == "NonResponder")

# Display the counts for verification
graft_resp
graft_nonresp

rm(response_data)
rm(response_data_repeated)
rm(prx_graft_case_ids)


## HOST ##
# Extract the Case_IDs from PRX_graft column names
prx_host_case_ids <- substr(colnames(PRX_host), 1, 7)

# Match the Case_IDs with the Definitive_Response_3_6wk in colData_filtered
response_data <- colData_filtered %>%
  filter(Case_ID %in% prx_host_case_ids) %>%
  select(Case_ID, Definitive_Response_3_6wk)

# Repeat the response data for each sample in PRX_graft
response_data_repeated <- response_data[rep(1:nrow(response_data), times = table(factor(prx_host_case_ids, levels = response_data$Case_ID))),]

# Add the response category for each sample in PRX_graft
response_data_repeated$Sample <- prx_host_case_ids

# Count the number of Responders and Non-Responders in PRX_graft
host_resp <- sum(response_data_repeated$Definitive_Response_3_6wk == "Responder")
host_nonresp <- sum(response_data_repeated$Definitive_Response_3_6wk == "NonResponder")

# Display the counts for verification
host_resp
host_nonresp

rm(response_data)
rm(response_data_repeated)
rm(prx_host_case_ids)


## HUMAN ##
# Extract the Case_IDs from PRX_graft column names
human_case_ids <- substr(colnames(Human), 1, 7)

# Match the Case_IDs with the Definitive_Response_3_6wk in colData_filtered
response_data <- colData_filtered %>%
  filter(Case_ID %in% human_case_ids) %>%
  select(Case_ID, Definitive_Response_3_6wk)

# Repeat the response data for each sample in PRX_graft
response_data_repeated <- response_data[rep(1:nrow(response_data), times = table(factor(human_case_ids, levels = response_data$Case_ID))),]

# Add the response category for each sample in PRX_graft
response_data_repeated$Sample <- human_case_ids

# Count the number of Responders and Non-Responders in PRX_graft
human_resp <- sum(response_data_repeated$Definitive_Response_3_6wk == "Responder")
human_nonresp <- sum(response_data_repeated$Definitive_Response_3_6wk == "NonResponder")

rm(response_data)
rm(response_data_repeated)
rm(human_case_ids)

# Display the counts for verification
human_resp
human_nonresp



########################################################## BOX PLOTS - MAP EXPRESSION SHH GRAFT ##########################################################

## SHH  genes in normalised counts matrix ##
gene_shh <- intersect(SHH$gene, rownames(PRX_graft))
SHH_GEX <- PRX_graft[gene_shh, ]

# Convert to long format and filter columns
shh_long <- as.data.frame(SHH_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>% 
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

p_values <- shh_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance))

p_values <- p_values %>% 
  add_xy_position(x = "Gene")

resp <- shh_long %>% 
  group_by(Gene, Definitive_Response_3_6wk) %>% 
  summarise(Count = n(), .groups = "drop")


plot <- ggplot(shh_long, aes(x = factor(Gene), y = Expression)) +
  geom_boxplot(aes(fill = Definitive_Response_3_6wk), 
               width = 0.5, 
               colour = "darkgrey", 
               position = position_dodge(width = 0.9), 
               alpha = 0.7) + 
  labs(title = paste0("SHH Biocarta Genes: PDX Graft"),
       subtitle = paste0(subtitle, "- Wilcox Test"),
       x = "Gene",
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", graft_resp, ")"),
                               "NonResponder" = paste0("NonResponder (", graft_nonresp, ")"))) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_SHH_PDX_Graft.pdf"), plot = plot, width = 11, height = 14)


# Extract legend from your plot
legend <- get_legend(plot)
ggsave(paste0(output, "Legend_Response.pdf"), plot = legend, width = 3, height = 2)


########################################################## BOX PLOTS PANEL - MAP EXPRESSION SHH_subset GRAFT ##########################################################

## SHH_subset  genes in normalised counts matrix ##
gene_shh <- intersect(SHH_subset, rownames(PRX_graft))
SHH_GEX <- PRX_graft[gene_shh, ]

# Convert to long format and filter columns
shh_long <- as.data.frame(SHH_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>% 
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

# Calculate max expression per gene
max_expression_per_gene <- shh_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

# Adjust p-values dataframe to use max expression for y.position
p_values <- shh_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>%  # Join max expression
  mutate(y.position = max_expression * 1.1)

resp <- shh_long %>% 
  group_by(Gene, Definitive_Response_3_6wk) %>% 
  summarise(Count = n(), .groups = "drop")


plot <- ggplot(shh_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ Gene, scales = "free_y") +  # Facet by gene
  labs(title = paste0("SHH Genes: PDX Graft"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", graft_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", graft_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14),  # Facet strip labels
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_SHHsubset_PDX_Graft.pdf"), plot = plot, width = 11, height = 14)

rm(gene_shh, SHH_GEX, shh_long)


########################################################## BOX PLOTS - MAP EXPRESSION EGFR EXTRA GRAFT ##########################################################

## EGFR  ligands in normalised counts matrix ##
gene_egfr <- intersect(EGFR_extra, rownames(PRX_graft))
EGFR_GEX <- PRX_graft[gene_egfr, ]

# Convert to long format and filter columns
egfr_long <- as.data.frame(EGFR_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>%
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

egfr_long$Gene <- factor(egfr_long$Gene, levels = c("EGFR", "AREG", "BTC", "EGF", "EPGN", "ERBB1", "EREG", "HBEGF", "HER1", "TGFA", "CAV1", "SOX2"))

# Calculate max expression per gene
max_expression_per_gene <- egfr_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

#Adjust p-values dataframe to use max expression for y.position
p_values <- egfr_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>%  # Join max expression
  mutate(y.position = max_expression * 1.1)

plot <- ggplot(egfr_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) + 
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = paste0("EGFR, Ligands CAV1 and SOX2: PDX Graft"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", graft_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", graft_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_EGFRExtraPanel_PDX_Graft.pdf"), plot = plot, width = 11, height = 14)

rm(gene_egfr, EGFR_GEX, egfr_long)


########################################################## BOX PLOTS - MAP EXPRESSION KRT Markers GRAFT ##########################################################

## KRT Marters in normalised counts matrix ##
gene_krt <- intersect(krt, rownames(PRX_graft))
KRT_GEX <- PRX_graft[gene_krt, ]

# Convert to long format and filter columns
krt_long <- as.data.frame(KRT_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>%
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

krt_long$Gene <- factor(krt_long$Gene, levels = c("KRT1", "KRT5", "KRT8", "KRT10", "KRT14", "KRT18"))

p_values <- krt_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance))

p_values <- p_values %>% 
  add_xy_position(x = "Gene")


plot <- ggplot(krt_long, aes(x = factor(Gene), y = Expression)) +
  geom_boxplot(aes(fill = Definitive_Response_3_6wk), 
               width = 0.5, 
               colour = "darkgrey", 
               position = position_dodge(width = 0.9), 
               alpha = 0.7) + 
  labs(
    title = paste0("Keratin Markers: PDX Graft"),
    subtitle = paste0(subtitle, " - Wilcox Test"),
    x = "Gene",
    y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", graft_resp, ")"),
                               "NonResponder" = paste0("NonResponder (", graft_nonresp, ")"))) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_KRT_PDX_Graft.pdf"), plot = plot, width = 11, height = 14)




########################################################## BOX PLOTS PANEL - MAP EXPRESSION KRT MARKERS GRAFT ##########################################################
# Calculate max expression per gene
max_expression_per_gene <- krt_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

#Adjust p-values dataframe to use max expression for y.position
p_values <- krt_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>% 
  mutate(y.position = max_expression * 1.1)


plot <- ggplot(krt_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) + 
  facet_wrap(~ Gene, ncol = 3, scales = "free_y") +
  labs(title = paste0("Keratin Markers: PDX Graft"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", graft_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", graft_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_KRTPanel_PDX_Graft.pdf"), plot = plot, width = 11, height = 14)

rm(gene_krt, KRT_GEX, krt_long)

rm(PRX_graft, graft_resp, graft_nonresp)

########################################################## BOX PLOTS - MAP EXPRESSION SHH HOST ##########################################################

## SHH  genes in normalised counts matrix ##
gene_shh <- intersect(SHH_mouse$gene, rownames(PRX_host))
SHH_GEX <- PRX_host[gene_shh, ]

# Convert to long format and filter columns
shh_long <- as.data.frame(SHH_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>%
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

p_values <- shh_long %>%
  group_by(Gene) %>%
  filter(!is.na(Expression) & !is.na(Definitive_Response_3_6wk)) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>%
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance))
p_values <- p_values %>% 
  add_xy_position(x = "Gene")

resp <- shh_long %>% 
  group_by(Gene, Definitive_Response_3_6wk) %>% 
  summarise(Count = n(), .groups = "drop")


plot <- ggplot(shh_long, aes(x = factor(Gene), y = Expression)) +
  geom_boxplot(aes(fill = Definitive_Response_3_6wk), 
               width = 0.5, 
               colour = "darkgrey", 
               position = position_dodge(width = 0.9), 
               alpha = 0.7) + 
  labs(
    title = paste0("SHH Biocarta Genes: PDX Host"),
    subtitle = paste0(subtitle, "- Wilcox Test"),
    x = "Gene",
    y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", host_resp, ")"),
                               "NonResponder" = paste0("NonResponder (", host_nonresp, ")"))) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_SHH_PDX_Host.pdf"), plot = plot, width = 11, height = 14)




########################################################## BOX PLOTS PANEL - MAP EXPRESSION SHH_subset HOST ##########################################################

## SHH_subset  genes in normalised counts matrix ##
gene_shh <- intersect(SHH_subset_mouse, rownames(PRX_host))
SHH_GEX <- PRX_host[gene_shh, ]

# Convert to long format and filter columns
shh_long <- as.data.frame(SHH_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>%
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

# Calculate max expression per gene
max_expression_per_gene <- shh_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))


#Adjust p-values dataframe to use max expression for y.position
p_values <- shh_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>%
  mutate(y.position = max_expression * 1.1)

resp <- shh_long %>% 
  group_by(Gene, Definitive_Response_3_6wk) %>% 
  summarise(Count = n(), .groups = "drop")


plot <- ggplot(shh_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = paste0("SHH Biocarta Genes: PDX Host"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", host_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", host_nonresp, " samples)")
                    )) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14), 
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_SHHsubset_PDX_Host.pdf"), plot = plot, width = 11, height = 14)

rm(gene_shh, SHH_GEX, shh_long)



########################################################## BOX PLOTS - MAP EXPRESSION EGFR EXTRA GRAFT ##########################################################

## EGFR  ligands in normalised counts matrix ##
gene_egfr <- intersect(EGFR_mouse_extra, rownames(PRX_host))
EGFR_GEX <- PRX_host[gene_egfr, ]

# Convert to long format and filter columns
egfr_long <- as.data.frame(EGFR_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>%
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

egfr_long$Gene <- factor(egfr_long$Gene, levels = c("Egfr", "Areg", "Btc", "Egf", "Epgn", "Ereg", "Hbegf", "Tgfa", "Cav1", "Sox2"))

# Calculate max expression per gene
max_expression_per_gene <- egfr_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

#Adjust p-values dataframe to use max expression for y.position
p_values <- egfr_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>%  # Join max expression
  mutate(y.position = max_expression * 1.1)

plot <- ggplot(egfr_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = paste0("EGFR, Ligands Cav1 and Sox2: PDX Host"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", host_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", host_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14), 
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_EGFRExtraPanel_PDX_Host.pdf"), plot = plot, width = 11, height = 14)

rm(gene_egfr, EGFR_GEX, egfr_long)



########################################################## BOX PLOTS - MAP EXPRESSION KRT Markers HOST ##########################################################

## KRT Marters in normalised counts matrix ##
gene_krt <- intersect(krt_mouse, rownames(PRX_host))
KRT_GEX <- PRX_host[gene_krt, ]

# Convert to long format and filter columns
krt_long <- as.data.frame(KRT_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>% 
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)  
krt_long$Gene <- factor(krt_long$Gene, levels = c("Krt1", "Krt5", "Krt8", "Krt10", "Krt14", "Krt18"))

p_values <- krt_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance))

p_values <- p_values %>% 
  add_xy_position(x = "Gene")


plot <- ggplot(krt_long, aes(x = factor(Gene), y = Expression)) + 
  geom_boxplot(aes(fill = Definitive_Response_3_6wk), 
               width = 0.5, 
               colour = "darkgrey", 
               position = position_dodge(width = 0.9), 
               alpha = 0.7) + 
  labs(
    title = paste0("Keratin Markers: PDX Host"),
    subtitle = paste0(subtitle, " - Wilcox Test"),
    x = "Gene",
    y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", host_resp, ")"),
                               "NonResponder" = paste0("NonResponder (", host_nonresp, ")"))) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_KRT_PDX_Host.pdf"), plot = plot, width = 11, height = 14)




########################################################## BOX PLOTS PANEL - MAP EXPRESSION KRT MARKERS HOST ##########################################################

# Calculate max expression per gene
max_expression_per_gene <- krt_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

#Adjust p-values dataframe to use max expression for y.position
p_values <- krt_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>% 
  mutate(y.position = max_expression * 1.1)




plot <- ggplot(krt_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ Gene, ncol = 3, scales = "free_y") + 
  labs(title = paste0("Keratin Markers: PDX Host"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", host_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", host_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14), 
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_KRTPanel_PDX_Host.pdf"), plot = plot, width = 11, height = 14)

rm(gene_krt, KRT_GEX, krt_long)

rm(PRX_host, host_resp, host_nonresp)








########################################################## BOX PLOTS - MAP EXPRESSION SHH HUMAN ##########################################################

## SHH  genes in normalised counts matrix ##
gene_shh <- intersect(SHH$gene, rownames(Human))
SHH_GEX <- Human[gene_shh, ]

# Convert to long format and filter columns
shh_long <- as.data.frame(SHH_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>% 
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

p_values <- shh_long %>%
  group_by(Gene) %>%
  filter(n_distinct(Expression) > 1) %>%  # Only keep groups with more than one unique value for Expression  
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance))

p_values <- p_values %>% 
  add_xy_position(x = "Gene")

resp <- shh_long %>% 
  group_by(Gene, Definitive_Response_3_6wk) %>% 
  summarise(Count = n(), .groups = "drop")


plot <- ggplot(shh_long, aes(x = factor(Gene), y = Expression)) +
  geom_boxplot(aes(fill = Definitive_Response_3_6wk), 
               width = 0.5, 
               colour = "darkgrey", 
               position = position_dodge(width = 0.9), 
               alpha = 0.7) + 
  labs(
    title = paste0("SHH Biocarta Genes: Human"),
    subtitle = paste0(subtitle, "- Wilcox Test"),
    x = "Gene",
    y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", human_resp, ")"),
                               "NonResponder" = paste0("NonResponder (", human_nonresp, ")"))) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_SHH_Human.pdf"), plot = plot, width = 11, height = 14)




########################################################## BOX PLOTS PANEL - MAP EXPRESSION SHH_subset Human ##########################################################

## SHH_subset  genes in normalised counts matrix ##
gene_shh <- intersect(SHH_subset, rownames(Human))
SHH_GEX <- Human[gene_shh, ]

# Convert to long format and filter columns
shh_long <- as.data.frame(SHH_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>%
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)

# Calculate max expression per gene
max_expression_per_gene <- shh_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

#Adjust p-values dataframe to use max expression for y.position
p_values <- shh_long %>%
  group_by(Gene) %>%
  filter(n_distinct(Expression) > 1) %>%  # Only keep groups with more than one unique value for Expression
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>%
  mutate(y.position = max_expression * 1.1)

resp <- shh_long %>% 
  group_by(Gene, Definitive_Response_3_6wk) %>% 
  summarise(Count = n(), .groups = "drop")


plot <- ggplot(shh_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = paste0("SHH Biocarta Genes: Human"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", human_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", human_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_SHHsubset_Human.pdf"), plot = plot, width = 11, height = 14)

rm(gene_shh, SHH_GEX, shh_long)




########################################################## BOX PLOTS - MAP EXPRESSION EGFR EXTRA HUMAN ##########################################################

## EGFR  ligands in normalised counts matrix ##
gene_egfr <- intersect(EGFR_extra, rownames(Human))
EGFR_GEX <- Human[gene_egfr, ]

# Convert to long format and filter columns
egfr_long <- as.data.frame(EGFR_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>% # Optional: Remove the temporary key column if no longer needed
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)  # Keep only the desired columns

egfr_long$Gene <- factor(egfr_long$Gene, levels = c("EGFR", "AREG", "BTC", "EGF", "EPGN", "ERBB1", "EREG", "HBEGF", "HER1", "TGFA", "CAV1", "SOX2"))

# Calculate max expression per gene
max_expression_per_gene <- egfr_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

#Adjust p-values dataframe to use max expression for y.position
p_values <- egfr_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>%
  mutate(y.position = max_expression * 1.1)

plot <- ggplot(egfr_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ Gene, scales = "free_y") + 
  labs(title = paste0("EGFR, Ligands CAV1 and SOX2: PDX Human"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", human_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", human_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_EGFRExtraPanel_PDX_Human.pdf"), plot = plot, width = 11, height = 14)

rm(gene_egfr, EGFR_GEX, egfr_long)



########################################################## BOX PLOTS - MAP EXPRESSION KRT Markers Human ##########################################################

## KRT Marters in normalised counts matrix ##
gene_krt <- intersect(krt, rownames(Human))
KRT_GEX <- Human[gene_krt, ]

# Convert to long format and filter columns
krt_long <- as.data.frame(KRT_GEX) %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Expression") %>% 
  mutate(Sample_Key = substr(Sample, 1, 7)) %>%  # Extract the first 7 characters of Sample
  left_join(as.data.frame(colData_filtered), by = c("Sample_Key" = "Case_ID")) %>% 
  dplyr::select(-Sample_Key)  %>%
  mutate(Definitive_Response_3_6wk = factor(Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))) %>%
  dplyr::select(Gene, Definitive_Response_3_6wk, Expression)
krt_long$Gene <- factor(krt_long$Gene, levels = c("KRT1", "KRT5", "KRT8", "KRT10", "KRT14", "KRT18"))

p_values <- krt_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance))

p_values <- p_values %>% 
  add_xy_position(x = "Gene")


plot <- ggplot(krt_long, aes(x = factor(Gene), y = Expression)) + 
  geom_boxplot(aes(fill = Definitive_Response_3_6wk), 
               width = 0.5, 
               colour = "darkgrey", 
               position = position_dodge(width = 0.9), 
               alpha = 0.7) + 
  labs(
    title = paste0("Keratin Markers: Human"),
    subtitle = paste0(subtitle, " - Wilcox Test"),
    x = "Gene",
    y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", human_resp, ")"),
                               "NonResponder" = paste0("NonResponder (", human_nonresp, ")"))) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_KRT_Human.pdf"), plot = plot, width = 11, height = 14)




########################################################## BOX PLOTS PANEL - MAP EXPRESSION KRT MARKERS Human ##########################################################

# Calculate max expression per gene
max_expression_per_gene <- krt_long %>%
  group_by(Gene) %>%
  summarise(max_expression = max(Expression, na.rm = TRUE))

#Adjust p-values dataframe to use max expression for y.position
p_values <- krt_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Definitive_Response_3_6wk, p.adjust.method = "BH") %>% 
  mutate(significance = case_when(p < 0.001 ~ "***",
                                  p < 0.01 ~ "**",
                                  p < 0.05 ~ "*",
                                  TRUE ~ "ns"),
         p_label = paste0("p = ", signif(p, 3), " ", significance)) %>% 
  left_join(max_expression_per_gene, by = "Gene") %>%  # Join max expression
  mutate(y.position = max_expression * 1.1)




plot <- ggplot(krt_long, aes(x = Definitive_Response_3_6wk, y = Expression, fill = Definitive_Response_3_6wk)) +  
  geom_boxplot(width = 0.5, colour = "darkgrey", position = position_dodge(width = 0.9), alpha = 0.7,
               outlier.shape = NA,
               outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ Gene, ncol = 3, scales = "free_y") +
  labs(title = paste0("Keratin Markers: Human"),
       subtitle = paste0(subtitle, " - Wilcox Test"),
       x = NULL,
       y = "Counts Per Milliom (CPM)") +
  scale_fill_manual(values = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
                    labels = c("Responder" = paste0("Responder (", human_resp, " samples)"),
                               "NonResponder" = paste0("NonResponder (", human_nonresp, " samples)"))) + 
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.position = "bottom") +
  stat_pvalue_manual(p_values, label = "p_label", tip.length = 0.02,
                     y.position = p_values$y.position, bracket.shorten = 0.05,
                     size = 4,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0(output, "BoxPlot_KRTPanel_Human.pdf"), plot = plot, width = 11, height = 14)

rm(gene_krt, KRT_GEX, krt_long)


rm(Human, human_resp, human_nonresp)
