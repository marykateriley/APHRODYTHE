##########################################################
# ssGSEA (GSVA package) TPM                              #
# PDX + HUMAN                                            #
# HALLMARKS                                              #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 19-11-2025                             #
##########################################################

setwd("/Users/mary-kateriley/Library/CloudStorage/OneDrive-Queen'sUniversityBelfast/Data Analysis/HNC/Basal PDX/Reanalysis Basal '24")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(broom)
library(clusterProfiler)
library(KEGGREST)
library(msigdbr)

library(GSVA)
library(decoupleR)
library(viper)

library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(reshape2) # Melting the data frame for ggplot2

library(ggExtra)
library(ggpubr)
library(rstatix) # t-test




########################################################## READ DATA ##########################################################
# Load filtered Graft Raw Counts Matrix - from DESeq_Data_Input.R
# As ssGSEA requires integers 
# Read in GEX for GSVA
GEX <- read.csv("./0_data/edgeR/PDXHuman/TPM/PDX_Human_TPM_GEX.csv", row.names = 1)
rownames(GEX) <- toupper(rownames(GEX))
GEX <- log2(GEX + 1)

metadata <- read.csv("./2_pipelines/PDX_Human/DESeq2/Unfiltered/PDXandHuman_Merged_colData_Filtered.csv", header = T, row.names = 1)

H_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

network <- H_t2g %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
  distinct(gs_name, gene_symbol) %>%
  rename(GeneSet = gs_name, Gene = gene_symbol)

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/PDX_Human/Enrichment_Analysis/ssGSEA/TPM_Hallmarks/"
output <- "./3_output/PDX_Human/Enrichment_Analysis/ssGSEA/TPM_Hallmarks/"
name <- "ssGSEA"
title <- "PDX HNSCC: ssGSEA scores Hallmarks"

dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

source("./1_code/RFunctions/new_heatmap.R")

########################################################## DATA CURATION ##########################################################

# Check which genes from the network are NOT in your data
missing_genes <- setdiff(unique(network$Gene), rownames(GEX))
print(missing_genes)

########################################################## ssGSEA DECOUPLER ##########################################################

# decoupler can allow ssGSEA calculation through the GSVA package! SET METHOD = ssGSEA !!
ssgsea <- run_gsva(
  mat = GEX,
  network = network,
  .source = "GeneSet",
  .target = "Gene",
  verbose = FALSE,
  method = "ssgsea")

# format into df
ssgsea_wide <- ssgsea %>%
  dplyr::select(source, condition, score) %>%
  tidyr::pivot_wider(names_from = source,
                     values_from = score) %>%
  tibble::column_to_rownames("condition")
ssgsea_1 <- as.data.frame(ssgsea_wide)   # samples x gene sets

#save single sample activity predictions
write.csv(ssgsea_1, file = paste0(pipelines,"ssGSEAScores_PDXHuman.csv"))

# Combine ssGSEA scores with metadata
ssgsea_1$SampleID <- rownames(ssgsea_1)

# Are all samples in colnames of GEX present rownames of colData and are they in the same order
all_samples_present <- all(rownames(metadata) %in% rownames(ssgsea_1))
are_same_order <- identical(rownames(metadata), rownames(ssgsea_1))
metadata <- metadata[match(rownames(ssgsea_1), rownames(metadata)), ] # order
all_samples_present <- all(rownames(metadata) %in% rownames(ssgsea_1)) # check again
are_same_order <- identical(rownames(metadata), rownames(ssgsea_1))


combined_data <- merge(ssgsea_1, metadata, by.x = "SampleID", by.y = "Sample_ID")
colnames(combined_data)
ssgsea_1$SampleID <- NULL # Remove SampleID column


########################################################## T- test ##########################################################

# Grouping column
group_col <- "Model"

# Identify ssGSEA columns
gs_cols <- setdiff(colnames(combined_data), c(colnames(metadata), "SampleID"))

# Long format for stats/plots
long_ssgsea <- combined_data %>%
  dplyr::select(SampleID, all_of(group_col), all_of(gs_cols)) %>%
  tidyr::pivot_longer(cols = all_of(gs_cols),
                      names_to = "GeneSet", values_to = "Score") %>%
  tidyr::drop_na(Score)

# Stats
ssgsea_stats <- long_ssgsea %>%
  dplyr::group_by(GeneSet) %>%
  rstatix::wilcox_test(as.formula(paste("Score ~", group_col))) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p.adj)

write.csv(ssgsea_stats, file = file.path(pipelines, paste0("ssGSEA_stats_by_", group_col, ".csv")), row.names = FALSE)

# FDR filter
sig_sets <- ssgsea_stats %>%
  dplyr::filter(p.adj < 0.1) %>%
  dplyr::pull(GeneSet) %>%
  unique()

# long_ssgsea 
long_sig <- long_ssgsea %>% dplyr::filter(GeneSet %in% sig_sets)

########################################################## HEATMAP ALL ANNOTATIONS ##########################################################

# Expression matrix: samples × gene sets
hm <- as.matrix(combined_data[, gs_cols, drop = FALSE])
rownames(hm) <- combined_data$SampleID

# Sample annotation
annotation_df <- data.frame(
  Model = combined_data$Model,
  row.names = combined_data$SampleID)

# Annotation colors
ann_colors <- list(Model = c("PDX" = "#f5b041", "Human" = "#1abc9c"))

# Significant gene sets
keep_cols <- intersect(colnames(hm), sig_sets)
hm <- hm[, keep_cols, drop = FALSE]

# Order samples by group 
ord <- order(combined_data[[group_col]])
hm_ord <- hm[ord, , drop = FALSE]
annotation_df <- annotation_df[rownames(hm_ord), , drop = FALSE]

# Pathways as rows (scale='row' per pathway)
mat_plot <- t(hm_ord)

breaks <- c(-3, 0, 3)
colors <- colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))

pdf(file.path(output, "ssGSEA_heatmap_all_annotations.pdf"), width = 8, height = 8)
new_heatmap(
  mat = mat_plot,
  scale = "row",
  color = colors,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  annotation_names_col = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  treeheight_row = 0,
  treeheight_col = 25,
  fontsize = 7,
  column_split = factor(annotation_df$Model, levels = c("PDX", "Human")),
  main = paste0(title, " (p.adjust <0.1)"),
  name = "Z-score")

dev.off()

