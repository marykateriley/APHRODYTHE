##########################################################
# ssGSEA (GSVA package) TPM                              #
# Graft                                                  #
# Responder v NonResponder                               #
# GO KRT                                                 #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 05-11-2025                             #
##########################################################

setwd("")

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


col <- msigdbr_collections()

########################################################## READ DATA ##########################################################
# Load filtered Graft Raw Counts Matrix - from DESeq_Data_Input.R
# As ssGSEA requires integers 
# Read in GEX for GSVA
GEX <- read.csv("./0_data/edgeR/PRX/TPM/PDX_TPM_GEX.csv", row.names = 1)
rownames(GEX) <- toupper(rownames(GEX))
GEX <- log2(GEX + 1)

metadata <- read.csv("./0_data/DESeq2/PRX/Metadata/All_RemovedSamples_Unfiltered_GraftHNSCC_colData_240724.csv", header = T, row.names = 1)

G_t2g <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>% 
  dplyr::select(gs_name, gene_symbol)

network <- G_t2g %>%
  filter(!is.na(gs_name), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
  distinct(gs_name, gene_symbol) %>%
  rename(GeneSet = gs_name, Gene = gene_symbol)

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/Enrichment_Analysis/ssGSEA/Filtered/TPM_GOBP/"
output <- "./3_output/Enrichment_Analysis/ssGSEA/Filtered/TPM_GOBP/"
name <- "ssGSEA"
title <- "PDX HNSCC: ssGSEA scores GO BP"

dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

source("./1_code/RFunctions/new_heatmap.R")

########################################################## DATA CURATION ##########################################################

# Remove samples with no Definitive_Response_3_6wk data
na_samples <- metadata$Sample_ID[is.na(metadata$Definitive_Response_3_6wk)]
GEX <- GEX[, !colnames(GEX) %in% na_samples]
metadata <- metadata %>% 
  filter(!Sample_ID %in% na_samples)

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
write.csv(ssgsea_1, file = paste0(pipelines,"ssGSEAScores_Graft.csv"))

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
ssgsea_1$SampleID <- NULL

########################################################## STATS ##########################################################

# Group
group_col <- "Definitive_Response_3_6wk"

# Identify ssGSEA columns
gs_cols <- setdiff(colnames(combined_data), c(colnames(metadata), "SampleID"))

# Long format
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

# Only filder for keratin / cornified pathways
ssgsea_stats_krt <- ssgsea_stats %>%
  dplyr::filter(grepl("keratin|cornified", GeneSet, ignore.case = TRUE))

write.csv(ssgsea_stats, file = file.path(pipelines, paste0("ssGSEA_stats_by_", group_col, ".csv")), row.names = FALSE)

write.csv(ssgsea_stats_krt, file = file.path(pipelines, paste0("ssGSEA_krt_stats_by_", group_col, ".csv")), row.names = FALSE)

# Basic FDR filter
sig_sets_krt <- ssgsea_stats_krt %>%
  dplyr::filter(p.adj < 0.05) %>%
  dplyr::arrange(p.adj) %>%
  dplyr::pull(GeneSet) %>%
  unique()

# long_ssgsea came from your earlier code
long_sig_krt <- long_ssgsea %>% dplyr::filter(GeneSet %in% sig_sets_krt)


########################################################## HEATMAP KRT ANNOTATIONS ##########################################################

# Expression matrix: samples × gene sets
hm <- as.matrix(combined_data[, gs_cols, drop = FALSE])
rownames(hm) <- combined_data$SampleID

# Sample annotation 
annotation_df <- data.frame(
  Definitive_Response_3_6wk = combined_data$Definitive_Response_3_6wk,
  Response_3_6wk = combined_data$Response_3_6wk,
  percentage_vol_change_3_6wk = combined_data$percentage_vol_change_3_6wk,
  Sex = combined_data$Sex,
  Site.of.Primary = combined_data$Site.of.Primary,
  row.names = combined_data$SampleID,
  check.names = FALSE)


# continuous mapper for % change
pal <- rev(RColorBrewer::brewer.pal(11, "PRGn"))
pct_col_fun <- circlize::colorRamp2(c(-90, 0, 335), c(pal[2], pal[6], pal[10]))

#  numeric
annotation_df$percentage_vol_change_3_6wk <- as.numeric(as.character(annotation_df$percentage_vol_change_3_6wk))

# Annotation colors
ann_colors <- list(
  Model = c("PDX" = "#f5b041", "Human" = "#1abc9c"),
  Definitive_Response_3_6wk = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
  Response_3_6wk = c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF"),
  percentage_vol_change_3_6wk = pct_col_fun,
  Sex = c("M" = "#A6CEE3", "F" = "#FB9A99"),
  Site.of.Primary = c("Larynx" = "#1B9E77FF", "Oropharynx" = "#66A61EFF", "Hypopharynx" = "#7570B3FF", "Oral Cavity" = "#E6AB02FF", "Larynx + Oral Cavity" = "#D95F02FF", "Oral Cavity + Oropharynx" = "#E7298AFF", "Unknown" = "#666666FF"))



# significant gene sets
keep_cols <- intersect(colnames(hm), sig_sets_krt)
hm <- hm[, keep_cols, drop = FALSE]

# order samples by group
ord <- order(combined_data[[group_col]])
hm_ord <- hm[ord, , drop = FALSE]
annotation_df <- annotation_df[rownames(hm_ord), , drop = FALSE]

mat_plot <- t(hm_ord) 

breaks <- c(-3, 0, 3)
colors <- colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))

pdf(file.path(output, "ssGSEA_heatmap_krt_annotations.pdf"), width = 9, height = 9)
ht <- new_heatmap(
  mat = mat_plot,
  scale = "row",
  color = colors,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  annotation_names_col = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  treeheight_row = 25,
  treeheight_col = 25,
  fontsize = 8,
  column_split = factor(annotation_df$Definitive_Response_3_6wk, levels = c("NonResponder", "Responder")),
  main = paste0(title, " Keratin GeneSets (p.adjust < 0.05)"),
  name = "Z-score",
  row_names_max_width = grid::unit(12, "cm"),
  annotation_show_legend = c(percentage_vol_change_3_6wk = TRUE),
  legend_orientation = "horizontal",
  annotation_legend_param = list(
    percentage_vol_change_3_6wk = list(
      title = "% volume change 3–6 wk",
      at = c(-90, 0, 100, 200, 335),
      labels = c("-90", "0", "100", "200", "335"),
      direction = "horizontal")))

ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

