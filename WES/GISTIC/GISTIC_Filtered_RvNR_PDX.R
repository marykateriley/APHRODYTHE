##########################################################
# GISTIC Copy Number PDX                                 #
# RESPONDERS + NONRESPONDERS                             #
# Filtered out HPV/Lymphoma                              #
#                                                        #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 23-08-2025                             #
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
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
#BiocManager::install("NMF")
library(NMF)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(gplots)
library(clusterProfiler)
library(RColorBrewer)


########################################################## READ DATA ##########################################################

# Read in metadata file
ClinicalInfo <- read.csv("./0_data/PDX/PDX_Rem_HPV_Lymphoma_Clinical_Information.csv", header = T, row.names = 1)
ClinicalInfo <- ClinicalInfo  %>%
  filter(!is.na(Definitive_Response_3_6wk))

ClinicalInfo$Definitive_Response_3_6wk <- factor(ClinicalInfo$Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))
levels(ClinicalInfo$Definitive_Response_3_6wk)

ClinicalInfo$Response_3_6wk <- factor(ClinicalInfo$Response_3_6wk, levels = c("PR", "SD", "PD", "No_data"))
levels(ClinicalInfo$Response_3_6wk)

# Set rownames in ClinicalInfo to Tumor_Sample_Barcode
head(ClinicalInfo)
rownames(ClinicalInfo)=ClinicalInfo[, ncol(ClinicalInfo)] # Set last column as rownames (Tumor_Sample_Barcode)
head(ClinicalInfo)

# Read data
all_tresh_genes <- read.delim("./0_data/GISTIC/PDX_HPVLymp_Rem/all_thresholded.by_genes.txt", header = TRUE, stringsAsFactors = FALSE)


# Read in SHH GMT and EGFR ligands Files and KRT genes of interest
EGFR <- read.csv("./0_data/GMT/EGFR/EGFR_ligands.csv", header = T)
EGFR <- EGFR$EGF_Ligands
EGFR <- c("EGFR",unique(EGFR) %>% 
            setdiff("EGFR") %>% 
            sort())
EGFR <- c(EGFR, "CAV1", "SOX2")
reccurant_cna <- c("TP63", "SOX2", "PIK3CA", "CCDN1", "FAT1", "CSMD1", "CDKN2A", "SMAD4", "CCND1")
########################################################## SET VARIABLES ##########################################################
data <- "./0_data/PDX/"
pipelines <- "./2_pipelines/GISTIC/PDX/NEW_Filtered_HPV_Lymphoma/"
output <- "./3_output/GISTIC/PDX/NEW_Filtered_HPV_Lymphoma/"
title <- "HN0039A_Filtered_PDX_RVNR_"
plot_title <- "GISTIC PDX:"

########################################################## ALL THRESHOLD BY GENES HEATMAP ##########################################################


# Make longer
merged_data <- all_tresh_genes %>%
  pivot_longer(cols = starts_with("log2"),
               names_to = "Sample_ID",
               values_to = "CopyNumber")
merged_data$Tumor_Sample_Barcode <- substr(merged_data$Sample_ID, 6, 17)

# Create list with all genes of interest
all_genes_of_interest <- list(EGFR, reccurant_cna)
all_genes_of_interest <- unlist(all_genes_of_interest)

# Filter the merged data to only include the genes of interest
long_gene_data <- merged_data %>%
  filter(Gene.Symbol %in% all_genes_of_interest) 
colnames(long_gene_data)

# Make longer so genes are rownames
long_gene_data <- long_gene_data %>% 
  dplyr::select(Tumor_Sample_Barcode, Gene.Symbol, CopyNumber) %>% 
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = CopyNumber, values_fill = 0) %>% 
  column_to_rownames("Gene.Symbol")

# Subset annotation data to include only samples present in gistic_output
annotations <- ClinicalInfo[colnames(long_gene_data), c("Definitive_Response_3_6wk", "Response_3_6wk", "Sex", "Site.of.Primary")]
annotations <- annotations[!is.na(annotations$Definitive_Response_3_6wk), ]
annotations$Site.of.Primary[annotations$Site.of.Primary == ""] <- "Unknown"

annotations$Site.of.Primary <- factor(annotations$Site.of.Primary)
annotations$Site.of.Primary

# Set column names of long_gene_data to match row names
long_gene_data <- long_gene_data[, rownames(annotations), drop = FALSE]

# Set column names of ind to match row names 
identical(rownames(annotations), colnames(long_gene_data)) # check if row names match

# Set annotation colours
ann_colors <- list(
  Definitive_Response_3_6wk = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF", "No_data" = "#666666FF"),
  Response_3_6wk = c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF", "No_data" = "#666666FF"),
  Sex = c("M" = "#A6CEE3", "F" = "#FB9A99"),
  Site.of.Primary = c("Larynx" = "#1B9E77FF", "Oropharynx" = "#66A61EFF", "Hypopharynx" = "#7570B3FF", "Oral Cavity" = "#E6AB02FF", "Larynx + Oral Cavity" = "#D95F02FF", "Oral Cavity + Oropharynx" = "#E7298AFF", "Unknown" = "#666666FF"),
  CNV = c("-2" = "blue", "-1" = "lightblue", "0" = "white", "1" = "pink","2" = "red"))

# Discrete colour mapping + labels
CNV_colors <- circlize::colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("blue", "lightblue", "white", "pink", "red"))

legend_breaks <- c(-2, -1, 0, 1, 2)
legend_labels <- c("Del (-2)", "Loss (-1)", "Stable (0)", "Gain (+1)", "Amp (+2)")


## EGFR ##
EGFR_long <- rownames(long_gene_data) %in% EGFR
EGFR_long <- long_gene_data[EGFR_long, ]
EGFR_long <- as.matrix(EGFR_long)
EGFR_long <- EGFR_long[EGFR[EGFR %in% rownames(EGFR_long)], ]

# Unique CNV values present in the data
present_values <- sort(unique(as.numeric(EGFR_long)))

pdf(paste0(output, title, "Copy_Number_Alterations_EGFR_Resp.pdf"), width = 10, height = 10)
ComplexHeatmap::pheatmap(
  EGFR_long,
  color = CNV_colors,
  cluster_cols = TRUE,  # Cluster samples
  cluster_rows = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotations,
  annotation_colors = ann_colors,  
  column_split = as.factor(annotations$Definitive_Response_3_6wk),
  annotation_names_col = FALSE,
  fontsize_row = 12,  # Adjust text size
  main = paste0(plot_title, " Copy Number Alterations: EGFR and Ligands"),
  name = "CNV",
  heatmap_legend_param = list(
    at = c(-2,-1,0,1,2),
    labels = c("Del (-2)","Loss (-1)","Stable (0)","Gain (+1)","Amp (+2)"),
    color_bar = "discrete"))
dev.off()


## RECCURANT CNAs ##
CNA_long <- rownames(long_gene_data) %in% reccurant_cna 
CNA_long <- long_gene_data[CNA_long, ]
CNA_long <- as.matrix(CNA_long)

# Get the unique CNV values present in the data
present_values <- sort(unique(as.numeric(CNA_long)))

pdf(paste0(output, title, "Copy_Number_Alterations_ReccurantCNA_Resp.pdf"), width = 10, height = 10)
ComplexHeatmap::pheatmap(
  CNA_long,
  color = CNV_colors,
  cluster_cols = TRUE,  # Cluster samples
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotations,
  annotation_colors = ann_colors,  
  column_split = as.factor(annotations$Definitive_Response_3_6wk),
  annotation_names_col = FALSE,
  fontsize_row = 12,  # Adjust text size
  main = paste0(plot_title, " Copy Number Alterations: Reccurant CNAs"),
  name = "CNV",
  heatmap_legend_param = list(
    at = c(-2,-1,0,1,2),
    labels = c("Del (-2)","Loss (-1)","Stable (0)","Gain (+1)","Amp (+2)"),
    color_bar = "discrete"))
dev.off()

