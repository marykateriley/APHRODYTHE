##########################################################
# GISTIC Copy Number                                     #
# HUMAN + PDX                                            #
# Filtered out HPV/Lymphoma                              #
#                                                        #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 28-10-2025                             #
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

########################################################## READ DATA ##########################################################

# Read in metadata file
ClinicalInfo_PDX <- read.csv("./0_data/PDX/PDX_Rem_HPV_Lymphoma_Clinical_Information.csv", header = T, row.names = 1)

# Set rownames in ClinicalInfo_PDX to Tumor_Sample_Barcode
head(ClinicalInfo_PDX)
rownames(ClinicalInfo_PDX)=ClinicalInfo_PDX[, ncol(ClinicalInfo_PDX)] # Set last column as rownames (Tumor_Sample_Barcode)
head(ClinicalInfo_PDX)


# Read in metadata file
ClinicalInfo_Human <- read.csv("./0_data/Human/Human_Rem_HPV_Lymphoma_Clinical_Information.csv", header = T, row.names = 1)

# Set rownames in ClinicalInfo_Human to Tumor_Sample_Barcode
head(ClinicalInfo_Human)
rownames(ClinicalInfo_Human)=ClinicalInfo_Human[, ncol(ClinicalInfo_Human)] # Set last column as rownames (Tumor_Sample_Barcode)
head(ClinicalInfo_Human)

# log2 PDX vs log2 NMH
gistic_output_pdx = read.table("./0_data/GISTIC/PDX_HPVLymp_Rem/broad_values_by_arm.txt", header=T, row.names = 1, sep="\t")
colnames(gistic_output_pdx)

# log2 PDX vs log2 NMH
gistic_output_human = read.table("./0_data/GISTIC/Human_HPVLymp_Rem/broad_values_by_arm.txt", header=T, row.names = 1, sep="\t")
colnames(gistic_output_human)


########################################################## SET VARIABLES ##########################################################
data <- "./0_data/"
pipelines <- "./2_pipelines/GISTIC/PDXHuman/NEW_Filtered_HPV_Lymphoma/"
output <- "./3_output/GISTIC/PDXHuman/"
title <- "HN0039A_Filtered_"
plot_title <- "GISTIC:"

########################################################## BROAD VALUEWS BY ARM PDX ##########################################################

# Rename column names so it matches the PDX sample names
colnames(gistic_output_pdx) <- substr(colnames(gistic_output_pdx), 6, 17)
colnames(gistic_output_pdx)

# Only include columns that have samples present in ClinicalInfo_PDX
gistic_output_pdx <- gistic_output_pdx[, colnames(gistic_output_pdx) %in% rownames(ClinicalInfo_PDX)]
colnames(gistic_output_pdx)

gistic_output_pdx <- as.matrix(gistic_output_pdx)

# Subset annotation data to include only samples present in gistic_output_pdx
annotations <- ClinicalInfo_PDX[colnames(gistic_output_pdx), c("Definitive_Response_3_6wk", "Response_3_6wk", "Sex", "Site.of.Primary")]
annotations <- data.frame(annotations)
annotations$Site.of.Primary[annotations$Site.of.Primary == ""] <- "Unknown"

annotations$Site.of.Primary <- factor(annotations$Site.of.Primary)
annotations$Site.of.Primary

# Set column names of ind to match row names 
colnames(gistic_output_pdx) <- rownames(annotations)
identical(rownames(annotations), colnames(gistic_output_pdx)) # check if row names match



# Set annotation colours
ann_colors <- list(
  Model = c("PDX" = "#f5b041", "Human" = "#1abc9c"),
  Definitive_Response_3_6wk = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF", "No_data" = "#666666FF"),
  Response_3_6wk = c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF", "No_data" = "#666666FF"),
  Sex = c("M" = "#A6CEE3", "F" = "#FB9A99"),
  Site.of.Primary = c("Larynx" = "#1B9E77FF", "Oropharynx" = "#66A61EFF", "Hypopharynx" = "#7570B3FF", "Oral Cavity" = "#E6AB02FF", "Larynx + Oral Cavity" = "#D95F02FF", "Oral Cavity + Oropharynx" = "#E7298AFF", "Unknown" = "#666666FF")
)
# Count values in 0.1 bins (incl. negatives)
g_vals <- as.numeric(na.omit(as.vector(gistic_output_pdx)))
bins <- cut(g_vals, breaks = seq(floor(min(g_vals)), ceiling(max(g_vals)), by = 0.1), right = FALSE)
table(bins)

# Find  min and max
min_val <- min(gistic_output_pdx, na.rm = TRUE) # -1.293
max_val <- max(gistic_output_pdx, na.rm = TRUE) # +3.657

# Define the number of color steps for negative and positive sides
n_colors_neg <- 50
n_colors_pos <- 50
total_colors <- n_colors_neg + n_colors_pos # 100

# Create breaks for the negative side (min_val to 0)
breaks_neg <- seq(min_val, 0, length.out = n_colors_neg + 1)

# Create breaks for the positive side (0 to max_val)
breaks_pos <- seq(0, max_val, length.out = n_colors_pos + 1)
# Combine, removing the duplicate 0
heatmap_breaks_asym <- unique(c(breaks_neg, breaks_pos))

# Generate the diverging color palette with the total number of colors
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(total_colors)

# Breaks centered at zero
pdf(paste0(output, title, "PDX_BroadValuesByArm_ThesisReady_AsymBreaks.pdf"), width = 10, height = 9)
ComplexHeatmap::pheatmap(
  mat = gistic_output_pdx,
  color = heatmap_colors,
  breaks = heatmap_breaks_asym,
  na_col = "grey",
  legend_breaks = c(round(min_val, 1), 0, round(max_val, 1)),
  annotation_col = annotations,
  annotation_colors = ann_colors,
  annotation_names_col = FALSE,
  cluster_cols = TRUE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 12,
  name = " ",
  main = paste0("PDX Broad Values by Arm"))
dev.off()


rm(gistic_output_pdx)
rm(ClinicalInfo_PDX)


########################################################## BROAD VALUEWS BY ARM HUMAN ##########################################################

# Rename column names so it matches the PDX sample names
colnames(gistic_output_human) <- substr(colnames(gistic_output_human), 6, 17)
colnames(gistic_output_human)

# Only include columns that have samples present in ClinicalInfo_Human
gistic_output_human <- gistic_output_human[, colnames(gistic_output_human) %in% rownames(ClinicalInfo_Human)]
colnames(gistic_output_human)

gistic_output_human <- as.matrix(gistic_output_human)

# Subset annotation data to include only samples present in gistic_output_human
annotations <- ClinicalInfo_Human[colnames(gistic_output_human), c("Definitive_Response_3_6wk", "Response_3_6wk", "Sex", "Site.of.Primary")]
annotations <- data.frame(annotations)
annotations$Site.of.Primary[annotations$Site.of.Primary == ""] <- "Unknown"

annotations$Site.of.Primary <- factor(annotations$Site.of.Primary)
annotations$Site.of.Primary

# Set column names of ind to match row names 
colnames(gistic_output_human) <- rownames(annotations)
identical(rownames(annotations), colnames(gistic_output_human)) # check if row names match

# Set annotation colours
ann_colors <- list(
  Model = c("PDX" = "#f5b041", "Human" = "#1abc9c"),
  Definitive_Response_3_6wk = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF", "No_data" = "#666666FF"),
  Response_3_6wk = c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF", "No_data" = "#666666FF"),
  Sex = c("M" = "#A6CEE3", "F" = "#FB9A99"),
  Site.of.Primary = c("Larynx" = "#1B9E77FF", "Oropharynx" = "#66A61EFF", "Hypopharynx" = "#7570B3FF", "Oral Cavity" = "#E6AB02FF", "Larynx + Oral Cavity" = "#D95F02FF", "Oral Cavity + Oropharynx" = "#E7298AFF", "Unknown" = "#666666FF"))


# Count values in 0.1 bins (incl. negatives)
g_vals <- as.numeric(na.omit(as.vector(gistic_output_human)))
bins <- cut(g_vals, breaks = seq(floor(min(g_vals)), ceiling(max(g_vals)), by = 0.1), right = FALSE)
table(bins)

# Find min and max
min_val <- min(gistic_output_human, na.rm = TRUE) # -1.293
max_val <- max(gistic_output_human, na.rm = TRUE) # +3.657

# Define the number of color steps for negative and positive sides
n_colors_neg <- 50
n_colors_pos <- 50
total_colors <- n_colors_neg + n_colors_pos # 100

# Create breaks for the negative side (min_val to 0)
breaks_neg <- seq(min_val, 0, length.out = n_colors_neg + 1)

# Create breaks for the positive side (0 to max_val)
breaks_pos <- seq(0, max_val, length.out = n_colors_pos + 1)

# Combine, removing the duplicate 0
heatmap_breaks_asym <- unique(c(breaks_neg, breaks_pos))

# Generate the diverging color palette with the total number of colors
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(total_colors)

# Breaks centered at zero
pdf(paste0(output, title, "Human_BroadValuesByArm_ThesisReady_AsymBreaks.pdf"), width = 10, height = 9)
ComplexHeatmap::pheatmap(
  mat = gistic_output_human,
  color = heatmap_colors,
  breaks = heatmap_breaks_asym,
  na_col = "grey",
  legend_breaks = c(round(min_val, 1), 0, round(max_val, 1)),
  annotation_col = annotations,
  annotation_colors = ann_colors,
  annotation_names_col = FALSE,
  cluster_cols = TRUE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 12,
  name = " ",
  main = paste0("Human Broad Values by Arm"))
dev.off()