##########################################################
# maftools compare PDX and Human (PRH)                   #
# JACCARD INDEX                                          #
# PDX VAF < 0.1                                          #
# HUMAN VAF < 0.05                                       #
#                                                        #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 28-10-25                               #
##########################################################

setwd("/Users/mary-kateriley/Library/CloudStorage/OneDrive-Queen'sUniversityBelfast/Data Analysis/HNC/Basal PDX/Reanlysis Exome Basal '24")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(dplyr)
#BiocManager::install("maftools")
library(maftools)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
#BiocManager::install("NMF")
library(NMF)
library(pheatmap)
library(ggplot2)
library(ggtext)

library(RColorBrewer) # For defining the heatmap colors
library(scales) # To show colors if needed

########################################################## READ DATA ##########################################################
# Read in MAF file
pdx_maf <- read.maf("./2_pipelines/MAFtools/PDX/Filtered_HPV_Lymphoma_VAF0.1/Filtered_HN0039A_Filtered_PDX_hardfiltered_maftools.maf")
human_maf <- read.maf("./2_pipelines/MAFtools/Human/Filtered_HPV_Lymphoma/Filtered_HN0039A_Filtered_Human_hardfiltered_maftools.maf")
#ClinicalInfo <- read.csv("./0_data/Human/Human_Complete_Curated_Clinical_Information.csv", header = T)

# Extract Case_ID from first 7 characters
human_maf@clinical.data$Case_ID <- substr(human_maf@clinical.data$Tumor_Sample_Barcode, 1, 7)
pdx_maf@clinical.data$Case_ID   <- substr(pdx_maf@clinical.data$Tumor_Sample_Barcode, 1, 7)

# Extract Case_ID from first 7 characters (ensure this is correct for your IDs)
human_maf@data$Case_ID <- substr(human_maf@data$Tumor_Sample_Barcode, 1, 7)
pdx_maf@data$Case_ID   <- substr(pdx_maf@data$Tumor_Sample_Barcode, 1, 7)

########################################################## SET VARIABLES ##########################################################
pipelines <- "./2_pipelines/MAFtools/PDX_Human_Comparison/Filtered_HPV_Lymphoma_VAF0.1/"
output <- "./3_output/MAFtools/PDX_Human_Comparison/Filtered_HPV_Lymphoma_VAF0.1/"
title <- "HN0039A_Filtered_PDX_Human_Comparison"

########################################################## JACCARD INDEX ##########################################################

## Collect mutated gene sets per Case_ID (union across all samples of a case) 
genes_by_case <- function(maf_obj) {
  cd <- maf_obj@clinical.data
  samples_in_maf <- unique(as.character(maf_obj@data$Tumor_Sample_Barcode)) # Get samples *actually* in the MAF mutation data
  cd_filtered <- cd[cd$Tumor_Sample_Barcode %in% samples_in_maf, ] # Filter clinical data to only those samples
  split(cd_filtered$Tumor_Sample_Barcode, cd_filtered$Case_ID) |> # split the *filtered* list
    lapply(function(tsbs) {
      if (length(tsbs) == 0) return(character(0)) # Safety check
      m <- subsetMaf(maf_obj, tsb = tsbs, mafObj = FALSE)
      unique(m$Hugo_Symbol)
    })
}
human_genes_list <- genes_by_case(human_maf)
pdx_genes_list   <- genes_by_case(pdx_maf)

# Shared Case_IDs present in BOTH cohorts
shared_cases <- intersect(names(human_genes_list), names(pdx_genes_list))
human_genes_list <- human_genes_list[shared_cases]
pdx_genes_list   <- pdx_genes_list[shared_cases]

## Oorder by decreasing diagonal Jaccard so the diagonal looks nicely graded
diag_jaccard <- sapply(shared_cases, function(cid) {
  h <- human_genes_list[[cid]]; p <- pdx_genes_list[[cid]]
  if (length(h) == 0 && length(p) == 0) return(NA_real_)
  length(intersect(h, p)) / length(union(h, p))
})
# Order cases by diagonal similarity (highest -> lowest); fall back to Case_ID if all NAs
ord_cases <- if (all(is.na(diag_jaccard))) shared_cases else shared_cases[order(diag_jaccard, decreasing = TRUE)]

##  Compute all pairwise Jaccards (PDX x Human)
all_pairs <- expand.grid(
  PDX_Case   = ord_cases,
  Human_Case = ord_cases,
  stringsAsFactors = FALSE
)

jacc_fun <- function(pdx_case, human_case) {
  g1 <- pdx_genes_list[[pdx_case]]
  g2 <- human_genes_list[[human_case]]
  if (length(g1) == 0 && length(g2) == 0) return(NA_real_)
  length(intersect(g1, g2)) / length(union(g1, g2))
}

all_pairs$Jaccard <- mapply(jacc_fun, all_pairs$PDX_Case, all_pairs$Human_Case)

## Keep only one triangle (upper, including diagonal)
all_pairs$PDX_idx   <- match(all_pairs$PDX_Case, ord_cases)
all_pairs$Human_idx <- match(all_pairs$Human_Case, ord_cases)

upper_tri <- dplyr::filter(all_pairs, PDX_idx >= Human_idx)

upper_tri$PDX_Case   <- factor(upper_tri$PDX_Case, levels = ord_cases)
upper_tri$Human_Case <- factor(upper_tri$Human_Case, levels = rev(ord_cases))

# Create the cutoff variable 
upper_tri$Cutoff <- ifelse(upper_tri$Jaccard >= 0.6, "≥ 0.6", NA)


# Plot
p_tri <- ggplot(upper_tri, aes(x = PDX_Case, y = Human_Case, fill = Jaccard)) +
  geom_tile(aes(color = Cutoff), size = 0.4) + 
  scale_fill_gradientn(limits = c(0,1),
                       colours = c("white", "#f7fbff", "#4393C3", "#053061"),
                       values  = c(0.00, 0.001, 0.5, 1),
                       oob = scales::squish,
                       name = "Jaccard Index") +
  scale_color_manual(name = "Similarity Cutoff", values = c("≥ 0.6" = "red"), na.value = NA, na.translate = FALSE) +
  coord_equal() +
  labs(title = "PDX vs Human Mutation Overlap (Jaccard Index)",
       x = "CaseID (PDX)",
       y = "CaseID (Human)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        legend.box = "vertical")

ggsave(filename = paste0(output, "Heatmap_PDX_vs_Human_Jaccard_upperTriangle.pdf"),
       plot = p_tri, width = 9, height = 9)

## full matrix (rows: Human, cols: PDX) 
mat <- matrix(NA_real_, nrow = length(ord_cases), ncol = length(ord_cases),
              dimnames = list(Human = ord_cases, PDX = ord_cases))
mat[cbind(all_pairs$Human_idx, all_pairs$PDX_idx)] <- all_pairs$Jaccard
write.csv(mat, file = paste0(pipelines, "PDX_vs_Human_Jaccard_matrix.csv"), row.names = TRUE)

