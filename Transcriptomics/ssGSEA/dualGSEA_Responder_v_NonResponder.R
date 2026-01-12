##########################################################
# dualGSEA                                               #
# https://github.com/MolecularPathologyLab/Bull-et-al    #
# Responders v Non-Responders                            #
# Avg 3-6 week response                                  #
# Uses limma for DGEA                                    #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 31-03-2025                             #
##########################################################

setwd("/Users/mary-kateriley/Library/CloudStorage/OneDrive-Queen'sUniversityBelfast/Data Analysis/HNC/Basal PDX/Reanalysis Basal '24")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

## List of packages
bio_pkgs <- c("limma", "msigdbr", "DOSE", "devtools","fgsea", "enrichplot","GSVA","circlize",
              "ComplexHeatmap","ggridges","pROC","cutpointr","ggalluvial","waterfalls",
              "randomForest","devtools","tibble","msigdbr","data.table","tidyverse","dplyr",
              "ggplot2","grid","gridExtra","ROCR","reshape2","stats","RColorBrewer","ggpubr",
              "ggbeeswarm","tidyr","caret","escape", "doRNG")
github_pkgs <- c("kassambara/easyGgplot2")


## If not already installed
#suppressMessages({
#  if (!require(bio_pkgs, quietly = TRUE))
#    BiocManager::install(bio_pkgs)
#  if (!require(github_pkgs, quietly = TRUE))
#    devtools::install_github(github_pkgs, dependencies = TRUE, force=TRUE)
#  })


## Load R packages
suppressPackageStartupMessages(suppressMessages(
  lapply(bio_pkgs, library, character.only = TRUE)
))
suppressPackageStartupMessages(suppressMessages(
  library(easyGgplot2)
))

library(clusterProfiler)
library(cowplot) 
library(DESeq2)
########################################################## READ DATA ##########################################################

# Recommendation - normalised counts
data <- read.csv("./0_data/DESeq2/PRX/GEX/All_RemovedSamples_Unfiltered_RawCounts_GraftHNSCC_270424.csv", row.names = 1)
group_data <- read.csv("./0_data/DESeq2/PRX/Metadata/All_RemovedSamples_Unfiltered_GraftHNSCC_colData_240724.csv", header = T, row.names = 1)

# Read in SHH GMT and EGFR ligands Files
SHH <- read.gmt("./0_data/GMT/SHH/BIOCARTA_SHH_PATHWAY.v2023.2.Hs.gmt")
EGFR <- read.csv("./0_data/GMT/EGFR/EGFR_ligands.csv", header = T)


## Geneset list - we will use MSigDB to extract Hallmark collection for this example.
geneset1 <- msigdbr::msigdbr(category = "H")
geneset_list <- geneset1 %>% split(f = .$gs_name, x = .$gene_symbol)


########################################################## DATA CURATION ##########################################################

# Remove any sample that does not have Definitive_Response_3wk (NAs) from group_data
group_data <- group_data  %>%
  filter(!is.na(Definitive_Response_3_6wk))

# Ensure all samples in data match those samples that have a 3-week response in group_data
data <- data[, colnames(data) %in% rownames(group_data)]

# Set row names of group_data to match the samples in data
rownames(group_data) <- colnames(data)

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/Enrichment_Analysis/dualGSEA/Filtered/Response/RvNR"
output <- "./3_output/Enrichment_Analysis/dualGSEA/Filtered/Response/RvNR"
title <- "HNSCC: "
subtitle <- "Basal HNSCC PDX: Filtered Technical; Lymphoma; HPV"


########################################################## dualGSEA  ##########################################################
##### GEX ##### 
dds <- DESeqDataSetFromMatrix(countData = data, 
                              colData = group_data, 
                              design = ~ 1)
ddsBlind <- DESeq(dds)


# Remove genes which have less than 10 counts
rows_less_than_10 <- which(rowSums(counts(ddsBlind)) < 10) # 21108 genes
keep <- rowSums(counts(ddsBlind)) >= 10 
ddsBlind <- ddsBlind[keep,]

# Normalised count matrix
data <- counts(ddsBlind)
data <- as.data.frame(data)

## Get dualGSEA() function
dualgsea_fun_file <- here::here("1_code", "Enrichment_Analysis", "dualGSEA", "dualGSEA.R")
source(dualgsea_fun_file)

## Run dualGSEA 
## NOTE: invisible() and capture.output() is being used here to hide messages and progressbar
## that would be shown in the document otherwise.
set.seed(127)
invisible(capture.output(
  
  res1 <- dualGSEA(data = data, 
                   group_data = group_data, 
                   group_colname = "Definitive_Response_3_6wk",
                   geneset_list = geneset_list)
  
))

## Outputs in a list
names(res1)



## 1. Differential Gene Ranking ----
print(head(res1[["Differential_GeneRanking"]])[1:5,])

## 2. Pairwise GSEA Result Table ----
print(head(res1[["Pairwise_ResultTable"]])[1:5, 1:6])

## 3. Pairwise GSEA Barplot ----
print(res1[["Pairwise_BarPlot"]])
ggsave(filename = file.path(output, "pairwise_barplot.png"), 
       plot = res1[["Pairwise_BarPlot"]], 
       width = 11, height = 7, dpi = 300) 

## 4. Pairwise Enrichments ----
print(res1[["Pairwise_EnrichmentPlots"]])

output_pdf <- file.path(output, "pairwise_EnrichmentPlots_combined.pdf")
pdf(output_pdf, width = 11, height = 7)
for (i in 1:length(res1[["Pairwise_EnrichmentPlots"]])) {
  print(res1[["Pairwise_EnrichmentPlots"]][[i]])
}
dev.off()


## 5. Single Sample GSEA Result Table ----
print(head(res1[["SingleSample_ResultTable"]])[1:5, 1:6])

## 6. Single Sample GSEA Density Plots ----
print(res1[["SingleSample_DensityPlots"]])
output_pdf <- file.path(output, "SingleSample_DensityPlots_combined.pdf")
pdf(output_pdf, width = 11, height = 7)
for (i in 1:length(res1[["SingleSample_DensityPlots"]])) {
  print(res1[["SingleSample_DensityPlots"]][[i]])
}
dev.off()


## 7. Single Sample GSEA Histogram ----
print(res1[["SingleSample_Histograms"]])
output_pdf <- file.path(output, "SingleSample_Histograms_combined.pdf")
pdf(output_pdf, width = 11, height = 7)
for (i in 1:length(res1[["SingleSample_Histograms"]])) {
  print(res1[["SingleSample_Histograms"]][[i]])
}
dev.off()

## 8. Single Sample GSEA ROC Plot ----
print(res1[["SingleSample_ROCPlots"]])
output_pdf <- file.path(output, "SingleSample_ROCPlots_combined.pdf")
pdf(output_pdf, width = 11, height = 7)
for (i in 1:length(res1[["SingleSample_ROCPlots"]])) {
  print(res1[["SingleSample_ROCPlots"]][[i]])
}
dev.off()

## 9. Single Sample GSEA Optimum Cutoff Labels ----
print(res1[["SingleSample_OptimumCutoff_Labels"]])

## 10. Single Sample Waterfall Plot ----
print(res1[["SingleSample_WaterfallPlots"]])
output_pdf <- file.path(output, "SingleSample_WaterfallPlots_combined.pdf")
pdf(output_pdf, width = 11, height = 7)
for (i in 1:length(res1[["SingleSample_WaterfallPlots"]])) {
  print(res1[["SingleSample_WaterfallPlots"]][[i]])
}
dev.off()

# ---- SESSION INFO ----
devtools::session_info()


