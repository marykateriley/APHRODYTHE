##########################################################
# SHH Gene Enrichment: GSEA                              #             
# Graft Treatment                                        #
#                                                        #
# Use DESeq result for preranked GSEA for gene sets:     #
# visit MsigDB for more genesets:                        #
# http://www.gsea-msigdb.org/gsea/index.jsp              #
#                                                        #
# Non-Responders v Responders                            #
# Avg 3-6 week response                                  #
# Technical, lymphoma + HPV samples removed              #
#                                                        #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 05-11-24                               #
##########################################################


setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(GSVA)
library(GSEABase)
library(pheatmap)

set.seed(1234)

########################################################## READ DATA ##########################################################
# Read in SHH GMT File
SHH <- read.csv("./0_data/GMT/SHH/Curated_SHH_GMT.csv", row.names = 1)

########################################################## READ DATA FOR GSEA ##########################################################
# Read DESeq results
DESeq_results <- read.csv("./2_pipelines/DESeq2/Filtered/Response/NRvR/DESeq2_Definitive_Response_3_6wkNonResponder_v_Responder.csv", header = T, row.names = 1) 
DESeq_results <- DESeq_results[!is.na(DESeq_results$padj),] # remove NA values for padj

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/Enrichment_Analysis/SHH/Filtered/NRvR/"
output <- "./3_output/Enrichment_Analysis/SHH/Filtered/NRvR/"
name <- "SHH_NRvR"
title <- "Basal HNSCC: Average 3-6 week NonResponders v Responders"

dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

########################################################## CURATE GSEA GENE LIST ##########################################################

######### STEP 1: create geneList (based on stat value from DESeq resuls)
## Feature 1: numeric vector(save _stat_ column from DESeq results)
geneList = DESeq_results[,4]

## feature 2: named vector - adding the genenames from DESeq results
names(geneList) = as.character(rownames(DESeq_results))

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)


########################################################## GSEA ##########################################################


gseaResult <-GSEA(geneList, exponent = 1, pvalueCutoff = 1,
                  pAdjustMethod = "BH", TERM2GENE = SHH, verbose = TRUE, seed = TRUE) #can set seed=True or just run once and save result for further plots! 
gseaRes_df <- as.data.frame(gseaResult)

write.csv(gseaRes_df, file=paste0(pipelines,name,"_GSEA.csv"))


#plot all Hallmarks pathways in order of enrichment, in BARPLOT
pdf(file = paste0(output,name,"_GSEA.pdf"), width = 14, height = 11)
ggplot(gseaRes_df, aes(x = reorder(Description, NES), y = NES, fill = p.adjust < 0.1)) +
  geom_col() +
  coord_flip() +  # Flip coordinates to create horizontal bars
  scale_fill_manual(values = c("TRUE" = "#29BF4E", "FALSE" = "#EBF1EC")) +  # Explicitly map TRUE and FALSE
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0(title),
       subtitle = "Enrichment of SHH Gene Set: ",
       labels = c("TRUE", "FALSE")) + 
  theme(axis.line = element_line(),
        plot.title = element_text(size = 20, face="bold"),
        plot.subtitle = element_text(size = 18),
        axis.title.y = element_text(size = 18,  colour="black", face="bold"),
        axis.text.y = element_text(size=16,colour="black"),
        axis.title.x = element_text(size = 18,  colour="black", face="bold"),
        axis.text.x = element_text(size=12,colour="black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank()) 
dev.off()


