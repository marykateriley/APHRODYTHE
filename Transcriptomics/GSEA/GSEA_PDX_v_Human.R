##########################################################
# GSEA: Hallmarks, GO, KEGG gene sets                    #             
# Graft Treatment                                        #
#                                                        #
# Use DESeq result for preranked GSEA for gene sets:     #
# visit MsigDB for more genesets:                        #
# http://www.gsea-msigdb.org/gsea/index.jsp              #
#                                                        #
# PDX v Humans                                           #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 07-05-25                               #
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

set.seed(1234)

########################################################## READ DATA ##########################################################

# Read DESeq results
DESeq_results <- read.csv("./2_pipelines/PDX_Human/DESeq2/Unfiltered/DESeq2_ModelPDX_v_Human.csv", header = T, row.names = 1) 
DESeq_results <- DESeq_results[!is.na(DESeq_results$padj),] # remove NA values for padj

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/PDX_Human/Enrichment_Analysis/GSEA/Unfiltered/"
output <- "./3_output/PDX_Human/Enrichment_Analysis/GSEA/Unfiltered/"
name <- "GSEA_PDXvHuman"
title <- "Basal HNSCC PDX v Human: Filtered Technical; Lymphoma; HPV"

########################################################## GSEA using Cluster Profiler ##########################################################

######### STEP 1: create geneList (based on stat value from DESeq resuls)
## Feature 1: numeric vector(save _stat_ column from DESeq results)
geneList = DESeq_results[,4]

## feature 2: named vector - adding the genenames from DESeq results
names(geneList) = as.character(rownames(DESeq_results))

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
view(geneList)

############### STEP 2: Gather HALLMARKS Genesets for GSEA analysis from msigdb
# Format TERM2GENE df with gene names and signature 
# Change category to change which gene sets you analyse

H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# Gather GO gene sets
G_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)


############### STEP 3: RUN GSEA USING VALUES FROM STEP 1 & 2

gseaResult <-GSEA(geneList, exponent = 1, pvalueCutoff = 1,
                  pAdjustMethod = "BH", TERM2GENE = H_t2g, verbose = TRUE, seed = TRUE) #can set seed=True or just run once and save result for further plots! 

gseaRes_df <- as.data.frame(gseaResult) #save results as df

#save results as csv file: Make sure to include the comparison for groups (from DESeq) in your results file name
write.csv(gseaRes_df, file=paste0(pipelines,name,"_Hallmarks.csv"))


###############
#STEP 4: plot overall results as a barplot to see which pathways are most up/down in your analysis 
###############

#save as pdf with specific file name
pdf(file = paste0(output,name,"_Hallmarks.pdf"), width = 14.5, height = 11)

#plot all Hallmarks pathways in order of enrichment, in BARPLOT
ggplot(gseaRes_df, aes(x = reorder(Description, NES), y = NES, fill = p.adjust < 0.1)) +
  geom_col() +
  coord_flip() +  # Flip coordinates to create horizontal bars
  scale_fill_manual(values = c("TRUE" = "#29BF4E", "FALSE" = "#EBF1EC")) +  # Explicitly map TRUE and FALSE
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0(title),
       subtitle = "Enrichment of The Hallmarks: ",
       labels = c("TRUE", "FALSE")) +     #change title to match your analysis
  theme(axis.line = element_line(),
        plot.title = element_text(size = 20, face="bold"),
        plot.subtitle = element_text(size = 18),
        axis.title.y = element_text(size = 18,  colour="black", face="bold"),
        axis.text.y = element_text(size=12,colour="black"),
        axis.title.x = element_text(size = 18,  colour="black", face="bold"),
        axis.text.x = element_text(size=12,colour="black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank()) 

dev.off()



############### GO ###############


gse <- gseGO(geneList = geneList, 
             ont = "ALL", 
             keyType = "SYMBOL", 
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")
write.csv(gse, file=paste0(pipelines,name,"_GO.csv"))

plot <- enrichplot::dotplot(gse, 
                            showCategory = 10, 
                            title = NULL,
                            split = ".sign",
                            font.size = 15) + 
  ggtitle(paste0(title)) + 
  labs(subtitle = "Enriched GO Pathways") +
  facet_grid(.~.sign)  +
  theme(plot.title = element_text(size = 22, face = "bold"),
        plot.subtitle = element_text(size = 20),
        strip.text = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
ggsave(file = paste0(output,name,"_GO.pdf"), plot, width = 15, height = 18)


