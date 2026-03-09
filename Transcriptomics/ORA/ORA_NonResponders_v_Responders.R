##########################################################
# ORA:                                                   #
# KEGG                                                   #                        
# Graft Treatment: 3 and 6 week R v NR                   #
#                                                        #
# Use DEGs for gene sets:                                #
# visit MsigDB for more genesets:                        #
# http://www.gsea-msigdb.org/gsea/index.jsp              #
#                                                        #
# Non-Responders v Responders                            #
# Avg 3-6 week response                                  #
# Technical, lymphoma + HPV samples removed              #
#                                                        #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 05-11-25                               #
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
library(ReactomePA)
library(msigdbr)
library(enrichplot)
library(org.Hs.eg.db)

set.seed(1234)

########################################################## 3 WEEK RESPONSE ##########################################################
########################################################## READ DATA ##########################################################

# Read DEG DESeq2 results, 
DESeq_results <- read.csv("./2_pipelines/DESeq2/Filtered/Response/NRvR/DEG_Definitive_Response_3_6wkNonResponder_v_Responder.csv", header = T, row.names = 1) 
DESeq_results <- DESeq_results[!is.na(DESeq_results$padj),] # remove NA values for padj

DESeq_results_up <- read.csv("./2_pipelines/DESeq2/Filtered/Response/NRvR/DEG_UP_Definitive_Response_3_6wkNonResponder_v_Responder.csv", header = T, row.names = 1) 
DESeq_results_up <- DESeq_results_up[!is.na(DESeq_results_up$padj),] # remove NA values for padj

DESeq_results_down <- read.csv("./2_pipelines/DESeq2/Filtered/Response/NRvR/DEG_DOWN_Definitive_Response_3_6wkNonResponder_v_Responder.csv", header = T, row.names = 1) 
DESeq_results_down <- DESeq_results_down[!is.na(DESeq_results_down$padj),] # remove NA values for padj

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/Enrichment_Analysis/ORA/Filtered/NRvR/"
output <- "./3_output/Enrichment_Analysis/ORA/Filtered/NRvR/"
name <- "ORA_NRvR"
title <- "Basal HNSCC: Average 3-6 week NonResponders v Responders"


dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

########################################################## ORA using Cluster Profiler ##########################################################

# Get gene names from sig DEG results 
gene_symbols = rownames(DESeq_results)
gene_symbols_up = rownames(DESeq_results_up)
gene_symbols_down = rownames(DESeq_results_down)

# Convert gene symbols to Entrez IDs.
gene_entrez <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez_up <- mapIds(org.Hs.eg.db, keys = gene_symbols_up, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
gene_entrez_down <- mapIds(org.Hs.eg.db, keys = gene_symbols_down, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")


########################################################## GO ALL DEGS ##########################################################


####### KEGG Analysis UP ####### 
kegg_ora_up <- enrichKEGG(gene = gene_entrez_up, 
                          organism = 'hsa', 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.1, 
                          qvalueCutoff = 0.1)

# Visualize KEGG results
pdf(file = paste0(output,"KEGG/",name,"_KEGG_Up_DotPlot.pdf"), width = 10, height = 10)
dotplot(kegg_ora_up, showCategory = 20,
        title = paste0("KEGG Up DEGs: ", title))
dev.off()

pdf(file = paste0(output,"KEGG/",name,"_KEGG_Up_BarPlot.pdf"), width = 12, height = 12)
barplot(kegg_ora_up, showCategory = 20,
        title = NULL,
        font.size = 16) + 
  ggtitle(paste0(title)) + 
  labs(subtitle = "KEGG Up DEGs:") +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
dev.off()

# 1. Create the plot
# Note: font.size inside barplot() sets the base for the whole object
p_kegg <- barplot(kegg_ora_up, 
                  showCategory = 20, 
                  title = NULL,
                  font.size = 7) + 
  # 3. Apply the clean publication theme
  theme_bw(base_size = 7) +
  # Wrap labels at 50 characters since the plot is wide (15cm)
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  theme(
    plot.subtitle = element_text(size = 7, face = "italic"),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 7, color = "black"),
    
    # Legend position remains right by default, but we'll ensure it's there
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    
    # Grid and Border Cleanup
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  )

# 4. Save to your requested 9 x 15 cm dimensions
ggsave(
  filename = paste0(output, "KEGG/", name, "_KEGG_Up_PAPER.pdf"), 
  plot = p_kegg, 
  width = 14, 
  height = 12, 
  units = "cm",
  device = cairo_pdf
)


####### KEGG Analysis DOWN ####### 
kegg_ora_down <- enrichKEGG(gene = gene_entrez_down, 
                            organism = 'hsa', 
                            pAdjustMethod = "BH", 
                            pvalueCutoff = 0.1, 
                            qvalueCutoff = 0.1)

# Visualize KEGG results
pdf(file = paste0(output,"KEGG/",name,"_KEGG_Down_DotPlot.pdf"), width = 10, height = 10)
dotplot(kegg_ora_down, showCategory = 20,
        title = paste0("KEGG Down DEGs: ", title))
dev.off()

pdf(file = paste0(output,"KEGG/",name,"_KEGG_Down_BarPlot.pdf"), width = 12, height = 12)
barplot(kegg_ora_down, showCategory = 20,
        title = NULL,
        font.size = 16) + 
  ggtitle(paste0(title)) + 
  labs(subtitle = "KEGG Down DEGs:") +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
dev.off()


# Note: font.size inside barplot() sets the base for the whole object
p_kegg <- barplot(kegg_ora_down, 
                  showCategory = 20, 
                  title = NULL,
                  font.size = 7) + 
  # 3. Apply the clean publication theme
  theme_bw(base_size = 7) +
  # Wrap labels at 50 characters since the plot is wide (15cm)
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  theme(
    plot.subtitle = element_text(size = 7, face = "italic"),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 7, color = "black"),
    
    # Legend position remains right by default, but we'll ensure it's there
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    
    # Grid and Border Cleanup
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  )

# 4. Save to your requested 9 x 15 cm dimensions
ggsave(
  filename = paste0(output, "KEGG/", name, "_KEGG_Down_PAPER.pdf"), 
  plot = p_kegg, 
  width = 14, 
  height = 12, 
  units = "cm",
  device = cairo_pdf
)
