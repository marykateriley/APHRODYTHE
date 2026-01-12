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

# Gather C6 (onogenic signature) gene sets
c6_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, gene_symbol)

# Gather KEGG gene sets
K_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
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



############### STEP 5: plot specific geneset as a GSEAPLOT
# DOWNREGULATED
#pdf(file = paste0(output,name,"_HALLMARK_OXIDATIVE_PHOSPHORYLATION.pdf"), width = 13, height = 8, units = "in", res = 500)
#gseaplot2(
#gseaResult,
#"HALLMARK_OXIDATIVE_PHOSPHORYLATION",                         #insert the name of the specific pathway to plot 
#title = paste0(title,": HALLMARK_OXIDATIVE_PHOSPHORYLATION"),                 #change title 
#color = "green",
#base_size = 15,
#rel_heights = c(1.5, 0.5, 1),
#subplots = 1:3,
#pvalue_table = TRUE,
#ES_geom = "line"
#)
#dev.off()


## C6 (onogenic signature) ##
c6_gseaResult <-GSEA(geneList, exponent = 1, pvalueCutoff = 1,
                     pAdjustMethod = "BH", TERM2GENE = c6_t2g, verbose = TRUE, seed = TRUE) #can set seed=True or just run once and save result for further plots! 
c6_gseaRes_df <- as.data.frame(c6_gseaResult)

write.csv(c6_gseaRes_df, file=paste0(pipelines,name,"_C6.csv"))

# Ridgeplot
pdf(file = paste0(output,name,"_C6_RidgePlot.pdf"), width = 13, height = 13)
ridgeplot(c6_gseaResult)
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



############### KEGG ###############

original_gene_list = DESeq_results[,4]
## adding the genenames from DESeq results
names(original_gene_list) = as.character(rownames(DESeq_results))

# Convert gene symbols to ENTREZ IDs
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Remove duplicate IDs
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
# Ensure that rownames are converted to a column in DESeq_results
DESeq_results$SYMBOL <- rownames(DESeq_results)

# New df2 containing genes successfully mapped
df2 <- DESeq_results[DESeq_results$SYMBOL %in% dedup_ids$SYMBOL, ]

# New column in df2 with the corresponding ENTREZ IDs
df2$ENTREZID <- dedup_ids$ENTREZID

# Vector
kegg_gene_list <- df2$stat

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$ENTREZID

# Omit NA 
kegg_gene_list<-na.omit(kegg_gene_list)

# Sort descending
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_enrichment <- gseKEGG(geneList = kegg_gene_list,
                           organism = 'hsa',
                           pvalueCutoff = 0.1,
                           verbose = TRUE,
                           pAdjustMethod = "BH")
write.csv(gse, file=paste0(pipelines,name,"_KEGG.csv"))

plot <- dotplot(kegg_enrichment, showCategory = 10, title = paste0("Enriched KEGG Pathways: ", title), split=".sign") + facet_grid(.~.sign)
ggsave(file = paste0(output,name,"_KEGG.pdf"), plot, width = 10, height = 13)
