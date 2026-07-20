##########################################################
# DESeq2                                                 #
# Non-Responders v Responders                            #
# Avg 3-6 week response                                  #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 04-11-25                               #
##########################################################

setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(circlize)
library(rstatix) # wilcox_test
library(ggpubr) # stat_pvalue_manual


########################################################## READ DATA ##########################################################

GEX <- read.csv("./0_data/DESeq2/PRX/GEX/All_RemovedSamples_Unfiltered_RawCounts_GraftHNSCC_270424.csv", row.names = 1)

colData <- read.csv("./0_data/DESeq2/PRX/Metadata/All_RemovedSamples_Unfiltered_GraftHNSCC_colData_240724.csv", header = T, row.names = 1)



########################################################## SET VARIABLES ##########################################################

design_formula <- "~ Definitive_Response_3_6wk" # Can't include patient as each sample from that patient will have the same response - will not be full rank
relevel_ref <- "Responder"
contrast_column <- "Definitive_Response_3_6wk"
primary_var <- "NonResponder"
secondary_var <- "Responder"
pipelines <- "./2_pipelines/DESeq2/Filtered/Response/NRvR/"
output <- "./3_output/DESeq2/Filtered/Response/NRvR/"
title <- "PDX HNSCC: "
subtitle <- "Basal HNSCC PDX: Filtered Technical; Lymphoma; HPV"

dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

########################################################## DATA CURATION ##########################################################

# Are all samples in colnames of GEX present rownames of colData and are they in the same order
all_samples_present <- all(colnames(GEX) %in% rownames(colData))
are_same_order <- identical(rownames(colData), colnames(GEX))
colData <- colData[match(colnames(GEX), rownames(colData)), ] # order
are_same_order <- identical(colData$Sample_ID, colnames(GEX)) # check again
all_samples_present <- all(colnames(GEX) %in% rownames(colData)) # check again

# Replace blanks in Site.of.Primary with unknown
colData$Site.of.Primary[colData$Site.of.Primary == ""] <- "Unknown"

# Set as factors for DESeq2
colData$Site.of.Primary <- factor(colData$Site.of.Primary)
colData$Definitive_Response_3_6wk <- factor(colData$Definitive_Response_3_6wk)


########################################################## NR v R DESEQ2 ##########################################################

# Remove any sample that does not have Definitive_Response_3wk (NAs) from colData
colData_filtered <- colData  %>%
  filter(!is.na(Definitive_Response_3_6wk))

# Ensure all samples in GEX match those samples that have a 3-week response in colData_filtered
GEX_filtered <- GEX[, colnames(GEX) %in% rownames(colData_filtered)]

# Set row names of colData_filtered to match the samples in GEX_filtered
rownames(colData_filtered) <- colnames(GEX_filtered)

# Create a DESeqDataSet to normalise and used for differential expression gene analysis 
dds <- DESeqDataSetFromMatrix(countData = GEX_filtered, 
                              colData = colData_filtered, 
                              design = as.formula(design_formula))

# Remove genes which have less than 10 counts
rows_less_than_10 <- which(rowSums(counts(dds)) < 10) # 24713 genes
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]


# Non-Responder as reference - this will look for genes enriched in response groups!
dds$Definitive_Response_3_6wk <- relevel(dds$Definitive_Response_3_6wk, ref = relevel_ref) 

var_1 <- sum(dds$Definitive_Response_3_6wk == primary_var) # 15
var_2 <- sum(dds$Definitive_Response_3_6wk == secondary_var) # 40

# Run DESeq
dds_sub <- DESeq(dds)
resultsNames(dds_sub)
summary(dds_sub)

# Perform differential expression analysis compared to Responder
res <- results(dds_sub, contrast = c(contrast_column, primary_var, secondary_var))
summary(res)

# Order results table by the smallest padj value
resOrdered <- res[order(res$padj),]

# Write results
res <- as.data.frame(res) 
write.csv(res, file = paste0(pipelines, "DESeq2_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)

normCount <- counts(dds_sub, normalized = TRUE)
write.csv(normCount, file = paste0(pipelines, "Norm_counts_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)


########################################################## MA PLOT ##########################################################

pdf(paste0(output,"MAPlot_",contrast_column,primary_var,"_v_",secondary_var,".pdf"))
par(cex.main = 0.85)
plotMA(resOrdered, main = paste0(contrast_column,": ",primary_var, " v ", secondary_var,"  
                                 " ,subtitle))
dev.off()

########################################################## PCA PLOT ##########################################################

vsd <- vst(dds_sub, blind = FALSE)
saveRDS(vsd, file = paste0(pipelines,"vsd.RData"))

vsd_blind <- vst(dds_sub, blind = TRUE)
saveRDS(vsd_blind, file = paste0(pipelines,"vsdBlind.RData"))

pcaData <- plotPCA(vsd_blind, intgroup = paste0(contrast_column), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  ggtitle("Definitive Response 3-6 week") +
  labs(subtitle = paste0(subtitle)) +
  scale_color_manual(name = "Response", values = c("Responder" = "#5AAE61FF", "NonResponder" = "#9970ABFF")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0, size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15))
ggsave(paste0(output, "BlindPCA_", contrast_column,"_",primary_var, "v", secondary_var,".pdf"), plot = p, width = 8, height = 8)

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  geom_text(aes(label = Case_ID), vjust = -1, size = 3, show.legend = FALSE) + 
  ggtitle("Definitive Response 3-6 week") +
  labs(subtitle = paste0(subtitle)) +
  scale_color_manual(name = "Response", values = c("Responder" = "#5AAE61FF", "NonResponder" = "#9970ABFF")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0, size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15))
ggsave(paste0(output, "BlindPCA_CaseID", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"), plot = p, width = 8, height = 8)


pcaData <- plotPCA(vsd_blind, intgroup = "Site.of.Primary", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  ggtitle("Primary Site") +
  labs(subtitle = paste0(subtitle)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Primary Site", values = c("Larynx" = "#1B9E77FF", "Oropharynx" = "#66A61EFF", "Hypopharynx" = "#7570B3FF", "Oral Cavity" = "#E6AB02FF", "Larynx + Oral Cavity" = "#D95F02FF", "Oral Cavity + Oropharynx" = "#E7298AFF", "Unknown" = "#666666FF")) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0, size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15))
ggsave(paste0(output, "BlindPCA_PrimarySite_",primary_var, "v", secondary_var,".pdf"), plot = p, width = 8, height = 8)




########################################################## VOLCANO PLOT ##########################################################

# padj < 0.05 #
res_filtered <- res %>% 
  filter(!is.na(padj)) %>% 
  mutate(threshold = padj < 0.05)

res_filtered$diffexpressed <- "NO"
res_filtered$diffexpressed[res_filtered$log2FoldChange > 0.585 & res_filtered$padj < 0.05] <- "UP"
res_filtered$diffexpressed[res_filtered$log2FoldChange < -0.585 & res_filtered$padj < 0.05] <- "DOWN"

mycolors <- c("#1f78B4", "black","#FF7F00")
names(mycolors) <- c("DOWN", "NO", "UP")

res_filtered$delabel <- NA
res_filtered$delabel[res_filtered$diffexpressed != "NO"] <- rownames(res_filtered)[res_filtered$diffexpressed != "NO"]


# Calculate total counts of up-regulated and down-regulated genes
total_up_genes <- sum(res_filtered$diffexpressed == "UP")
total_down_genes <- sum(res_filtered$diffexpressed == "DOWN")
total_up_down_genes <- sum(res_filtered$diffexpressed == "UP" | res_filtered$diffexpressed == "DOWN")

# Create a ranking score for top genes
res_filtered <- res_filtered %>%
  mutate(combined_score = abs(log2FoldChange) / padj)  # Higher score = stronger fold change & significance

top_up_genes <- res_filtered %>%
  filter(diffexpressed == "UP" & padj < 0.05) %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = 20)
write.csv(top_up_genes, file = paste0(pipelines, "TOP_DEG_UP_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)

top_down_genes <- res_filtered %>%
  filter(diffexpressed == "DOWN" & padj < 0.05) %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = 20)
write.csv(top_down_genes, file = paste0(pipelines, "TOP_DEG_DOWN_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)


# Calculate thresholds
ymax <- max(-log10(res_filtered$padj), na.rm = TRUE) + 1
xmax <- (0 + max(res_filtered$`log2FoldChange`)) / 2 #choose the middle value
xmin <- (0 + min(res_filtered$`log2FoldChange`)) / 2 #choose the middle value

p <- ggplot(data = res_filtered, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_point() + 
  geom_text_repel(data = rbind(top_up_genes, top_down_genes), aes(label = delabel), max.overlaps = 17, force_pull = 5, show.legend = FALSE) +
  theme_bw() +
  xlab("log2FoldChange") +
  scale_color_manual(values = mycolors) +
  geom_hline(yintercept = -log10(0.05), col = "#0BDA51") +
  geom_vline(xintercept = c(-0.585, 0.585), col = "#0BDA51") +
  ggtitle(paste0(contrast_column,": ",primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") ")) +
  labs(subtitle = paste0(subtitle, " padj <0.05, Fold Change 1.5")) +
  theme(axis.line = element_line(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size = 18,  colour = "black", face = "bold"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 18,  colour = "black", face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.title.position = "plot") +
  annotate("text", x = xmax, y = ymax, label = ifelse(total_up_genes > 0, paste("UP DEGs:", total_up_genes), ""), vjust = -1, col = "#FF7F00", size = 6, fontface = "bold") +
  annotate("text", x = xmin, y = ymax, label = ifelse(total_down_genes > 0, paste("DOWN DEGs:", total_down_genes), ""), vjust = -1, col = "#1f78B4", size = 6, fontface = "bold") +
  annotate("text", x = Inf, y = -Inf, label = paste("Total DEGs:", total_up_down_genes),
           hjust = 1, vjust = -0.5,size = 6, col = "black", fontface = "bold") 
ggsave(paste0(output, "Volcano_", contrast_column,"_",primary_var, "v", secondary_var,"_padj0.05.pdf"), plot = p, width = 11, height = 11)


