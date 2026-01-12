##########################################################
# DESeq2                                                 #
# PDX v Humans                                           #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 05-05-25                               #
##########################################################

setwd("/Users/mary-kateriley/Library/CloudStorage/OneDrive-Queen'sUniversityBelfast/Data Analysis/HNC/Basal PDX/Reanalysis Basal '24")

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
library(circlize)
library(rstatix) # wilcox_test
library(ggpubr) # stat_pvalue_manual
library(RColorBrewer)

########################################################## SET VARIABLES ##########################################################

design_formula <- "~ Case_ID + Model"
relevel_ref <- "Human"
contrast_column <- "Model"
primary_var <- "PDX"
secondary_var <- "Human"
pipelines <- "./2_pipelines/PDX_Human/DESeq2/Unfiltered/"
output <- "./3_output/PDX_Human/DESeq2/Unfiltered/"
title <- "HNSCC: "
subtitle <- "Basal HNSCC PDX v Human: Filtered Technical; Lymphoma; HPV"


########################################################## READ DATA ##########################################################

GEX_pdx <- read.csv("./0_data/DESeq2/PRX/GEX/All_RemovedSamples_Unfiltered_RawCounts_GraftHNSCC_270424.csv", row.names = 1)
colData_pdx <- read.csv("./0_data/DESeq2/PRX/Metadata/All_RemovedSamples_Unfiltered_GraftHNSCC_colData_240724.csv", header = T, row.names = 1)
colData_pdx$Model <- "PDX"

GEX_human <- read.csv("./0_data/DESeq2/Human/GEX/All_RemovedSamples_Unfiltered_RawCounts_HumanHNSCC_051124.csv", row.names = 1)
colData_human <- read.csv("./0_data/DESeq2/Human/Metadata/All_RemovedSamples_Unfiltered_HumanHNSCC_colData_051124.csv", header = T, row.names = 1)
colData_human$Model <- "Human"


# Combine the counts 
GEX <- cbind(GEX_pdx, GEX_human)
write.csv(GEX, file = paste0(pipelines, "PDXandHuman_Merged_RawCounts_Filtered.csv"), quote=F, row.names = TRUE)


# Combine the metadata 
colData_pdx <- colData_pdx[, !(colnames(colData_pdx) == "Replicate")] # this is not present in human colData so need to remove first
colData <- rbind(colData_pdx, colData_human)
write.csv(colData, file = paste0(pipelines, "PDXandHuman_Merged_colData_Filtered.csv"), quote=F, row.names = TRUE)

########################################################## DATA CURATION ##########################################################

# Are all samples in colnames of GEX present rownames of colData and are they in the same order
all_samples_present <- all(colnames(GEX) %in% rownames(colData))
are_same_order <- identical(rownames(colData), colnames(GEX))
colData <- colData[match(colnames(GEX), rownames(colData)), ] # order
are_same_order <- identical(colData$Sample_ID, colnames(GEX)) # check again
all_samples_present <- all(colnames(GEX) %in% rownames(colData)) # check again

colData$Site.of.Primary[colData$Site.of.Primary == ""] <- "Unknown"
colData$Definitive_Response_3_6wk[is.na(colData$Definitive_Response_3_6wk)] <- "No_data"
colData$Response_3_6wk[is.na(colData$Response_3_6wk)] <- "No_data"
colData$percentage_vol_change_3_6wk[is.na(colData$percentage_vol_change_3_6wk)] <- "No_data"

# Set as factors for DESeq2
colData$Model <- factor(colData$Model)
colData$Case_ID <- factor(colData$Case_ID)


########################################################## PDX v Human DESEQ2 ##########################################################

# Create a DESeqDataSet to normalise and used for differential expression gene analysis 
dds <- DESeqDataSetFromMatrix(countData = GEX, 
                              colData = colData, 
                              design = as.formula(design_formula))

# Remove genes which have less than 10 counts
rows_less_than_10 <- which(rowSums(counts(dds)) < 10) # 19679 genes
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]


# Human as reference - this will look for genes enriched in PDX groups!
dds$Model <- relevel(dds$Model, ref = relevel_ref) 

var_1 <- sum(dds$Model == primary_var) # 67
var_2 <- sum(dds$Model == secondary_var) # 47

# Run DESeq
dds_sub <- DESeq(dds)
resultsNames(dds_sub)
summary(dds_sub)

# Perform differential expression analysis compared to Human
res <- results(dds_sub, contrast = c(contrast_column, primary_var, secondary_var))
summary(res)

# Order results table by the smallest padj value
resOrdered <- res[order(res$padj),]


# Write results
res <- as.data.frame(res) 
write.csv(res, file = paste0(pipelines, "DESeq2_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)
#res <- read.csv(file = paste0(pipelines, "DESeq2_", contrast_column, primary_var, "_v_", secondary_var, ".csv"), row.names = 1)

normCount <- counts(dds_sub, normalized = TRUE)
write.csv(normCount, file = paste0(pipelines, "Norm_counts_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)
#normCount <- read.csv(file = paste0(pipelines, "Norm_counts_", contrast_column, primary_var, "_v_", secondary_var, ".csv"), row.names = 1)

########################################################## MA PLOT ##########################################################

pdf(paste0(output,"MAPlot_",contrast_column,primary_var,"_v_",secondary_var,".pdf"))
par(cex.main = 0.85)
DESeq2::plotMA(resOrdered, main = paste0(contrast_column,": ",primary_var, " v ", secondary_var,"  
                                 " ,subtitle))
dev.off()

########################################################## PCA PLOT ##########################################################

vsd <- vst(dds_sub, blind = FALSE)
saveRDS(vsd, file = paste0(pipelines,"vsd.RData"))
#vsd <- readRDS(file = paste0(pipelines, "vsd.RData"))

vsd_blind <- vst(dds_sub, blind = TRUE)
saveRDS(vsd_blind, file = paste0(pipelines,"vsdBlind.RData"))
#vsd_blind <- readRDS(file = paste0(pipelines, "vsdBlind.RData"))

pcaData <- plotPCA(vsd_blind, intgroup = paste0(contrast_column), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  ggtitle(paste0(contrast_column, ": ", primary_var, " (", var_1, ") v ", secondary_var, " (", var_2, ")")) +
  labs(subtitle = subtitle) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Model", values = c("PDX" = "#f5b041", "Human" = "#1abc9c")) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0, size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15))
ggsave(paste0(output, "BLIND_PCA_", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"),plot = p, width = 8, height = 8)

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  geom_text(aes(label = Case_ID), vjust = -1, size = 3, show.legend = FALSE) + 
  ggtitle(paste0(contrast_column, ": ", primary_var, " (", var_1, ") v ", secondary_var, " (", var_2, ")")) +
  labs(subtitle = subtitle) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Model", values = c("PDX" = "#f5b041", "Human" = "#1abc9c")) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0, size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15))
ggsave(paste0(output, "BLIND_PCA_CaseID", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"), plot = p, width = 8, height = 8)




# Get Case_ID pairs
paired_lines <- merge(pcaData, pcaData, by = "Case_ID") %>%
  filter(group.x != group.y) %>%
  mutate(pair_id = paste(pmin(Sample_ID.x, Sample_ID.y),
                         pmax(Sample_ID.x, Sample_ID.y), sep = "_")) %>%
  distinct(pair_id, .keep_all = TRUE)

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_segment(data = paired_lines,
               aes(x = PC1.x, y = PC2.x, xend = PC1.y, yend = PC2.y),
               color = "gray50", alpha = 0.6, linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = Case_ID), vjust = -1, size = 3, show.legend = FALSE) +
  ggtitle(paste0(contrast_column, ": ", primary_var, " (", var_1, ") v ", secondary_var, " (", var_2, ")")) +
  labs(subtitle = subtitle) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(name = "Model", values = c("PDX" = "#f5b041", "Human" = "#1abc9c")) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0, size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15))
ggsave(paste0(output, "BLIND_PCA_LinkedCaseIDs_", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"),plot = p, width = 8, height = 8)


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

res_filtered <- res_filtered %>%
  mutate(
    padj = ifelse(padj < 1e-300, 1e-300, padj),
    pvalue = ifelse(pvalue < 1e-300, 1e-300, pvalue))  # need to change as a padjust/pvalue of 0.000 will return Inf score
  

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
ymax <- max(-log10(res_filtered$padj), na.rm = TRUE) + 10
xmax <- (2 + max(res_filtered$`log2FoldChange`)) / 2 #choose the middle value
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
  annotate("text", x = xmax / 0.7, y = ymax, label = ifelse(total_up_genes > 0, paste("UP DEGs:", total_up_genes), ""), vjust = -1, col = "#FF7F00", size = 6, fontface = "bold") +
  annotate("text", x = xmin, y = ymax, label = ifelse(total_down_genes > 0, paste("DOWN DEGs:", total_down_genes), ""), vjust = -1, col = "#1f78B4", size = 6, fontface = "bold") +
  annotate("text", x = Inf, y = -Inf, label = paste("Total DEGs:", total_up_down_genes),
           hjust = 1, vjust = -0.5,size = 6, col = "black", fontface = "bold") 
p <- p + scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)))
ggsave(paste0(output, "Volcano_", contrast_column,"_",primary_var, "v", secondary_var,"_padj0.05.pdf"), plot = p, width = 12, height = 14)




########################################################## HEAT MAPS ##########################################################

# Get a list of differential expressed genes
sig <- res_filtered %>% 
  filter(diffexpressed != "NO")
write.csv(sig, file = paste0(pipelines, "DEG_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)

sig_up <- res_filtered %>% 
  filter(diffexpressed == "UP")
write.csv(sig_up, file = paste0(pipelines, "DEG_UP_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)

sig_down <- res_filtered %>% 
  filter(diffexpressed == "DOWN")
write.csv(sig_down, file = paste0(pipelines, "DEG_DOWN_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)

# Transform the data to only include significant genes
ind <- assay(vsd)[rownames(assay(vsd)) %in% row.names(sig), ]


######## ALL DEGS ######## 

# Subset annotation data to include only samples present in ind
annotations <- colData[colnames(ind), c("Model", "Definitive_Response_3_6wk", "Response_3_6wk", "Sex", "Site.of.Primary")]
annotations <- data.frame(annotations)

# Set column names of ind to match row names 
colnames(ind) <- rownames(annotations)
identical(rownames(annotations), colnames(ind)) # check if row names match



# Set annotation colours
ann_colors <- list(
  Model = c("PDX" = "#f5b041", "Human" = "#1abc9c"),
  Definitive_Response_3_6wk = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF", "No_data" = "#666666FF"),
  Response_3_6wk = c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF", "No_data" = "#666666FF"),
  Sex = c("M" = "#A6CEE3", "F" = "#FB9A99"),
  Site.of.Primary = c("Larynx" = "#1B9E77FF", "Oropharynx" = "#66A61EFF", "Hypopharynx" = "#7570B3FF", "Oral Cavity" = "#E6AB02FF", "Larynx + Oral Cavity" = "#D95F02FF", "Oral Cavity + Oropharynx" = "#E7298AFF", "Unknown" = "#666666FF")
)


# Generate colour palette for heatmap
breaks <- c(-3, 0, 3)
colors <- colorRamp2(breaks = breaks, colors = c("#1f78B4", "#F7F7F7FF", "#FF7F00"))


# Save legend horizontally 
legend_list <- lapply(names(ann_colors), function(name) {
  Legend(title = name,
         labels = names(ann_colors[[name]]),
         legend_gp = gpar(fill = unname(ann_colors[[name]])),
         direction = "horizontal")
})

# Combine legends into one horizontal block
annotation_legend <- do.call(packLegend, c(
  legend_list,
  list(direction = "horizontal", gap = unit(1.5, "cm"))))

# Estimate better dimensions
pdf(paste0(output, "AnnotationLegend_Horizontal.pdf"), width = 10, height = 10)
draw(annotation_legend, x = unit(0.5, "npc"), just = "center")
dev.off()

########################################################## HEAT MAP ALL DEGS #########################################################

# Natural Clustering
pdf(paste0(output, "Heatmap_", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"))
ComplexHeatmap::pheatmap(ind,
                         scale = 'row',
                         main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
                         fontsize = 10,
                         color = colors,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()

# No labels Natural Clustering
pdf(paste0(output, "NoLabel_Heatmap_", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"))
ComplexHeatmap::pheatmap(ind,
                         scale = 'row',
                         main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
                         fontsize = 10,
                         color = colors,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         annotation_legend = FALSE,
                         annotation_names_col = FALSE,
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()


# Natural Clustering
pdf(paste0(output, "Heatmap_CaseIDs", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"), width = 12, height = 8)
ComplexHeatmap::pheatmap(ind,
                         scale = 'row',
                         #main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
                         fontsize = 10,
                         color = colors,
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()


# Column Split
pdf(paste0(output, "ColumnSplit_Heatmap", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"))
ComplexHeatmap::pheatmap(ind,
                         scale = 'row',
                         main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
                         fontsize = 10,
                         color = colors,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         column_split = as.factor(annotations$Definitive_Response_3_6wk),
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()

# Column Split Model
pdf(paste0(output, "ColumnSplit_Model_Heatmap", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"))
ComplexHeatmap::pheatmap(ind,
                         scale = 'row',
                         main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
                         fontsize = 10,
                         color = colors,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         column_split = as.factor(annotations$Model),
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()


