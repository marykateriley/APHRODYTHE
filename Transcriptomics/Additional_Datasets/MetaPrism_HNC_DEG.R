##########################################################
# DESeq2 - METAPRISM HNSC                                #
# PEARSON CORR                                           #
# SHH and KRT genes                                      #
#                                                        #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 16-02-26                               #
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
library(AnnotationDbi)
library(org.Hs.eg.db)

########################################################## READ DATA ##########################################################

# GES is RSEM from cBioPortal
GEX <- read.delim("./0_data/MetaPrism/Data_Table_1.rna_gene_counts.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
GEX <- as.matrix(GEX)

# TCGA metadata
meta <- read.delim("./0_data/MetaPrism/metaprism_2023_clinical_data_HNC_subset.txt", header = TRUE, stringsAsFactors = FALSE) # 6 samples do not have cna

# Read in SHH GMT and EGFR ligands Files
SHH <- read.gmt("./0_data/GMT/SHH/BIOCARTA_SHH_PATHWAY.v2023.2.Hs.gmt")
SHH_subset <- SHH$gene[SHH$gene %in% c("SHH", "PTCH1", "GLI1", "GLI2", "GLI3", "SMO","SUFU")]
KRT <- c("KRT8", "KRT18", "KRT5", "KRT14", "KRT1", "KRT10")
all_genes <- c(SHH_subset, KRT)
EGFR <- read.csv("./0_data/GMT/EGFR/EGFR_ligands.csv", header = T)
EGFR_subset <- EGFR$EGF_Ligands
all_genes <- c(SHH_subset, KRT, EGFR)

########################################################## SET VARIABLES ##########################################################
design_formula <- "~ Treatment" 
relevel_col <- "Treatment"
relevel_ref <- "No_Cetuximab"
contrast_column <- "Treatment"
primary_var <- "Cetuximab"
secondary_var <- "No_Cetuximab"
title <- "MetaPrism HNSC: "
subtitle <- "Comparing Treatments"
pipelines <- "./2_pipelines/MetaPrism/"
output <- "./3_output/MetaPrism/"

dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

########################################################## DATA CURATION ##########################################################
# GEX - is kallisto counts so non-integers
mx <- max(abs(GEX - round(GEX)), na.rm = TRUE)
mx

GEX <- round(GEX)
storage.mode(GEX) <- "integer"

# sanity check
stopifnot(all(GEX == floor(GEX)))

# Curate samples
colData <- meta %>%
  dplyr::mutate(
    Treatment = ifelse(
      grepl("Cetuximab", Treatments.Before.Biopsy, ignore.case = TRUE),
      "Cetuximab",
      "No_Cetuximab"
    )
  )

rownames(colData) <- colData$Patient.ID


## Only keep samples that were treated with cetuximab
#colData <- meta %>% 
#  dplyr::filter(grepl("cetuximab", Treatments.Before.Biopsy, ignore.case = TRUE))
#rownames(colData) <- colData$Patient.ID

# Change colnames from . to -
colnames(GEX) <- gsub("\\.", "-", colnames(GEX))
colnames(GEX)

# Only keep samples in GEX that are in meta
ids <- substr(colnames(GEX), 1, 7)
keep <- ids %in% colData$Patient.ID
GEX <- GEX[, keep]


# Only keep samples in colData that are in filtered GEX
ids <- substr(colnames(GEX), 1, 7)
colData <- colData[match(ids, colData$Patient.ID), , drop = FALSE]
identical(rownames(colData), ids)

rownames(colData) <- colnames(GEX)
stopifnot(identical(rownames(colData), colnames(GEX)))

# Are all samples in colnames of GEX present rownames of colData and are they in the same order
print(all(colnames(GEX) %in% rownames(colData)))
print(identical(rownames(colData), colnames(GEX)))
print(colData[match(colnames(GEX), rownames(colData)), ]) # order
print(identical(colData$Patient.ID, colnames(GEX))) # check again
print(all(colnames(GEX) %in% rownames(colData))) # check again


########################################################## DESEQ2 ##########################################################

# Create a DESeqDataSet to normalise and used for differential expression gene analysis 
dds <- DESeqDataSetFromMatrix(countData = GEX, 
                              colData = colData, 
                              design = as.formula(design_formula))

# Remove genes which have less than 10 counts
rows_less_than_10 <- which(rowSums(counts(dds)) < 10) # 17599 genes
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]


# Non-Responder as reference - this will look for genes enriched in response groups!
dds[[relevel_col]] <- relevel(dds[[relevel_col]], ref = relevel_ref) 

var_1 <- sum(dds[[relevel_col]] == primary_var) # 38
var_2 <- sum(dds[[relevel_col]] == secondary_var) # 11


# Run DESeq
dds_sub <- DESeq(dds)
resultsNames(dds_sub)
summary(dds_sub)

# Perform differential expression analysis compared to Non-Cetuximab
res <- results(dds_sub, contrast = c(contrast_column, primary_var, secondary_var))
summary(res)



# Write results
res <- as.data.frame(res) 

ensg <- sub("\\..*$", "", rownames(res))

symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensg,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res$ENSEMBL <- ensg
res$SYMBOL  <- unname(symbols)
write.csv(res, file = paste0(pipelines, "DESeq2_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)

# Order results table by the smallest padj value
resOrdered <- res[order(res$padj),]



normCount <- as.data.frame(counts(dds_sub, normalized = TRUE))

normCount$ENSEMBL <- sub("\\..*$", "", rownames(normCount))
normCount$SYMBOL <- unname(mapIds(
  org.Hs.eg.db,
  keys = normCount$ENSEMBL,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
))

normCount <- normCount %>%
  dplyr::select(ENSEMBL, SYMBOL, everything())

write.csv(
  normCount,
  file = paste0(pipelines, "Norm_counts_", contrast_column, primary_var, "_v_", secondary_var, "_annotated.csv"),
  quote = FALSE,
  row.names = FALSE
)



########################################################## MA PLOT ##########################################################

# Open a pdf device
pdf(paste0(output,"MAPlot_",contrast_column,primary_var,"_v_",secondary_var,".pdf"))

# Plot the MA plot
# Set the title size
par(cex.main = 0.85)
# Plot the MA plot with the correct main title
plotMA(resOrdered, main = paste0(contrast_column,": ",primary_var, " v ", secondary_var,"  
                                 " ,subtitle))

# Close the pdf device
dev.off()


########################################################## PCA PLOT ##########################################################

vsd <- vst(dds_sub, blind = FALSE)
saveRDS(vsd, file = paste0(pipelines,"vsd.RData"))

pcaData <- plotPCA(vsd, intgroup = paste0(contrast_column), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  ggtitle("Treatment") +
  labs(subtitle = paste0(subtitle)) +
  scale_color_manual(name = "Treatment", values = c("Cetuximab" = "#5AAE61FF", "No_Cetuximab" = "#9970ABFF")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0, size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15))
ggsave(paste0(output, "PCA_", contrast_column,"_",primary_var, "v", secondary_var,".pdf"), plot = p, width = 8, height = 8)


########################################################## VOLCANO PLOT ##########################################################

# padj < 0.05 #
res_filtered <- res %>% 
  filter(!is.na(padj)) %>% 
  mutate(threshold = padj < 0.1)

res_filtered$diffexpressed <- "NO"
res_filtered$diffexpressed[res_filtered$log2FoldChange > 0.585 & res_filtered$padj < 0.1] <- "UP"
res_filtered$diffexpressed[res_filtered$log2FoldChange < -0.585 & res_filtered$padj < 0.1] <- "DOWN"

mycolors <- c("#1f78B4", "black","#FF7F00")
names(mycolors) <- c("DOWN", "NO", "UP")

res_filtered$delabel <- NA
res_filtered$delabel[res_filtered$diffexpressed != "NO"] <- res_filtered$SYMBOL[res_filtered$diffexpressed != "NO"]


# Calculate total counts of up-regulated and down-regulated genes
total_up_genes <- sum(res_filtered$diffexpressed == "UP")
total_down_genes <- sum(res_filtered$diffexpressed == "DOWN")
total_up_down_genes <- sum(res_filtered$diffexpressed == "UP" | res_filtered$diffexpressed == "DOWN")

# Create a ranking score for top genes
res_filtered <- res_filtered %>%
  mutate(combined_score = abs(log2FoldChange) / padj)  # Higher score = stronger fold change & significance

# Identify the top 20 up-regulated and down-regulated genes
#top_up_genes <- res_filtered %>%
#  filter(diffexpressed == "UP" & padj < 0.05) %>%
#  top_n(20, abs(log2FoldChange))

top_up_genes <- res_filtered %>%
  filter(diffexpressed == "UP" & padj < 0.1) %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = 20)
write.csv(top_up_genes, file = paste0(pipelines, "TOP_DEG_UP_",contrast_column,primary_var,"_v_",secondary_var,".csv"), quote=F, row.names = TRUE)

#top_down_genes <- res_filtered %>%
#  filter(diffexpressed == "DOWN" & padj < 0.05) %>%
#  top_n(20, abs(log2FoldChange))

top_down_genes <- res_filtered %>%
  filter(diffexpressed == "DOWN" & padj < 0.1) %>%
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
  geom_hline(yintercept = -log10(0.1), col = "#0BDA51") +
  geom_vline(xintercept = c(-0.585, 0.585), col = "#0BDA51") +
  ggtitle(paste0(contrast_column,": ",primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") ")) +
  labs(subtitle = paste0(subtitle, " padj <0.1, Fold Change 1.5")) +
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
  # Add total samples label if total_up_genes and total_down_genes are not zero
  annotate("text", x = xmax, y = ymax, label = ifelse(total_up_genes > 0, paste("UP DEGs:", total_up_genes), ""), vjust = -1, col = "#FF7F00", size = 6, fontface = "bold") +
  annotate("text", x = xmin, y = ymax, label = ifelse(total_down_genes > 0, paste("DOWN DEGs:", total_down_genes), ""), vjust = -1, col = "#1f78B4", size = 6, fontface = "bold") +
  # Add total samples label
  annotate("text", x = Inf, y = -Inf, label = paste("Total DEGs:", total_up_down_genes),
           hjust = 1, vjust = -0.5,size = 6, col = "black", fontface = "bold") 
ggsave(paste0(output, "Volcano_", contrast_column,"_",primary_var, "v", secondary_var,"_padj0.05.pdf"), plot = p, width = 11, height = 11)






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

sym_map <- res %>%
  dplyr::select(ENSEMBL, SYMBOL) %>%
  dplyr::distinct()

vsd_mat <- assay(vsd)
# Transform the data to only include significant genes
ensg <- sub("\\..*$", "", rownames(vsd_mat))

vsd_df <- as.data.frame(vsd_mat)
vsd_df$ENSEMBL <- ensg

vsd_df <- dplyr::left_join(vsd_df, sym_map, by = "ENSEMBL")

# remove genes without symbol
vsd_df <- vsd_df %>% dplyr::filter(!is.na(SYMBOL), SYMBOL != "")

# collapse duplicate symbols by mean
mat <- as.matrix(vsd_df[, colnames(vsd_mat)])
sym <- vsd_df$SYMBOL

vsd_bySymbol <- sapply(split(seq_along(sym), sym), function(idx) {
  colMeans(mat[idx, , drop = FALSE])
})

ind <- t(vsd_bySymbol)



######## ALL DEGS ######## 

stopifnot(all(colnames(ind) %in% rownames(colData)))

# Subset annotation data to include only samples present in ind
annotations <- colData[match(colnames(ind), rownames(colData)), "Treatment", drop = FALSE]
annotations <- data.frame(annotations)

# Set column names of ind to match row names 
colnames(ind) <- rownames(annotations)
print(identical(rownames(annotations), colnames(ind))) # check if row names match

# annotation colors: include the function for this column
ann_colors <- list(
  Treatment = c("Cetuximab" = "purple", "No_Cetuximab" = "black"))



range(ind)

# Generate color palette for heatmap
breaks <- c(-3, 0, 3)
colors <- colorRamp2(breaks = breaks, colors = c("#1f78B4", "#F7F7F7FF", "#FF7F00"))

source("./1_code/RFunctions/new_heatmap.R")

########################################################## HEAT MAP SHH GENES ##########################################################

# Extract the gene names for the SHH pathway
SHH_genes <- SHH[["gene"]]  # Adjust if your GMT structure is different

# Ensure it's a character vector
SHH_genes <- unlist(SHH_genes)

common_genes <- intersect(SHH_genes, rownames(ind))

SHH_GLI2_data <- ind[common_genes, , drop = FALSE]

pdf(paste0(output, "Heatmap_BIOCARTA_SHH_PATHWAY_", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"), width = 9, height = 8)
new_heatmap(SHH_GLI2_data,
            scale = 'row',
            main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
            fontsize = 8,
            color = colors,
            annotation_col = annotations,
            annotation_colors = ann_colors,  
            annotation_names_col = FALSE,
            cluster_cols = TRUE, 
            cluster_rows = TRUE,  
            show_rownames = TRUE,  
            show_colnames = FALSE,
            name = "Z-score")
dev.off()




########################################################## HEAT MAP KERATIN MARKERS ##########################################################

# Extract  keratin markers
KRT_markers <- grepl("^KRT", rownames(ind))
KRT_genes <- ind[KRT_markers, , drop = FALSE]


pdf(paste0(output, "Heatmap_Keratin_Markers_", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"), width = 8, height = 8)
new_heatmap(KRT_genes,
            scale = 'row',
            main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
            fontsize = 10,
            color = colors,
            annotation_col = annotations,
            annotation_colors = ann_colors,  
            annotation_names_col = FALSE,
            cluster_cols = TRUE, 
            cluster_rows = TRUE,  
            show_rownames = TRUE,  
            show_colnames = FALSE,  
            name = "Z-score")

dev.off()



########################################################## HEAT MAP SHH and KERATIN MARKERS ##########################################################
## Not significant only ##
# Extract  keratin markers
SHH_KRT <- c("KRT8", "KRT18", "KRT5", "KRT14", "KRT1", "KRT10",
             "SHH", "SMO", "PTCH1", "GLI1", "GLI2", "GLI3", "SUFU")
SHH_KRT_genes <- ind[SHH_KRT, , drop = FALSE]


pdf(paste0(output, "Heatmap_SHH_KRT_Genes", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"), width = 8, height = 8)
new_heatmap(SHH_KRT_genes,
            scale = 'row',
            main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
            fontsize = 10,
            color = colors,
            annotation_col = annotations,
            annotation_colors = ann_colors,  
            annotation_names_col = FALSE,
            cluster_cols = TRUE, 
            cluster_rows = TRUE,  
            show_rownames = TRUE,  
            show_colnames = FALSE,  
            name = "Z-score")

dev.off()

########################################################## HEAT MAP SHH and EGFR MARKERS ##########################################################
## Not significant only ##
# Extract  keratin markers
# Ensure it's a character vector
plot_all_genes <- unlist(all_genes)

common_genes <- intersect(plot_all_genes, rownames(ind))

plot_all_genes <- ind[common_genes, , drop = FALSE]

pdf(paste0(output, "Heatmap_SHH_KRT_EGFR_Genes", contrast_column, "_", primary_var, "v", secondary_var, ".pdf"), width = 8, height = 8)
new_heatmap(plot_all_genes,
            scale = 'row',
            main = paste0(title, primary_var," (",var_1,") ", "v ", secondary_var," (",var_2,") "),
            fontsize = 10,
            color = colors,
            annotation_col = annotations,
            annotation_colors = ann_colors,  
            annotation_names_col = FALSE,
            cluster_cols = TRUE, 
            cluster_rows = TRUE,  
            show_rownames = TRUE,  
            show_colnames = FALSE,  
            name = "Z-score")

dev.off()

pdf(paste0(output, "PAPER_SHH_KRT_EGFR_Genes.pdf"),  width = 3.4, height = 4.5)
ComplexHeatmap::pheatmap(plot_all_genes,
                         scale = 'row',
                         main = paste0(title, primary_var,"(",var_1,") ", "v ", secondary_var,"(",var_2,") "),
                         fontsize = 6.5,
                         color = colors,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         annotation_names_col = FALSE,
                         show_rownames = TRUE,
                         show_colnames = FALSE,
                         treeheight_row = 10,
                         treeheight_col = 0,
                         legend = FALSE,
                         annotation_legend = FALSE)

dev.off()


site_col_vector <- ann_colors$Treatment

# 2. Create the Legend object
site_legend <- Legend(
  at = names(site_col_vector),
  title = "Treatment",
  legend_gp = gpar(fill = site_col_vector),
  direction = "horizontal",
  nrow = 2, # Creates a 2-column layout for 7 items
  title_position = "topleft",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 7)
)

ggsave(
  filename = paste0(output, "MetaPrism_CetuxTreatment.pdf"), 
  plot = site_legend, 
  width = 2.8, 
  height = 1.5, 
  units = "cm",
  device = cairo_pdf
)

