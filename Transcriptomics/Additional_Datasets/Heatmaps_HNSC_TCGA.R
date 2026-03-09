##########################################################
# HEATMAPS - TCGA HNSC                                   #
# DGEA                                                   #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 16-06-25                               #
##########################################################

setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(dplyr)
library(matrixStats)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(circlize)
library(rstatix) # wilcox_test
library(ggpubr) # stat_pvalue_manual

########################################################## READ DATA ##########################################################

# GES is RSEM from cBioPortal
GEX <- read.csv("./0_data/hnsc_tcga_pan_can_atlas_2018/HPVRem/rna_filtered.csv")
GEX_hugo <- GEX[!is.na(GEX$Hugo_Symbol) & GEX$Hugo_Symbol != "", ] # Remove rows with missing or empty Hugo symbols
GEX_hugo$Hugo_Symbol_unique <- make.unique(as.character(GEX_hugo$Hugo_Symbol)) # Create unique gene symbols by appending a suffix to duplicates
rownames(GEX_hugo) <- GEX_hugo$Hugo_Symbol_unique # Set rownames
GEX_hugo <- GEX_hugo[, !(colnames(GEX_hugo) %in% c("Hugo_Symbol", "Entrez_Gene_Id", "Hugo_Symbol_unique"))] # Remove the Hugo_Symbol + Entrez_ID columns
GEX_hugo[] <- lapply(GEX_hugo, function(x) as.numeric(as.character(x)))
GEX_hugo <- as.matrix(GEX_hugo)

# TCGA metadata
tcga_meta <- read.csv("./0_data/hnsc_tcga_pan_can_atlas_2018/HPVRem/curated_sitesrem_hpvrem_hnsc_tcga_pan_can_atlas_2018_clinical_data.csv", header = TRUE, stringsAsFactors = FALSE) # 6 samples do not have cna

# Read in SHH GMT and EGFR ligands Files
SHH <- read.gmt("./0_data/GMT/SHH/BIOCARTA_SHH_PATHWAY.v2023.2.Hs.gmt")
EGFR <- read.csv("./0_data/GMT/EGFR/EGFR_ligands.csv", header = T)


SHH_subset <- SHH$gene[SHH$gene %in% c("SHH", "PTCH1", "GLI1", "GLI2", "GLI3", "SMO","SUFU")]
KRT <- c("KRT8", "KRT18", "KRT5", "KRT14", "KRT1", "KRT10")
EGFR_subset <- EGFR$EGF_Ligands
all_genes <- c(SHH_subset, KRT, EGFR)

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/TCGA/"
output <- "./3_output/TCGA/"



########################################################## DATA CURATION ##########################################################

table(tcga_meta$cases.primary_site)

group_map <- list(
  "Oral Cavity" = c("Floor of mouth", "Gum", "Palate", "Other and unspecified parts of mouth", "Other and unspecified parts of tongue"),
  "Oropharynx" = c("Base of tongue", "Tonsil", "Oropharynx"),
  "Hypopharynx" = c("Hypopharynx"),
  "Larynx" = c("Larynx"),
  "Other" = c("Other and ill-defined sites in lip, oral cavity and pharynx")
)

# Function to map primary sites into broader categories
map_primary_site <- function(site) {
  for (group in names(group_map)) {
    if (site %in% group_map[[group]]) return(group)
  }
  return("Unclassified")
}

# Apply to your data
tcga_meta$Primary_Site_Group <- sapply(tcga_meta$cases.primary_site, map_primary_site)


# Fix expression matrix column names
mat_ids <- substr(colnames(GEX_hugo), 1, 12)  # Get first 12 chars
mat_ids <- gsub("\\.", "-", mat_ids)         # Replace . with -

# Update colnames of GEX_hugo to match tcga_meta format
colnames(GEX_hugo) <- mat_ids
GEX_hugo <- as.matrix(GEX_hugo)


# Filter metadata to only samples in GEX
GEX_meta <- tcga_meta[tcga_meta$Patient.ID %in% mat_ids, ]
# Match order to expression matrix
meta <- GEX_meta[match(colnames(GEX_hugo), GEX_meta$Patient.ID), ]
rownames(meta) <- meta$Patient.ID

identical(colnames(GEX_hugo), rownames(meta))  # should be TRUE

# Top 500 most variable genes across all samples
top_var_genes <- head(order(rowVars(as.matrix(GEX_hugo)), decreasing = TRUE), 500)
heat_mat <- GEX_hugo[top_var_genes, ]
dim(heat_mat)

########################################################## HEAT MAP TOP 500 VARIABLE GENES ##########################################################

# Subset annotation data to include only samples present in ind
annotations <- meta[colnames(heat_mat), c("Primary_Site_Group")]
annotations <- data.frame(Primary_Site_Group = annotations)

# Set column names of ind to match row names 
colnames(heat_mat) <- rownames(annotations)
identical(rownames(annotations), colnames(heat_mat)) # check if row names match



# Set annotation colours
ann_colors <- list(
  Primary_Site_Group = c(
    "Larynx" = "#1B9E77FF",
    "Oropharynx" = "#66A61EFF",
    "Hypopharynx" = "#7570B3FF",
    "Oral Cavity" = "#E6AB02FF",
    "Other" = "#666666FF"))

# Generate color palette for heatmap
breaks <- c(-3, 0, 3)
colors <- colorRamp2(breaks = breaks, colors = c("#1f78B4", "#F7F7F7FF", "#FF7F00"))


pdf(paste0(output, "Heatmap_Top500VariableGenes.pdf"))
ComplexHeatmap::pheatmap(heat_mat,
                         scale = 'row',
                         main = paste0("Pan Cancer Atlas HNSC - Top 500 Variable Genes"),
                         fontsize = 8,
                         color = colors,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()

pdf(paste0(output, "HeatmapPrimarySite_Top500VariableGenes.pdf"))
ComplexHeatmap::pheatmap(heat_mat,
                         scale = 'row',
                         main = paste0("Pan Cancer Atlas HNSC - Top 500 Variable Genes"),
                         fontsize = 8,
                         color = colors,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()



########################################################## HEAT MAP TOP SHH, KRT AND EGFR GENES ##########################################################

all_genes <- unlist(all_genes)
sum(all_genes %in% rownames(GEX_hugo))  # How many match?

# Keep only genes present in GEX_hugo
matched_genes <- intersect(all_genes, rownames(GEX_hugo))

# Subset expression matrix
gene_subset_expr <- GEX_hugo[matched_genes, ]


pdf(paste0(output, "HeatmapPrimarySite_SHH_KRT_EGFR_Genes.pdf"), width =7.6, height = 8)
ComplexHeatmap::pheatmap(gene_subset_expr,
                         scale = 'row',
                         main = paste0("Pan Cancer Atlas HNSC - SHH, KRT and EGFR Genes"),
                         fontsize = 10,
                         color = colors,
                         annotation_names_col = FALSE,
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         show_rownames = TRUE,
                         show_colnames = FALSE,
                         treeheight_col = 0,
                         name = "Z-score")
dev.off()


pdf(paste0(output, "Heatmap_SplitPrimarySite_SHH_KRT_EGFR_Genes.pdf"))
ComplexHeatmap::pheatmap(gene_subset_expr,
                         scale = 'row',
                         main = paste0("Pan Cancer Atlas HNSC - SHH, KRT and EGFR Genes"),
                         fontsize = 8,
                         color = colors,
                         column_split = as.factor(annotations$Primary_Site_Group),
                         annotation_col = annotations,
                         annotation_colors = ann_colors,
                         show_rownames = TRUE,
                         show_colnames = FALSE,
                         treeheight_row = 0,
                         name = "Z-score")
dev.off()

pdf(paste0(output, "PAPER_SHH_KRT_EGFR_Genes.pdf"),  width = 3.4, height = 4.5)
ComplexHeatmap::pheatmap(gene_subset_expr,
                         scale = 'row',
                         main = paste0("TCGA HNSC Pan Cancer Atlas"),
                         fontsize = 7,
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


site_col_vector <- ann_colors$Primary_Site_Group

# 2. Create the Legend object
site_legend <- Legend(
  at = names(site_col_vector),
  title = "Site of Primary",
  legend_gp = gpar(fill = site_col_vector),
  direction = "horizontal",
  nrow = 5, # Creates a 2-column layout for 7 items
  title_position = "topleft",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 7)
)

ggsave(
  filename = paste0(output, "TCGA_Legend_Primary_Site.pdf"), 
  plot = site_legend, 
  width = 3, 
  height = 2.5, 
  units = "cm",
  device = cairo_pdf
)

# Z-SCORE LEGEND ---
z_legend <- Legend(
  col_fun = colors, 
  at = c(-4, -2, 0, 2, 4), 
  title = "Z-score", 
  direction = "vertical",
  title_position = "topleft",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 7)
)


ggsave(
  filename = paste0(output, "TCGA_Legend_Zscore.pdf"), 
  plot = z_legend, 
  width = 3, 
  height = 3, 
  units = "cm",
  device = cairo_pdf
)
