##########################################################
# maftools compare PDX and Human (PRH)                   #
# (NMH used as baseline in GATK pipeline to detect germline mutations)   #
# Comparing with TCGA                                    #
# Filtered out HPV/Lymohoma                              #
# PDX VAF < 0.1                                          #
# HUMAN VAF < 0.05                                       #
#                                                        #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 13-02-25                               #
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

########################################################## READ DATA ##########################################################
# Read in MAF file
pdx_maf <- read.maf("./2_pipelines/MAFtools/PDX/Filtered_HPV_Lymphoma_VAF0.1/Filtered_HN0039A_Filtered_PDX_hardfiltered_maftools.maf")
human_maf <- read.maf("./2_pipelines/MAFtools/Human/Filtered_HPV_Lymphoma/Filtered_HN0039A_Filtered_Human_hardfiltered_maftools.maf")
#ClinicalInfo <- read.csv("./0_data/Human/Human_Complete_Curated_Clinical_Information.csv", header = T)

# Extract Case_ID from first 7 characters
human_maf@clinical.data$Case_ID <- substr(human_maf@clinical.data$Tumor_Sample_Barcode, 1, 7)
pdx_maf@clinical.data$Case_ID   <- substr(pdx_maf@clinical.data$Tumor_Sample_Barcode, 1, 7)


########################################################## SET VARIABLES ##########################################################
pipelines <- "./2_pipelines/MAFtools/PDX_Human_Comparison/Filtered_HPV_Lymphoma_VAF0.1/"
output <- "./3_output/MAFtools/PDX_Human_Comparison/Filtered_HPV_Lymphoma_VAF0.1/"
title <- "HN0039A_Filtered_PDX_Human_Comparison"

onco_colours = c(
  "Missense_Mutation" = "#33A02CFF",
  "Nonsense_Mutation" = "#FB9A99FF",
  "Frame_Shift_Del" = "#1F78B4FF",
  "Frame_Shift_Ins" = "#6A3D9AFF" ,
  "Splice_Site" = "#FF7F00FF",
  "In_Frame_Del" = "#FFFF99FF",
  "In_Frame_Ins" = "#B15928FF",
  "Multi_Hit" = "black")

tcga_primary <- c(
  "Larynx" = "#1B9E77FF",
  "Oropharynx" = "#66A61EFF",
  "Hypopharynx" = "#7570B3FF",
  "Oral Cavity" = "#E6AB02FF",
  "Other" = "#666666FF") # RColorBrewer Dark2 pallet

primary_site_colors <- c(
  "Floor of mouth" = "#FFF7BCFF", # Oral Cavity (YlOrBr derived)
  "Gum" = "#FEE391FF",
  "Palate" = "#FE9929FF",
  "Other and unspecified parts of mouth" = "#EC7014FF",
  "Other and unspecified parts of tongue" = "#CC4C02FF",
  
  "Base of tongue" = "#C7E9C0FF", # Oropharynx (RColorBrewer Greens pallet)
  "Tonsil" = "#238B45FF",
  "Oropharynx" = "#66A61EFF",

  "Hypopharynx" = "#7570B3FF",
  
  "Larynx" = "#1B9E77FF",
  
  "Other and ill-defined sites in lip, oral cavity and pharynx" = "#666666FF")
########################################################## DEGs ##########################################################

# Detecting differentially mutated genes
# Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
# Below shows which genes are highly mutated in human compared to pdx
human_vs_pdx <- mafCompare(m1 = human_maf, 
                           m2 = pdx_maf, 
                           m1Name = 'Human', 
                           m2Name = 'PDX', 
                           minMut = 5)
print(human_vs_pdx)

res <- human_vs_pdx$results
write.table(res,paste0(pipelines,title,"Human_vs_PDX.txt"),sep = "\t",row.names = F)

# Forest plot based on the highly mutated genes from the results
pdf(paste0(output,title,"_forestplot.pdf"),width = 7.5, height = 6)
forestPlot(mafCompareRes = human_vs_pdx, pVal = 0.5)
dev.off()


########################################################## Oncoplots ##########################################################
# Add cohort label
human_maf@clinical.data$Cohort <- "Human"
pdx_maf@clinical.data$Cohort <- "PDX"

# Merge MAFs
combined_maf <- merge_mafs(maf = list(human_maf, pdx_maf))

# Plot all genes that appear in either
all_genes <- unique(combined_maf@data$Hugo_Symbol)

oncoplot(maf = combined_maf,
         genes = all_genes,
         clinicalFeatures = "Cohort",
         sortByAnnotation = TRUE,
         removeNonMutated = FALSE,
         titleText = "All Mutations: Human vs PDX")
dev.off()

genes = c("TP53", "CSMD3","FAT1", "KMT2D", "CDKN2A","NOTCH1","FSIP2","PCLO", "RYR2", "ATM", "LRP1B", "CSMD1")

# Plot co oncoplot
pdf(paste0(output,title,"_coOncoplot.pdf"),width = 7.5, height = 6)
coOncoplot(m1 = human_maf, 
           m2 = pdx_maf, 
           m1Name = 'Human', 
           m2Name = 'PDX', 
           genes = genes,
           removeNonMutated = TRUE,
           showSampleNames = F,
           color = onco_colours)
dev.off()

# Plot co bar plot
pdf(paste0(output,title,"_coBarplot.pdf"),width = 7.5, height = 6)
coBarplot(m1 = human_maf,
          m2 = pdx_maf, 
          m1Name = "Human", 
          m2Name = "PDX",
          genes = genes,
          geneSize = 0.75,
          colors = onco_colours)
dev.off()

########################################################## Forest Plots ##########################################################

# New forest plot only containing "genes"
human_vs_pdx_filtered <- human_vs_pdx
human_vs_pdx_filtered$results <- human_vs_pdx_filtered$results %>% 
  filter(Hugo_Symbol %in% genes)

pdf(paste0(output,title,"_forestplot_genelist.pdf"),width = 7.5, height = 6)
forestPlot(mafCompareRes = human_vs_pdx_filtered,
           titleSize = 1,
           pVal = 1)
dev.off()

source("./1_code/RFunctions/forestPlot_pvalue.R")
pdf(paste0(output,title,"_forestplot_genelist_pvalue.pdf"),width = 7.5, height = 6)
forestPlot_pvalue(
  mafCompareRes = human_vs_pdx_filtered,
  titleSize = 1,
  pVal = 1,
  pval_scientific = FALSE,
  pval_digits = 3)
dev.off()


########################################################## Lollipop plots ##########################################################

# Plot lollipop plot FAT1
pdf(paste0(output,title,"_lollipopPlot_TP53.pdf"),width = 7.5, height = 6)
lollipopPlot2(m1 = human_maf,
              m2 = pdx_maf, 
              gene = "TP53", 
              m1_name = "Human", 
              m2_name = "PDX")
dev.off()

# Plot lollipop plot FAT1
pdf(paste0(output,title,"_lollipopPlot_FAT1.pdf"),width = 7.5, height = 6)
lollipopPlot2(m1 = human_maf, 
              m2 = pdx_maf, 
              gene = "FAT1",
              m1_name = "Human", 
              m2_name = "PDX")
dev.off()


########################################################## Compare with TCGA ##########################################################

# Data downloaded directly from PanCancer Atlas on cBioPortal
mut_data <- read.csv("./0_data/hnsc_tcga_pan_can_atlas_2018/HPVRem/mutation_filtered.csv", stringsAsFactors = FALSE)
hnsc_tcga <- read.maf(maf = mut_data, clinicalData = "./0_data/hnsc_tcga_pan_can_atlas_2018/HPVRem/clinicalData_mutation_filtered.csv")

# flags() is an internal function & is not exported. Access it by `:::`
flag_genes = maftools:::flags(top = 20)

plot_comparisons <- function(model, genes_to_plot){
  
  if (model == "Human"){
    outdir <- "Human_vs_TCGA/"
    maf <- human_maf
  } else {
    outdir <- "PDX_vs_TCGA/"
    maf <- pdx_maf
  }
  
  # Detecting differentially mutated genes
  comparison <- mafCompare(m1 = hnsc_tcga, 
                           m2 = maf, 
                           m1Name = 'TCGA HNSC', 
                           m2Name = model, 
                           minMut = 5)
  print(comparison)
  
  res <- comparison$results
  write.table(res, paste0(pipelines, model, "_vs_TCGA_HNSC.txt"), sep = "\t", row.names = F)
  
  # Forest plot
  source("./1_code/RFunctions/forestPlot_pvalue.R")
  pdf(paste0(output, model, "_vs_TCGA_HNSC_forestplot.pdf"), width = 7.5, height = 6)
  forestPlot_pvalue(mafCompareRes = comparison, pVal = 0.05, fdr = 0.05)
  dev.off()
  
  # Filter to include only the genes in genes_to_plot
  comparison_filtered <- comparison
  comparison_filtered$results <- comparison_filtered$results %>% 
    filter(Hugo_Symbol %in% genes)
  
  # Forest plot for filtered genes
  pdf(paste0(output, model, "_vs_TCGA_HNSC_forestplot_genelist.pdf"), width = 7.5, height = 6)
  forestPlot_pvalue(mafCompareRes = comparison_filtered, pVal = 1)
  dev.off()
  
  # Co-oncoplot
  pdf(paste0(output, model, "_TCGA_coOncoplot.pdf"), width = 10, height = 6)
  coOncoplot(m1 = hnsc_tcga, 
             m2 = maf, 
             m1Name = 'TCGA HNSC HPV, Lip', 
             m2Name = model,
             genes = genes, 
             removeNonMutated = TRUE,
             geneNamefont = 0.35, 
             titleFontSize = 1,
             colors = onco_colours)
  dev.off()
  
  # Co-barplot
  pdf(paste0(output, model, "_TCGA_coBarplot.pdf"), width = 7.5, height = 6)
  coBarplot(m1 = hnsc_tcga, 
            m2 = maf, 
            m1Name = "TCGA PanCancer Atlas", 
            m2Name = model,
            genes = genes,
            geneSize = 0.75,
            colors = onco_colours)
  dev.off()
  
  # Loop through genes for lollipop plots
  for (gene in genes_to_plot) {
    lollipop_file <- paste0(output, model, "_TCGA_lollipopPlot_", gene, ".pdf")
    pdf(lollipop_file, width = 7.5, height = 6)
    lollipopPlot2(m1 = hnsc_tcga, 
                  m2 = maf, 
                  gene = gene, 
                  m1_name = "TCGA HNSC", 
                  m2_name = model)
    dev.off()
  }
}

# Remove FSIP2 from list
# Error during wrapup: Structure for protein FSIP2 not found.
genes = c("TP53", "CSMD3","FAT1", "KMT2D", "CDKN2A","NOTCH1","PCLO", "RYR2", "ATM", "LRP1B", "CSMD1", "FSIP2")

# Run function for both models
plot_comparisons(model = "Human", genes_to_plot = genes)
plot_comparisons(model = "PDX", genes_to_plot = genes)


########################################################## TCGA PIE CHART PRIMARY SITE ##########################################################

# Access the clinical data
clinical_data <- hnsc_tcga@clinical.data

# Prepare data
pie_data <- clinical_data %>%
  dplyr::count(cases.primary_site) %>%
  dplyr::mutate(percentage = n / sum(n) * 100,
                label = paste0(round(percentage, 1), "%"),
                legend_label = paste0(cases.primary_site, " (n=", n,")"))

# Custom legend labels
names(pie_data$legend_label) <- pie_data$cases.primary_site
pie_data$legend_label_wrapped <- str_wrap(pie_data$legend_label, width = 35) # Wrap after 40 cahracters, for long labels

total_samples <- sum(pie_data$n)

# Pie chart
p <- ggplot(pie_data, aes(x = "", y = n, fill = cases.primary_site)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = primary_site_colors, labels = pie_data$legend_label_wrapped) +
geom_text(aes(label = ifelse(percentage >= 5, label, "")),  # only label >=5% - skip smaller slices
          position = position_stack(vjust = 0.5), size = 7) +
  labs(title = paste0("Mutational Distribution of Primary Site (n=", total_samples,")"),
       subtitle = "TCGA PanCancer Atlas: HPV+, Lip and Bone/Cartilage Samples Removed",
       fill = "Primary Site") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"), # larger color boxes
        legend.spacing.y = unit(0.4, "cm"), # extra vertical spacing between legend rows
        legend.box.spacing = unit(0.5, "cm")) # space between legend and plot
ggsave(paste0(output, "PieChart_PanCancerAtlas.pdf"), plot = p, width = 10.5, height = 10)



# Grouped using classifications in:
# Johnson, D.E., Burtness, B., Leemans, C.R. et al. Head and neck squamous cell carcinoma. Nat Rev Dis Primers 6, 92 (2020). 
# https://doi.org/10.1038/s41572-020-00224-3
group_map <- list("Oral Cavity" = c("Floor of mouth","Gum","Palate","Other and unspecified parts of mouth","Other and unspecified parts of tongue"),
                  "Oropharynx" = c("Base of tongue","Tonsil","Oropharynx"),
                  "Hypopharynx" = c("Hypopharynx"),
                  "Larynx" = c("Larynx"),
                  "Other" = c("Other and ill-defined sites in lip, oral cavity and pharynx"))


# Map values to custom groups
map_primarysite <- function(x) {
  for (group in names(group_map)) {
    if (x %in% group_map[[group]]) return(group)
  }
  return("Other")
}

# Prepare and group data
pie_data <- clinical_data %>%
  mutate(grouped_site = sapply(cases.primary_site, map_primarysite)) %>%
  dplyr::count(grouped_site) %>%
  mutate(percentage = n / sum(n) * 100,
         label = paste0(round(percentage, 1), "%"),
         legend_label = paste0(grouped_site, " (n=", n,")"))

# Total sample count
total_samples <- sum(pie_data$n)

# Plot pie chart
p <- ggplot(pie_data, aes(x = "", y = n, fill = grouped_site)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = tcga_primary, labels = pie_data$legend_label) +
geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 7) +
  labs(title = paste0("Mutational Grouped Primary Site (n=", total_samples,")"),
       subtitle = "TCGA PanCancer Atlas: HPV+, Lip and Bone/Cartilage Samples Removed",
       fill = "Primary Site") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),        # larger color boxes
        legend.spacing.y = unit(0.4, "cm"),    # extra vertical spacing between legend rows
        legend.box.spacing = unit(0.5, "cm"))    # space between legend and plot
ggsave(paste0(output, "PieChart_PanCancerAtlas_Grouped.pdf"), plot = p, width = 10.5, height = 10)





