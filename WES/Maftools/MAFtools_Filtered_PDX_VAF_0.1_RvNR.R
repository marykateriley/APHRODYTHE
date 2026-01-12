##########################################################
# maftools compare PDX Responders and NonResponders      #
#                                                        #
# Filtered out HPV/Lymohoma                              #
# VAF < 0.1                                              #
#                                                        #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 13-02-25                               #
##########################################################

setwd("")

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
library(patchwork) # combine plots
library(reshape2) # melt
library(RColorBrewer)
library(ggpubr)
library(viridis)
library(MutationalPatterns)

########################################################## SET VARIABLES ##########################################################

data <- "./0_data/PDX/"
pipelines <- "./2_pipelines/MAFtools/PDX_RvNR/Filtered_HPV_Lymphoma_VAF0.1/"
output <- "./3_output/MAFtools/PDX_RvNR/Filtered_HPV_Lymphoma_VAF0.1/"
title <- "HN0039A_Filtered_PDX_RvNR_"

onco_colours = c(
  "Missense_Mutation" = "#33A02CFF",
  "Nonsense_Mutation" = "#FB9A99FF",
  "Frame_Shift_Del" = "#1F78B4FF",
  "Frame_Shift_Ins" = "#6A3D9AFF" ,
  "Splice_Site" = "#FF7F00FF",
  "In_Frame_Del" = "#FFFF99FF",
  "In_Frame_Ins" = "#B15928FF",
  "Multi_Hit" = "black")

########################################################## READ DATA ##########################################################

ClinicalInfo <- read.csv(paste0(data,"PDX_Complete_Curated_Clinical_Information.csv"), header = T)

# Replace NA values with "No_data" and relevel
ClinicalInfo$Definitive_Response_3_6wk <- ifelse(is.na(ClinicalInfo$Definitive_Response_3_6wk), "No_data",
                                                 ClinicalInfo$Definitive_Response_3_6wk)
ClinicalInfo <- ClinicalInfo[ClinicalInfo$Definitive_Response_3_6wk != "No_data", ] # Remove samples that dont have reponse
ClinicalInfo <- ClinicalInfo[ClinicalInfo$Remove != "Yes", ] # Remove HPV/Lymphoma samples
ClinicalInfo$Definitive_Response_3_6wk <- factor(ClinicalInfo$Definitive_Response_3_6wk, levels = c("Responder", "NonResponder"))
levels(ClinicalInfo$Definitive_Response_3_6wk)

# Convert to numeric (from character/factor)
ClinicalInfo$percentage_vol_change_3_6wk <- as.character(ClinicalInfo$percentage_vol_change_3_6wk)
ClinicalInfo$percentage_vol_change_numeric <- as.numeric(ClinicalInfo$percentage_vol_change_3_6wk)

# Reorder dataframe from high to low
ClinicalInfo <- ClinicalInfo[order(-ClinicalInfo$percentage_vol_change_numeric), ]

# Reset factor levels in the same high-to-low order
ordered_levels <- unique(ClinicalInfo$percentage_vol_change_3_6wk)
ClinicalInfo$percentage_vol_change_3_6wk <- factor(ClinicalInfo$percentage_vol_change_3_6wk, levels = ordered_levels)
ClinicalInfo$percentage_vol_change_numeric <- NULL

ClinicalInfo$Site.of.Primary[ClinicalInfo$Site.of.Primary == "" | is.na(ClinicalInfo$Site.of.Primary)] <- "Unknown"
ClinicalInfo$Site.of.Primary <- factor(ClinicalInfo$Site.of.Primary, levels = c("Larynx","Oropharynx","Hypopharynx","Oral Cavity", "Larynx + Oral Cavity", "Oral Cavity + Oropharynx", "Unknown"))

ClinicalInfo$T_Stage[ClinicalInfo$T_Stage == "" | is.na(ClinicalInfo$T_Stage)] <- "No_data"
ClinicalInfo$T_Stage <- factor(ClinicalInfo$T_Stage, levels = c("T2", "T3", "T4", "T4a", "No_data"))


ClinicalInfo$N_Stage[ClinicalInfo$N_Stage == "" | is.na(ClinicalInfo$N_Stage)] <- "No_data"
ClinicalInfo$N_Stage <- factor(ClinicalInfo$N_Stage, levels = c("N0", "N1", "N2a", "N2b", "N2c", "N3", "N3b", "No_data"))

ClinicalInfo$HISTOLOGICAL.GRADE[ClinicalInfo$HISTOLOGICAL.GRADE == "" | is.na(ClinicalInfo$HISTOLOGICAL.GRADE)] <- "No_data"
ClinicalInfo$HISTOLOGICAL.GRADE <- factor(ClinicalInfo$HISTOLOGICAL.GRADE, 
                                          levels = c("g1", "g1-g2", "g2", "g2-g3", "g3", "No_data"))

ClinicalInfo$Sex[ClinicalInfo$Sex == "" | is.na(ClinicalInfo$Sex)] <- "No_data"
ClinicalInfo$Sex <- factor(ClinicalInfo$Sex, 
                           levels = c("M", "F", "Unknown"))

ClinicalInfo$Age_Group <- factor(ClinicalInfo$Age_Group, levels = c("50-59", "60-69", "70-79", "80-89", "90-99"))
levels(ClinicalInfo$Age_Group)


# Read in maf file
maf_filtered_clin <- read.maf(maf = "./2_pipelines/MAFtools/PDX/Filtered_HPV_Lymphoma_VAF0.1/Filtered_HN0039A_Filtered_PDX_hardfiltered_maftools.maf", clinicalData = ClinicalInfo)

# Subset into responders
resp_maf <- subsetMaf(maf=maf_filtered_clin, clinQuery="Definitive_Response_3_6wk == 'Responder'")
write.mafSummary(resp_maf,basename=paste0(pipelines,title,"hardfiltered_Responders"))

# Subset into nonresponders
nonResp_maf <- subsetMaf(maf=maf_filtered_clin, clinQuery="Definitive_Response_3_6wk == 'NonResponder'")
write.mafSummary(nonResp_maf,basename=paste0(pipelines,title,"hardfiltered_NonResponders"))




########################################################## DEGs ##########################################################

# Detecting differentially mutated genes
# Considering only genes which are mutated in at-least in 3
# 5 samples is better to avoid bias due to genes mutated in single sample.
# But we're comparing 13 NR to 32 R - smaller cohort
# Below shows which genes are highly mutated in responder v nonresponder
resp_vs_nonResp <- mafCompare(m1 = resp_maf, 
                              m2 = nonResp_maf, 
                              m1Name = 'Responders', 
                              m2Name = 'NonResponders',
                              minMut = 3)

print(resp_vs_nonResp)

res <- resp_vs_nonResp$results
write.table(res,paste0(pipelines,title,"Responders_vs_NonResponders.txt"),sep = "\t",row.names = F)
########################################################## FOREST PLOTS ##########################################################

# Forest plot based on the highly mutated genes from the results
pdf(paste0(output,title,"_forestplot.pdf"),width = 7.5, height = 6)
forestPlot(mafCompareRes = resp_vs_nonResp, 
           titleSize = 1,
           pVal = 0.1)
dev.off()

# Significant genes in forest plot
forest_genelist <- c("DCHS1", "FAM135B","PCDHA6","POLD1", "SIPA1L3", "CSMD1")


########################################################## COONCOPLOTS ##########################################################

# Plot co oncoplot
pdf(paste0(output,title,"_coOncoplot.pdf"),width = 7.5, height = 6)
coOncoplot(m1 = resp_maf, 
           m2 = nonResp_maf,
           m1Name = 'Responders', 
           m2Name = 'NonResponders', 
           geneNamefont = 0.7,
           gene_mar = 3,
           legend_height = 4,
           titleFontSize = 1.05,
           removeNonMutated = TRUE,
           showSampleNames = FALSE)
dev.off()

# Plot co bar plot
pdf(paste0(output,title,"_coBarplot_genelist.pdf"),width = 7.5, height = 6)
coBarplot(m1 = resp_maf, 
          m2 = nonResp_maf, 
          m1Name = "Responders", 
          m2Name = "NonResponders")
dev.off()



genes = c("TP53", "CSMD3","FAT1", "KMT2D", "NOTCH1", "PCLO", "CDKN2A", "FSIP2", "ATM", "RYR2")

# Plot co oncoplot
pdf(paste0(output,title,"_coOncoplot_genelist.pdf"),width = 7.5, height = 6)
coOncoplot(m1 = resp_maf, 
           m2 = nonResp_maf,
           m1Name = 'Responders', 
           m2Name = 'NonResponders', 
           genes = genes, 
           geneNamefont = 0.7,
           gene_mar = 3,
           legend_height = 4,
           titleFontSize = 1.32,
           removeNonMutated = TRUE,
           showSampleNames = FALSE)
dev.off()

# Plot co bar plot
pdf(paste0(output,title,"_coBarplot_genelist.pdf"),width = 7.5, height = 6)
coBarplot(m1 = resp_maf, 
          m2 = nonResp_maf, 
          m1Name = "Responders", 
          m2Name = "NonResponders",
          genes = genes)
dev.off()

# New forest plot only containing "genes"
resp_vs_nonResp_filtered <- resp_vs_nonResp
resp_vs_nonResp_filtered$results <- resp_vs_nonResp_filtered$results %>% 
  filter(Hugo_Symbol %in% genes)

pdf(paste0(output,title,"_forestplot_genelist.pdf"),width = 7.5, height = 6)
forestPlot(mafCompareRes = resp_vs_nonResp_filtered,
           titleSize = 1,
           pVal = 1)
dev.off()

source("./1_code/RFunctions/forestPlot_pvalue.R")
pdf(paste0(output,title,"_forestplot_genelist_pvalue.pdf"),width = 7.5, height = 6)
forestPlot_pvalue(
  mafCompareRes = resp_vs_nonResp_filtered,
  titleSize = 1,
  pVal = 1,
  pval_scientific = FALSE,
  pval_digits = 3
)
dev.off()

########################################################## COONCOPLOTS + CLINICAL DATA ##########################################################

## SET UP COLOURS

# Define PRGn-anchored continuous mapper
pal <- rev(brewer.pal(11, "PRGn"))
col_fun <- circlize::colorRamp2(
  c(-89.98, 0, 334.80), # Colour breaks
  c(pal[2], pal[6], pal[10]))

# Build a named colour vector for the factor levels
levs <- levels(ClinicalInfo$percentage_vol_change_3_6wk)
lev_num <- as.numeric(as.character(levs)) # preserves decimals

# order levels numerically
o <- order(lev_num)
levs <- levs[o]
lev_num <- lev_num[o]

vol_colors_vec <- col_fun(lev_num)
names(vol_colors_vec) <- levs

# Create list
vol_colors <- list(percentage_vol_change_3_6wk = vol_colors_vec)


site_colours = c("#1B9E77FF", "#66A61EFF", "#7570B3FF","#E6AB02FF", "#D95F02FF", "#E7298AFF", "#666666FF")
names(site_colours) <- c("Larynx", "Oropharynx", "Hypopharynx", "Oral Cavity", 
                        "Larynx + Oral Cavity","Oral Cavity + Oropharynx", "Unknown") 
site_colours <- list(Site.of.Primary=site_colours)

sex_colours = c("green", "violet", "black")
names(sex_colours) <- c("M", "F", "Unknown")
sex_colours <- list(Sex=sex_colours)

t_colours <- c("#1abc9c", "#f5b041", "#95dc85", "#4158d1", "black")
names(t_colours) <- c("T2", "T3", "T4", "T4a", "No_data")
t_colours <- list(T_Stage = t_colours)

n_colours <- c("#FF5733", "#d141c5", "#ec1549",  "#dfdf7f", "#0a9f70", "#126f81", "#ef9219", "black")
names(n_colours) <- c("N0", "N1", "N2a", "N2b", "N2c", "N3", "N3b", "No_data")
n_colours <- list(N_Stage = n_colours)

grade_colours <- c("#f8b3b3", "#f8e3b3", "#c6f8b3",  "#b3e1f8", "#ebb3f8", "black")
names(grade_colours) <- c("g1", "g1-g2", "g2", "g2-g3", "g3", "No_data")
grade_colours <- list(HISTOLOGICAL.GRADE = grade_colours)

age_colours <- c("#160d7f", "#0d7f65", "#777f0d",  "#7f0d75", "#7f0d0d")
names(age_colours) <- c("50-59", "60-69", "70-79", "80-89", "90-99")
age_colours <- list(Age_Group = age_colours)


#### S0PT SAMPLES ####

# Convert to numeric for sorting
resp_clin <- resp_maf@clinical.data
nonresp_clin <- nonResp_maf@clinical.data

resp_clin$perc_num <- as.numeric(as.character(resp_clin$percentage_vol_change_3_6wk))
nonresp_clin$perc_num <- as.numeric(as.character(nonresp_clin$percentage_vol_change_3_6wk))

# Sort from low to high percentage volume change
sorted_samples_resp <- resp_clin$Tumor_Sample_Barcode[order(resp_clin$perc_num)]
sorted_samples_nonResp <- nonresp_clin$Tumor_Sample_Barcode[order(nonresp_clin$perc_num)]

########## update the colour to work with maftools

source("./1_code/RFunctions/CoOncoplot_AnnRem.R")
pdf(paste0(output,title,"_CoOncoplot_ThesisReady.pdf"), width = 17, height = 11)
coOncoplot(m1 = resp_maf, 
           m2 = nonResp_maf,
           m1Name = 'Responders', 
           m2Name = 'NonResponders', 
           legend_height = 0.1,
           legendFontSize = 0.1,
           annotationFontSize = 1,
           geneNamefont = 1.2,
           outer_mar = 3.5,
           titleFontSize = 2,
           removeNonMutated = FALSE,
           sortByAnnotation1 = TRUE, 
           sortByAnnotation2 = TRUE, 
           clinicalFeatures1 = c('percentage_vol_change_3_6wk', 'Site.of.Primary'),
           clinicalFeatures2 = c('percentage_vol_change_3_6wk', 'Site.of.Primary'),
           annotationColor1 = c(vol_colors, site_colours),
           annotationColor2 = c(vol_colors, site_colours),
           sampleOrder1 = sorted_samples_resp,
           sampleOrder2 = sorted_samples_nonResp,
           colors = onco_colours)
dev.off()

pdf(paste0(output,title,"_CoOncoplot_ThesisReady_genelist.pdf"), width = 16.5, height = 11)
coOncoplot(m1 = resp_maf, 
           m2 = nonResp_maf,
           m1Name = 'Responders', 
           m2Name = 'NonResponders', 
           genes = genes,
           legend_height = 0.1,
           legendFontSize = 0.1,
           annotationFontSize = 0.0001,
           geneNamefont = 1.2,
           outer_mar = 3.5,
           titleFontSize = 2,
           removeNonMutated = FALSE,
           sortByAnnotation1 = TRUE, 
           sortByAnnotation2 = TRUE, 
           clinicalFeatures1 = c('percentage_vol_change_3_6wk', 'Site.of.Primary'),
           clinicalFeatures2 = c('percentage_vol_change_3_6wk', 'Site.of.Primary'),
           annotationColor1 = c(vol_colors, site_colours),
           annotationColor2 = c(vol_colors, site_colours),
           sampleOrder1 = sorted_samples_resp,
           sampleOrder2 = sorted_samples_nonResp,
           colors = onco_colours)
dev.off()



########################################################## LEGENDS ##########################################################

### Primary Site ###
site_col_vector <- site_colours$Site.of.Primary
site_legend <- Legend(
  at = names(site_col_vector),
  title = "Site of Primary",
  legend_gp = gpar(fill = site_col_vector),
  direction = "horizontal",
  nrow = 4,
  title_position = "topleft",
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 10))
pdf(paste0(output, "Site_of_Primary_Legend.pdf"), width = 3.5, height = 4.5)
draw(site_legend)
dev.off()


### Percentage Volume Change ###
rng <- range(lev_num, na.rm = TRUE)
lo <- min(rng[1])
hi <- max(rng[2])

# Ticks at cut-offs
ticks <- c(lo, 0, hi)
labels <- c(paste0(round(lo, 2), "%"),
            "0%",
            paste0(round(hi, 2), "%"))

leg <- Legend(
  col_fun = col_fun,
  at = ticks,
  labels = labels,
  title = "Percentage Volume Change 3-6 wk",
  direction = "horizontal",
  title_position = "topleft",
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 10),
  grid_width = unit(3.4, "in"),
  grid_height = unit(0.28, "in"))
pdf(file = paste0(output, "percentage_vol_change_3_6wk_legend.pdf"), width = 4.6, height = 1.6)
draw(leg)
dev.off()


## Mutational ##
onco_legend <- Legend(
  at = names(onco_colours),
  title = NULL,
  legend_gp = gpar(fill = unname(onco_colours), col = NA),
  direction = "horizontal",
  nrow = 4,
  labels_gp = gpar(fontsize = 10),
  grid_width  = unit(4, "mm"),
  grid_height = unit(4, "mm"),
  gap = unit(1.5, "mm"))
pdf(paste0(output, "OncoLegend.pdf"), width = 4.5, height = 1.8)
draw(onco_legend)
dev.off()


########################################################## MUTATIONAL SIGNATURE ANALYSIS ##########################################################
# Mutational Signatures - effects of carcinogens, DNA repair defects etc.

## RESPONDERS ##

# Creates a trinucleotide matrix that considers each mutation in its surrounding nucleotide environment (e.g A[C>T]G)
resp_tnm = trinucleotideMatrix(maf = resp_maf, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
saveRDS(resp_tnm$nmf_matrix, file = paste0(pipelines,title,"_Responders_TNM_matrix.rds")) 

# Estimate no. of mutational signatures
resp_sign = estimateSignatures(mat = resp_tnm, nTry = 7) 
dev.off()

# Draw elbow plot to decide on optimal number of signatures from above (point where y-axis drops significantly)
# Plots as part of estimateSignatures
plotCophenetic(res = resp_sign)

# Extract signatures using information from the elbow plot
# Uses non-negative matrix factorization to decompose the matrix into n signatures
resp.sig = extractSignatures(mat = resp_tnm, n = 4)
saveRDS(resp.sig, file = paste0(pipelines,title,"_Responders_SignatureMatrix.rds"))  # or "Human_signature_matrix.rds"


# Compare against COSMIC SBS signature database
resp.v3.cosm = compareSignatures(nmfRes = resp.sig, sig_db = "SBS")
saveRDS(resp.v3.cosm, file = paste0(pipelines,title,"_Responders_cosmic_match.rds"))  # or "Human_signature_matrix.rds"



pdf(paste0(output,title,"_COSMIC_SBS_signature_responders.pdf"),width = 7.5, height = 6)
resp_plot <- plotSignatures(nmfRes = resp.sig, title_size = 0.8, sig_db = "SBS")# single base substitutions
#mtext("SBS Signatures Responders", side = 3, line = 0.5, cex = 1.2)
dev.off()

pdf(paste0(output,title,"_COSMIClegacy_signature_responders.pdf"),width = 7.5, height = 6)
plotSignatures(nmfRes = resp.sig, title_size = 0.8, sig_db = "legacy")
dev.off()




## NON-RESPONDERS ##


# Creates a trinucleotide matrix that considers each mutation in its surrounding nucleotide environment (e.g A[C>T]G)
nonResp_tnm = trinucleotideMatrix(maf = nonResp_maf, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
saveRDS(nonResp_tnm$nmf_matrix, file = paste0(pipelines,title,"_NonResponders_TNM_matrix.rds")) 

# Estimate no. of mutational signatures
nonResp_sign = estimateSignatures(mat = nonResp_tnm, nTry = 7) 
dev.off()

# Draw elbow plot to decide on optimal number of signatures from above (point where y-axis drops significantly)
# Plots as part of estimateSignatures
plotCophenetic(res = nonResp_sign)

# Extract signatures using information from the elbow plot
# Uses non-negative matrix factorization to decompose the matrix into n signatures
nonResp.sig = extractSignatures(mat = nonResp_tnm, n = 3)
saveRDS(nonResp.sig, file = paste0(pipelines,title,"_NonResponders_SignatureMatrix.rds"))  # or "Human_signature_matrix.rds"


# Compare against COSMIC SBS signature database
nonresp.v3.cosm = compareSignatures(nmfRes = nonResp.sig, sig_db = "SBS")
saveRDS(nonresp.v3.cosm, file = paste0(pipelines,title,"_NonResponders_cosmic_match.rds"))  # or "Human_signature_matrix.rds"


pdf(paste0(output,title,"_COSMIC_SBS_signature_nonresponders.pdf"),width = 7.5, height = 6)
non_resp_plot <- plotSignatures(nmfRes = nonResp.sig, title_size = 0.8, sig_db = "SBS") # single base substitutions
dev.off()

pdf(paste0(output,title,"_COSMIClegacy_signature_nonresponders.pdf"),width = 7.5, height = 6)
plotSignatures(nmfRes = nonResp.sig, title_size = 0.8, sig_db = "legacy")
dev.off()




# Load the saved trinucleotide mutation matrix for full cohort
full_tnm_matrix <- readRDS("./2_pipelines/MAFtools/PDX/Filtered_HPV_Lymphoma_VAF0.1/HN0039A_Filtered_PDX_TNM_matrix.rds")

full_tnm_matrix = trinucleotideMatrix(maf = maf_filtered_clin, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")




########################################################## TRINUCLEOTIDE FREUQNCY ##########################################################
# Define a color palette for mutation contexts
color_palette <- c(
  "C>A" = "#954535",   
  "C>G" = "grey",  
  "C>T" = "#F33A6A", 
  "T>A" = "#FF7F50",
  "T>C" = "darkgreen",
  "T>G" = "cornflowerblue")

## RESPONDERS ##
# Make matrix
tnm_counts <- resp_tnm$nmf_matrix  

#Calculate the average trinucleotide counts
avg_tnm <- colMeans(tnm_counts, na.rm = TRUE)  # Use na.rm = TRUE to ignore NAs

# Normalize to get frequencies
total_count <- sum(avg_tnm)
avg_frequencies <- avg_tnm / total_count

# Create a data frame for plotting
tnm_df <- data.frame(Trinucleotide = names(avg_frequencies), Frequency = avg_frequencies)

# Extract the mutation context C>A, C>G, etc. from the trinucleotide string
tnm_df$MutationContext <- gsub(".*?\\[(.*?)\\].*", "\\1", tnm_df$Trinucleotide)  # Extract mutation context without brackets


# Print the MutationContext values to check
print(tnm_df$MutationContext)

# Map colors based on mutation contexts
tnm_df$Colour <- color_palette[tnm_df$MutationContext]  # Directly use the color_palette mapping

# Ensure any unmatched MutationContext values are set to "grey"
tnm_df$Colour[is.na(tnm_df$Colour)] <- "grey"

# Convert MutationContext to a factor to maintain order
tnm_df$MutationContext <- factor(tnm_df$MutationContext, levels = unique(tnm_df$MutationContext))

# Reorder Trinucleotide based on the levels of MutationContext
tnm_df$Trinucleotide <- factor(tnm_df$Trinucleotide, levels = tnm_df$Trinucleotide[order(tnm_df$MutationContext)])

# Plot the average trinucleotide frequencies with colors reflecting MutationContext
p_resp <- ggplot(tnm_df, aes(x = Trinucleotide, y = Frequency, fill = MutationContext)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use position = "dodge" for side-by-side bars
  scale_fill_manual(values = color_palette, name = "Mutation Context") +  # Use manual color mapping with a legend
  theme_bw() +
  labs(title = "PDX Average Trinucleotide Mutational Frequencies: Responders",
       x = "Trinucleotide",
       y = "Average Frequency") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Adjust text size here
    axis.title.x = element_text(size = 10),  # Adjust x-axis title size
  ) +
  coord_cartesian(clip = 'off')  # Optional: prevent clipping of axis text
ggsave(paste0(output,title,"_AvgTrinucleotideFreq_responders.pdf"), plot = p_resp, width = 10, height = 7)




## NONRESPONDERS ##
# Make matrixs
tnm_counts <-nonResp_tnm$nmf_matrix  
avg_tnm <- colMeans(tnm_counts, na.rm = TRUE)  # Use na.rm = TRUE to ignore NAs

# Normalize to get frequencies
total_count <- sum(avg_tnm)
avg_frequencies <- avg_tnm / total_count

# Create a data frame for plotting
tnm_df <- data.frame(Trinucleotide = names(avg_frequencies), Frequency = avg_frequencies)

# Extract the mutation context C>A, C>G, etc. from the trinucleotide string
tnm_df$MutationContext <- gsub(".*?\\[(.*?)\\].*", "\\1", tnm_df$Trinucleotide)  # Extract mutation context without brackets

# Print the MutationContext values to check
print(tnm_df$MutationContext)

# Map colors based on mutation contexts
tnm_df$Colour <- color_palette[tnm_df$MutationContext]  # Directly use the color_palette mapping

# Ensure any unmatched MutationContext values are set to "grey"
tnm_df$Colour[is.na(tnm_df$Colour)] <- "grey"

# Convert MutationContext to a factor to maintain order
tnm_df$MutationContext <- factor(tnm_df$MutationContext, levels = unique(tnm_df$MutationContext))

# Reorder Trinucleotide based on the levels of MutationContext
tnm_df$Trinucleotide <- factor(tnm_df$Trinucleotide, levels = tnm_df$Trinucleotide[order(tnm_df$MutationContext)])

# Plot the average trinucleotide frequencies with colors reflecting MutationContext
p_nonResp <- ggplot(tnm_df, aes(x = Trinucleotide, y = Frequency, fill = MutationContext)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use position = "dodge" for side-by-side bars
  scale_fill_manual(values = color_palette, name = "Mutation Context") +  # Use manual color mapping with a legend
  theme_bw() +
  labs(title = "PDX Average Trinucleotide Mutational Frequencies: Non-Responders",
       x = "Trinucleotide",
       y = "Average Frequency") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Adjust text size here
    axis.title.x = element_text(size = 10),  # Adjust x-axis title size
  ) +
  coord_cartesian(clip = 'off')  # Optional: prevent clipping of axis text
ggsave(paste0(output,title,"_AvgTrinucleotideFreq_nonresponders.pdf"), plot = p_nonResp, width = 10, height = 7)


## COMBINE ##
# Combine using patchwork
combined_plot <- p_resp / p_nonResp +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom') +
  labs(x = "Trinucleotide", y = "Average Frequency")
# Save the combined plot
ggsave(paste0(output,title,"_AvgTrinucleotideFreq_combined.pdf"), plot = combined_plot, width = 10, height = 8)

