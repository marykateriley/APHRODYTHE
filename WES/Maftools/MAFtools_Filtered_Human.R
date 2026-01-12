##########################################################
# maftools Human (PRH)                                   #
# (NMH used as baseline in GATK pipeline to detect germline mutations)   #
# Filtered out HPV/Lymohoma                              #
#                                                        #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 14-09-24                               #
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

########################################################## READ DATA ##########################################################
# Read in MAF file
MAF_file <- read.table("./0_data/Human/AllSamples_Rem_HPV_Lymphoma.maf",header=TRUE,row.names = 1)
ClinicalInfo <- read.csv("./0_data/Human/Human_Rem_HPV_Lymphoma_Clinical_Information.csv", header = T)

# Replace NA values with "No_data" and relevel
ClinicalInfo$Definitive_Response_3_6wk <- ifelse(is.na(ClinicalInfo$Definitive_Response_3_6wk), "No_data",
                                                 ClinicalInfo$Definitive_Response_3_6wk)
ClinicalInfo$Definitive_Response_3_6wk <- factor(ClinicalInfo$Definitive_Response_3_6wk, levels = c("Responder", "NonResponder", "No_data"))
levels(ClinicalInfo$Definitive_Response_3_6wk)


########################################################## SET VARIABLES ##########################################################
data <- "./0_data/Human/"
pipelines <- "./2_pipelines/MAFtools/Human/Filtered_HPV_Lymphoma/"
output <- "./3_output/MAFtools/Human/Filtered_HPV_Lymphoma/"
title <- "HN0039A_Filtered_Human"

onco_colours = c(
  "Missense_Mutation" = "#33A02CFF",
  "Nonsense_Mutation" = "#FB9A99FF",
  "Frame_Shift_Del" = "#1F78B4FF",
  "Frame_Shift_Ins" = "#6A3D9AFF" ,
  "Splice_Site" = "#FF7F00FF",
  "In_Frame_Del" = "#FFFF99FF",
  "In_Frame_Ins" = "#B15928FF",
  "Multi_Hit" = "black")

maf_summary = c(
  "Missense_Mutation" = "#33A02CFF",
  "Nonsense_Mutation" = "#FB9A99FF",
  "Frame_Shift_Del" = "#1F78B4FF",
  "Frame_Shift_Ins" = "#6A3D9AFF" ,
  "Splice_Site" = "#FF7F00FF",
  "In_Frame_Del" = "#FFFF99FF",
  "In_Frame_Ins" = "#B15928FF",
  "Multi_Hit" = "black", # black is not in paired pallet
  "Translation_Start_Site" = "#E31A1CFF", # ok to use red as not shown alongside any green
  "Nonstop_Mutation" = "#A6CEE3FF") 


################################ Adding variant allele freqeuncies to MAF file ################################

# Calculate variant allele frequencies and add as column
MAF_file_with_vaf <- MAF_file
MAF_file_with_vaf$VAF <- (MAF_file_with_vaf$t_alt_count)/(MAF_file_with_vaf$t_alt_count+MAF_file_with_vaf$t_ref_count)

################################## Hard filtering within MAFtools #######################################

# Remove 20 gene flags (**This is different from Shuana's)
MAF_file_filtered <- read.maf(MAF_file_with_vaf, rmFlags = TRUE) 

# Remove introns
maf_nointrons <- subsetMaf(maf=MAF_file_filtered, query="Variant_Type != 'Intron'")

# Filtering MAF file for minor allele frequencies (filter out if not blank, i.e. SNP is not seen in 1000GP)
maf_nomafs <- subsetMaf(maf=maf_nointrons, query="is.na(AF)")

# Test if AF filtered out 
# subsetMaf(maf=maf_nointrons, query="is.na(AF)",mafObj = FALSE)[1:100,"AF"]
# subsetMaf(maf=maf_nointrons, query="is.na(AF)",mafObj = FALSE)[1:100,"Hugo_Symbol"]

# Remove known common variants from FILTER column
maf_filtered <- subsetMaf(maf=maf_nomafs, query="FILTER == 'PASS'")

maf_filtered <- subsetMaf(maf=maf_filtered, query="Consequence != 'synonymous_variant'")

# Remove VAF <= 10%
maf_filtered <- subsetMaf(maf=maf_filtered, query="VAF > 0.05")
maf_filtered 
getSampleSummary(maf_filtered)
getGeneSummary(maf_filtered)


# Save filtered MAF file
write.mafSummary(maf_filtered,basename=paste0(pipelines,"Filtered_",title,"_hardfiltered"))

###################################### Analysing filtered MAF file ########################################
# Plotting MAF summary
# Plots number of variants in each sample as a stacked barplot and variant types as a boxplot summarised by Variant_Classification
pdf(paste0(output,title,"_Filered_VariantSummary.pdf"),width = 7.5, height = 6)
plotmafSummary(maf_filtered,
               rmOutlier=TRUE,
               addStat = "median",
               dashboard = TRUE,
               titvRaw=FALSE,
               color = maf_summary)
dev.off()

# Oncoplot/waterfall plot
pdf(paste0(output,title,"_Filtered_Oncoplot.pdf"),width = 7.5, height = 6)
oncoplot(maf_filtered,top=10)
dev.off()

pdf(paste0(output,title,"_ThesisReady_Filtered_Oncoplot.pdf"),width = 7.5, height = 6)
oncoplot(maf_filtered,top=10, titleText = "Altered in 57 (93.44%) of 61 Human Samples.", # title took from normal oncoplot function, but added in PDX
         drawColBar = FALSE, # Removes top TMB bar chart
         drawRowBar = FALSE,
         colors = onco_colours,) # Removes right % mutation bar chart
dev.off()

# Read in the clinical information
maf_filtered_clin <- read.maf(maf = paste0(pipelines,"Filtered_",title,"_hardfiltered_maftools.maf"), clinicalData = ClinicalInfo)
#clinical_data <- maf_filtered_clin@clinical.data


# Plot lollipop plots for top mutated genes
genes = c("TP53", "CSMD3","FAT1", "KMT2D", "CDKN2A","NOTCH1","FSIP2","PCLO", "RYR2", "ATM", "LRP1B")

for (gene in genes) {
  tryCatch({
    
    pdf(paste0(output,title,"lollipop_",gene,".pdf"), width = 7.5, height = 6)
    lollipopPlot(maf = maf_filtered, gene = gene, AACol = "HGVSp_Short", showMutationRate = TRUE)
    dev.off()
    
  }, error = function(e) {
    message(paste("Skipping gene:", gene, "- Error:", e$message))
    dev.off()
  })  
}

###################################### TNM ########################################
# Mutational Signatures
hnscc.tnm = trinucleotideMatrix(maf = maf_filtered, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
saveRDS(hnscc.tnm$nmf_matrix, file = paste0(pipelines,title,"_TNM_matrix.rds")) 

# Signature analysis
hnscc.sign = estimateSignatures(mat = hnscc.tnm, nTry = 7)
dev.off()
# Draw elbow plot to decide on optimal number of signatures from above (point where y-axis drops significantly)
# Plots as part of estimateSignatures
#plotCophenetic(res = hnscc.sign)

hnscc.sig = extractSignatures(mat = hnscc.tnm, n = 6)
saveRDS(hnscc.sig, file = paste0(pipelines,title,"_SignatureMatrix.rds")) 

# Compare against COSMIC legacy signature database
# compareSignatures returns full table of cosine similarities against COSMIC signatures, which can be further analysed.
hnscc.og30.cosm = compareSignatures(nmfRes = hnscc.sig, sig_db = "legacy")

# Compare against COSMIC legacy signature database
hnscc.v3.cosm = compareSignatures(nmfRes = hnscc.sig, sig_db = "SBS")
saveRDS(hnscc.v3.cosm, file = paste0(pipelines,title,"_cosmic_match.rds"))

# Plot signature similarities
pdf(paste0(output,title,"_COSMIClegacy_sig_similarity.pdf"),width = 7.5, height = 6)
pheatmap::pheatmap(mat = hnscc.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
dev.off()

# Plot signature similarities
pdf(paste0(output,title,"_COSMIC_SBS_sig_similarity.pdf"),width = 7.5, height = 6)
pheatmap::pheatmap(mat = hnscc.v3.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
dev.off()

# Plot signatures
pdf(paste0(output,title,"_COSMIC_SBS_signature.pdf"),width = 7.5, height = 6)
plotSignatures(nmfRes = hnscc.sig, title_size = 0.8, sig_db = "SBS")
dev.off()

pdf(paste0(output,title,"_COSMIClegacy_signature.pdf"),width = 7.5, height = 6)
plotSignatures(nmfRes = hnscc.sig, title_size = 0.8, sig_db = "legacy")
dev.off()



########################################################## TRINUCLEOTIDE FREUQNCY ##########################################################
# Define a color palette for mutation contexts
color_palette <- c(
  "C>A" = "#954535",   
  "C>G" = "grey",  
  "C>T" = "#F33A6A", 
  "T>A" = "#FF7F50",
  "T>C" = "darkgreen",
  "T>G" = "cornflowerblue")



## HUMAN ##
# Make matrix
tnm_counts <- hnscc.tnm$nmf_matrix  

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
  labs(title = "Human Average Trinucleotide Mutational Frequencies",
       x = "Trinucleotide",
       y = "Average Frequency") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  # Adjust text size here
    axis.title.x = element_text(size = 10),  # Adjust x-axis title size
  ) +
  coord_cartesian(clip = 'off')  # Optional: prevent clipping of axis text
ggsave(paste0(output,title,"_AvgTrinucleotideFreq.pdf"), plot = p_resp, width = 10, height = 7)





# Create subset of maf file with responders only and with non-responders only

# Read in the clinical information
maf_filtered_clin <- read.maf(maf = paste0(pipelines,"Filtered_",title,"_hardfiltered_maftools.maf"), clinicalData = ClinicalInfo)
maf_responders <- subsetMaf(maf=maf_filtered_clin, clinQuery="Definitive_Response_3_6wk == 'Responder'")
write.mafSummary(maf_responders,basename=paste0(pipelines,title,"_hardfiltered_Responders"))

maf_nonresponders <- subsetMaf(maf=maf_filtered_clin, clinQuery="Definitive_Response_3_6wk == 'NonResponder'")
write.mafSummary(maf_nonresponders,basename=paste0(pipelines,title,"_hardfiltered_NonResponders"))

