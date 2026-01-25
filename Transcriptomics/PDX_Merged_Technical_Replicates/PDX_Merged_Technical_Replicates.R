################################################################
# Removing Technical Replicates                                #
#                                                              #
# PN0039 - Basal HNSCC                                         #
#                                                              #
# Mary-Kate Riley 24-07-24                                     #
################################################################

setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(dplyr)
library(DESeq2)

########################################################## READ DATA ##########################################################
data <- "./0_data/DESeq2/PRX/"

# Full merged htcounts
GEX <- read.csv("./0_data/Merged_Counts/Unfiltered/PDX_Graft/HNC_PRX_graft_raw_counts_matrix_080124.csv", header = T, row.names = 1)

# Subsetted GEX of only technical replicates
tech_GEX <- read.csv("./2_pipelines/Replicates/HNSCC_technical_replicate_matrix_10012024.csv", header = T, row.names = 1)

# Load the clinical information metadata data for all cases
metadata <- read.csv("./0_data/Complete_Curated_Clinical_Information.csv", header = TRUE)

# Load the sample data representing all sample batches, identifiers and timepoints
sample_names <- read.csv("./0_data/Names_Sample_Sheet_PRX_Graft.csv", header = TRUE)

# One Sample_ID originally has a "-" for a technical replicate, but this sample will be changed to a . in GEX (column names cannot contain -)
# Change the name now so it matches with the GEX
grep("-", sample_names$Sample_ID) # Checks what row has a "-" in Sample_ID
sample_names$Sample_ID <- gsub("-", ".", sample_names$Sample_ID) # Replaces "-" with a "."



########################################################## colData ##########################################################

########################################################## DATA CURATION colData ##########################################################
## Full colData ##
colData <- merge(sample_names, metadata, by = "Case_ID", all.x = TRUE)
row.names(colData) <- colData$Sample_ID
write.csv(colData, file = paste0(data, "Metadata/Complete_GraftHNSCC_colData_240724.csv"), quote=F, row.names = TRUE)
length(unique(colData$Case_ID))


########################################################## DATA CURATION colData MERGE TECHNICAL ##########################################################

## Merge Technical ##
tech_names <- colnames(tech_GEX)
tech_names_prefix <- substr(tech_names, 1, 7) # Extract the Case_ID: HNC000

# Subset the metadata to include only rows where Sample_ID is in tech_colnames
merged_metadata <- colData[colData$Sample_ID %in% tech_names, ]

# Group by Case_ID and summarize metadata by taking the first occurrence for each group
merged_metadata <- merged_metadata %>%
  group_by(Case_ID) %>%
  summarise(across(everything(), ~ .[1])) %>%
  as.data.frame()
rownames(merged_metadata) <- merged_metadata$Sample_ID

# Subset the original metadata to get rows that were not part of the merging process
remaining_metadata <- colData[!(colData$Sample_ID %in% tech_names), ]

# Combine the merged metadata with the remaining biological samples
colData_tech <- bind_rows(merged_metadata, remaining_metadata)
write.csv(colData_tech, file = paste0(data, "Metadata/TechMerged_Unfiltered_GraftHNSCC_colData_240724.csv"), quote=F, row.names = TRUE)
length(unique(colData_tech$Case_ID))



########################################################## GEX ##########################################################

########################################################## DATA CURATION GEX ##########################################################

##### GEX_raw (all counts) ##### 
GEX_raw <- GEX %>%
  select(all_of(colData$Sample_ID)) # only keep Sample_IDs present in graft GEX
write.csv(GEX_raw, file = paste0(data, "GEX/RawCounts_GraftHNSCC_270424.csv"), quote=F, row.names = TRUE)

case_id <- unique(colData$Case_ID) # How many Case_Ids remaining
print(length(case_id))



##### GEX_norm (all counts) ##### 
dds <- DESeqDataSetFromMatrix(countData = GEX_raw, 
                              colData = colData, 
                              design = ~ 1)
ddsBlind <- DESeq(dds)

# Normalised count matrix
GEX_norm <- counts(ddsBlind, normalized = TRUE)
write.csv(GEX_norm, file = paste0(data, "GEX/NormCounts_GraftHNSCC_270424.csv"), quote=F, row.names = TRUE)





########################################################## DATA CURATION GEX MERGE TECHNICAL ##########################################################
# Create a data frame with prefixes
merge_df <- data.frame(
  Original_Name = tech_names,
  Prefix = tech_names_prefix,
  stringsAsFactors = FALSE)

# Get the first occurrence of each Prefix
first_occurrence <- merge_df %>%
  distinct(Prefix, .keep_all = TRUE) %>%
  pull(Original_Name)

#Create a list to store merged columns
merged_columns <- list()

# Loop through each unique prefix and average the columns
for (prefix in unique(tech_names_prefix)) {
  cols_to_merge <- tech_names[tech_names_prefix == prefix]
  merged_columns[[prefix]] <- rowMeans(tech_GEX[, cols_to_merge, drop = FALSE])
}

# Combine the merged columns into a new matrix
merged_GEX_subset <- do.call(cbind, merged_columns)

# Ensure the column names match the first occurrence names
# Map each prefix to the first occurrence name
prefix_to_name <- setNames(first_occurrence, unique(tech_names_prefix))

# Set the column names for the merged matrix
colnames(merged_GEX_subset) <- prefix_to_name[names(merged_columns)]

merged_GEX_subset <- as.data.frame(merged_GEX_subset)
merged_GEX_subset <- round(merged_GEX_subset)

# Extract columns from GEX that are not technical replicates
non_tech_columns <- setdiff(colnames(GEX), tech_names)

# Extract the remaining part of GEX that is not part of technical replicates
remaining_GEX <- GEX[, non_tech_columns]

# Combine the remaining columns with the merged technical replicates
GEX_tech <- cbind(remaining_GEX, merged_GEX_subset)


##### GEX_tech_rem ##### 
GEX_tech_rem <- GEX_tech %>%
  select(all_of(colData_tech$Sample_ID))  # only keep Sample_IDs present in graft GEX
write.csv(GEX_tech_rem, file = paste0(data, "GEX/TechMerged_Unfiltered_RawCounts_GraftHNSCC_270424.csv"), quote=F, row.names = TRUE)

rem_samples <- setdiff(colData$Sample_ID, colData_tech$Sample_ID) # What samples were removed 
print(rem_samples)

case_id <- unique(colData_tech$Case_ID) # How many Case_Ids remaining
print(length(case_id))



##### GEX_norm (all counts) ##### 
dds <- DESeqDataSetFromMatrix(countData = GEX_tech_rem, 
                              colData = colData_tech, 
                              design = ~ 1)
ddsBlind <- DESeq(dds)

# Normalised count matrix
GEX_norm <- counts(ddsBlind, normalized = TRUE)
write.csv(GEX_norm, file = paste0(data, "GEX/TechMerged_Unfiltered_NormCounts_GraftHNSCC_270424.csv"), quote=F, row.names = TRUE)

