##########################################################
# Transcription Factor Activity + PROGENy                #
# decoupleR                                              #
#                                                        #
# NonResponders v Responders                             #
# Avg 3-6 week response                                  #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 05-11-25                               #
##########################################################

setwd("/Users/mary-kateriley/Library/CloudStorage/OneDrive-Queen'sUniversityBelfast/Data Analysis/HNC/Basal PDX/Reanalysis Basal '24")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
set.seed(1234)

# Load packages
library(tidyverse)
library(dplyr)
library(tidyr)
library(decoupleR)
library(OmnipathR) # decoupleR depends on to retrieve DoRothEA network https://omnipathdb.org/
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)



########################################################## READ DATA ##########################################################

# Dont set genes as row names - keep as first column
DESeq <- read.csv("./2_pipelines/DESeq2/Filtered/Response/NRvR/DESeq2_Definitive_Response_3_6wkNonResponder_v_Responder.csv", header = T)
DESeq[1:5,1:5] 

GEX <- read.csv("./2_pipelines/DESeq2/Filtered/Response/NRvR/Norm_counts_Definitive_Response_3_6wkNonResponder_v_Responder.csv", header = T, row.names = 1)


colData <- read.csv("./0_data/DESeq2/PRX/Metadata/All_RemovedSamples_Unfiltered_GraftHNSCC_colData_240724.csv", header = T, row.names = 1)

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/TF_Activity/NRvR/"
output <- "./3_output/TF_Activity/NRvR/"
title <- "Basal HNSCC PDX: Average 3-6 Week NonResponder v Responder"
subtitle <- "Filtered Technical; Lymphoma; HPV"
save_title <- "HNC_NRvR"

dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

########################################################## TRANSCRIPTION FACTORS ##########################################################

# CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources. 
# Increased coverage of transcription factors and a superior performance in identifying perturbed TFs compared to our previous DoRothEA network
net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE) # splits into sub units, recommended to keep them together


#format these results before using decoupleR
DESeq_new <- DESeq %>% 
  dplyr::select(X, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "X") %>%
  as.matrix()

colnames(DESeq_new)[1] <- "t" #re-name the column "t" instead of stat



# Run ulm - Univariate Linear Model (ULM)
# Infers TF activity using the DESeq2 stat 
# Between Responders and NonResponders
# Creates TF enrichment "score"
contrast_acts <- decoupleR::run_ulm(mat = DESeq_new[, 't', drop = FALSE], 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor='mor', 
                                    minsize = 5)
contrast_acts


# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  dplyr::mutate(rnk = NA) # creates a "rnk" column and fills it with NA values

# Creates a logical mask: TRUE = positive score
msk <- f_contrast_acts$score > 0 

# Positive msk scores - ranks are assigned by taking the negative of the score. Larger scores will get smaller (better/higher) ranks
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score']) 

# Negative msk scores - ranks assigned taking the absolute value of the score. All negative scores are treated as positive (for ranking)
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
write.csv(f_contrast_acts, file=paste0(pipelines,save_title,"_TF_Activities.csv"), row.names = FALSE)

# Top 25 TF factors
n_tfs <- 25

tfs <- f_contrast_acts %>%
  dplyr::arrange(rnk) %>% # sorts it by smallest rank first
  head(n_tfs) %>% # selects top 25 rows after sorting
  dplyr::pull(source) # takes source column into a vector

# Only keep TF that are in top 25
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)
write.csv(f_contrast_acts, file=paste0(pipelines,save_title,"_Top25_TF_Activities.csv"), row.names = FALSE)


colours <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])


p <- ggplot(data = f_contrast_acts, 
            mapping = aes(x = stats::reorder(source, score), y = score)) + 
  geom_bar(mapping = aes(fill = score), colour = "black",  stat = "identity") +
  scale_fill_gradient2(low = colours[1], mid = "whitesmoke", high = colours[2], midpoint = 0) + 
  theme_bw() +
  labs(title = paste0(title),
       subtitle = paste0(subtitle, " Top 25 Transcription Factors")) +
  theme(axis.line = element_line(),
        plot.title = element_text(size = 20, face="bold"),
        plot.subtitle = element_text(size = 18),
        axis.title.y = element_text(size = 18,  colour="black", face="bold"),
        axis.text.y = element_text(size=16,colour="black"),
        axis.title.x = element_text(size = 18,  colour="black", face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=12,colour="black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank()) +
  xlab("TFs") +
  theme(panel.grid.major.y = element_line(colour = "gray80"))
ggsave(paste0(output, "VolcanoPlot_Top25_TF_", save_title, ".pdf"), plot = p, width = 11, height = 11)

