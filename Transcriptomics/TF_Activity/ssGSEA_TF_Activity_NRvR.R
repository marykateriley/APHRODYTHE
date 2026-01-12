##########################################################
# Single Sample Transcription Factor Activity            #
# decoupleR                                              #
#                                                        #
# NonResponders v Responders                             #
# Avg 3-6 week response                                  #
#                                                        #
# Technical, lymphoma + HPV samples removed              #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 14-12-25                               #
##########################################################

setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
set.seed(1234)

# Load packages
library(tidyverse)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(decoupleR)
library(OmnipathR) # decoupleR depends on to retrieve DoRothEA network https://omnipathdb.org/
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(rstatix)



########################################################## READ DATA ##########################################################

GEX <- read.csv("./0_data/edgeR/PRX/TPM/PDX_TPM_GEX.csv", row.names = 1)
rownames(GEX) <- toupper(rownames(GEX))
GEX <- log2(GEX + 1)

metadata <- read.csv("./0_data/DESeq2/PRX/Metadata/All_RemovedSamples_Unfiltered_GraftHNSCC_colData_240724.csv", header = T, row.names = 1)

########################################################## SET VARIABLES ##########################################################

pipelines <- "./2_pipelines/TF_Activity/VIPER/"
output <- "./3_output/TF_Activity/VIPER/"
title <- "VIPER Single Sample TF Activity Scores"
subtitle <- "Filtered Technical; Lymphoma; HPV"
save_title <- "HNC_TF_ssGSEA"
tfs <- "GLI2"

dir.create(pipelines, recursive = TRUE, showWarnings = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

source("./1_code/RFunctions/new_heatmap.R")

########################################################## DATA CURATION ##########################################################

# Remove samples with no Definitive_Response_3_6wk data
na_samples <- metadata$Sample_ID[is.na(metadata$Definitive_Response_3_6wk)]
GEX <- GEX[, !colnames(GEX) %in% na_samples]
metadata <- metadata %>% 
  filter(!Sample_ID %in% na_samples)

########################################################## TRANSCRIPTION FACTORS ##########################################################

# CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources. 
# Increased coverage of transcription factors and a superior performance in identifying perturbed TFs compared to our previous DoRothEA network
net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE) # splits into sub units, recommended to keep them together

# VIPER
viper <- decoupleR::run_viper(
  mat = GEX,
  network = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5)


# format into df
viper_wide <- viper %>%
  dplyr::select(source, condition, score) %>%
  tidyr::pivot_wider(names_from = source,
                     values_from = score) %>%
  tibble::column_to_rownames("condition")
viper_1 <- as.data.frame(viper_wide)   # samples x gene sets

#save single sample activity predictions
write.csv(viper_1, file = paste0(pipelines,"VIPER_Scores_Graft.csv"))

# Combine ssGSEA scores with metadata
viper_1$SampleID <- rownames(viper_1)

# Are all samples in colnames of GEX present rownames of colData and are they in the same order
all_samples_present <- all(rownames(metadata) %in% rownames(viper_1))
are_same_order <- identical(rownames(metadata), rownames(viper_1))
metadata <- metadata[match(rownames(viper_1), rownames(metadata)), ] # order
all_samples_present <- all(rownames(metadata) %in% rownames(viper_1)) # check again
are_same_order <- identical(rownames(metadata), rownames(viper_1))


combined_data <- merge(viper_1, metadata, by.x = "SampleID", by.y = "Sample_ID")
colnames(combined_data)
viper_1$SampleID <- NULL 

viper_wide <- viper %>%
  dplyr::select(condition, source, score) %>%
  tidyr::pivot_wider(names_from = source, values_from = score) %>%
  tibble::column_to_rownames("condition")

viper_2 <- as.data.frame(viper_wide)


########################################################## STATS ##########################################################

# Group column
group_col <- "Definitive_Response_3_6wk"   # <- change only this if needed

# Identify ssGSEA columns
gs_cols <- setdiff(colnames(combined_data), c(colnames(metadata), "SampleID"))

# Long format
long_viper <- combined_data %>%
  dplyr::select(SampleID, all_of(group_col), all_of(gs_cols)) %>%
  tidyr::pivot_longer(cols = all_of(gs_cols),
                      names_to = "GeneSet", values_to = "Score") %>%
  tidyr::drop_na(Score)

# Stats
viper_stats <- long_viper %>%
  dplyr::group_by(GeneSet) %>%
  rstatix::wilcox_test(as.formula(paste("Score ~", group_col))) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p.adj)
tf_positions <- match(tfs, viper_stats$GeneSet)

write.csv(viper_stats,file = file.path(pipelines, paste0("VIPER_stats_by_", group_col, ".csv")), row.names = FALSE)


# FDR
sig_sets <- viper_stats %>%
  dplyr::filter(p.adj < 0.05) %>%
  slice_head(n = 25) %>%
  dplyr::arrange(p.adj) %>%
  dplyr::pull(GeneSet) %>%
  unique()

# long_ssgsea
long_sig <- long_viper %>% dplyr::filter(GeneSet %in% sig_sets)


sig_sets_and_tfs <- viper_stats %>%
  dplyr::filter(GeneSet %in% unique(c(sig_sets, tfs))) %>%
  dplyr::filter(p.adj < 0.05) %>%
  dplyr::arrange(p.adj) %>%
  dplyr::pull(GeneSet) %>%
  unique()

########################################################## HEATMAP TOP 25 ANNOTATIONS ##########################################################

# Expression matrix: samples × gene sets
hm <- as.matrix(combined_data[, gs_cols, drop = FALSE])
rownames(hm) <- combined_data$SampleID

# Sample annotation
annotation_df <- data.frame(
  Definitive_Response_3_6wk = combined_data$Definitive_Response_3_6wk,
  Response_3_6wk = combined_data$Response_3_6wk,
  percentage_vol_change_3_6wk = combined_data$percentage_vol_change_3_6wk,
  Sex = combined_data$Sex,
  Site.of.Primary = combined_data$Site.of.Primary,
  row.names = combined_data$SampleID,
  check.names = FALSE)


# continuous mapper for % change
pal <- rev(RColorBrewer::brewer.pal(11, "PRGn"))
pct_col_fun <- circlize::colorRamp2(c(-90, 0, 335), c(pal[2], pal[6], pal[10]))

# numeric
annotation_df$percentage_vol_change_3_6wk <- as.numeric(as.character(annotation_df$percentage_vol_change_3_6wk))
# Annotation colors
ann_colors <- list(
  Model = c("PDX" = "#f5b041", "Human" = "#1abc9c"),
  Definitive_Response_3_6wk = c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF"),
  Response_3_6wk = c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF"),
  percentage_vol_change_3_6wk = pct_col_fun,  # function => continuous legend
  Sex = c("M" = "#A6CEE3", "F" = "#FB9A99"),
  Site.of.Primary = c("Larynx" = "#1B9E77FF", "Oropharynx" = "#66A61EFF", "Hypopharynx" = "#7570B3FF", "Oral Cavity" = "#E6AB02FF", "Larynx + Oral Cavity" = "#D95F02FF", "Oral Cavity + Oropharynx" = "#E7298AFF", "Unknown" = "#666666FF"))



# significant gene sets 
keep_cols <- intersect(colnames(hm), sig_sets)
hm <- hm[, keep_cols, drop = FALSE]

# order samples by group 
ord <- order(combined_data[[group_col]])
hm_ord <- hm[ord, , drop = FALSE]
annotation_df <- annotation_df[rownames(hm_ord), , drop = FALSE]

mat_plot <- t(hm_ord)  

breaks <- c(-3, 0, 3)
colors <- colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))

pdf(file.path(output, "VIPER_Top20_SingleSample_TFAvtiviy.pdf"), width = 9, height = 9)
new_heatmap(
  mat = mat_plot,
  scale = "row",
  color = colors,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  annotation_names_col = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  treeheight_row = 25,
  treeheight_col = 25,
  fontsize = 8,
  column_split = factor(annotation_df$Definitive_Response_3_6wk, levels = c("NonResponder", "Responder")),
  main = paste0(title, " Top 25 (p.adjust < 0.05)"),
  name = "Z-score",
  annotation_show_legend = c(percentage_vol_change_3_6wk = TRUE),
  annotation_legend_param = list(
    percentage_vol_change_3_6wk = list(
      title = "% volume change 3–6 wk",
      at = c(-90, 0, 100, 200, 335),
      labels = c("-90", "0", "100", "200", "335"),
      direction = "horizontal")))
dev.off()



########################################################## HEATMAP TOP 25 + TFS ANNOTATIONS ##########################################################
hm <- as.matrix(combined_data[, gs_cols, drop = FALSE])
rownames(hm) <- combined_data$SampleID

# significant gene sets
keep_cols <- intersect(colnames(hm), sig_sets_and_tfs)
hm <- hm[, keep_cols, drop = FALSE]

# order samples by group
ord <- order(combined_data[[group_col]])
hm_ord  <- hm[ord, , drop = FALSE]
annotation_df <- annotation_df[rownames(hm_ord), , drop = FALSE]


mat_plot <- t(hm_ord)

breaks <- c(-3, 0, 3)
colors <- colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))

pdf(file.path(output, "VIPER_Top20_AndTFs_SingleSample_TFAvtiviy.pdf"), width = 10, height = 9)
new_heatmap(
  mat = mat_plot,
  scale = "row",
  color = colors,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  annotation_names_col = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = NA,
  treeheight_row = 25,
  treeheight_col = 25,
  fontsize = 8,
  column_split = factor(annotation_df$Definitive_Response_3_6wk, levels = c("NonResponder", "Responder")),
  main = paste0(title, " Top 25 and GLI2 (p.adjust < 0.05)"),
  name = "Z-score",
  annotation_show_legend = c(percentage_vol_change_3_6wk = TRUE),
  annotation_legend_param = list(
    percentage_vol_change_3_6wk = list(
      title = "% volume change 3–6 wk",
      at = c(-90, 0, 100, 200, 335),
      labels = c("-90", "0", "100", "200", "335"),
      direction = "horizontal")))
dev.off()



plot_df <- viper_stats %>%
  arrange(p.adj) %>%
  mutate(rank = row_number(),
         sig  = p.adj < 0.05)

# TFs to label
label_tfs <- unique(c(plot_df$GeneSet[1:10], tfs))

plot_df <- plot_df %>%
  mutate(label = ifelse(GeneSet %in% label_tfs, GeneSet, NA_character_))

ggplot(plot_df, aes(x = rank, y = -log10(p.adj))) +
  geom_point(aes(colour = sig), alpha = 0.6) +
  geom_text(aes(label = label), vjust = -0.7, size = 3, check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "TF differential activity significance (VIPER)",
       subtitle = "Red points indicate FDR < 0.05; labels show top 10 TFs + selected TFs",
       x = "Rank by FDR", y = "-log10(FDR)") +
  theme_bw() +
  theme(legend.position = "none")



########################################################## TF VIOLIN PLOTS ##########################################################

resp_cols <- c("Responder" = "#1B7837FF", "NonResponder" = "#762A83FF")

df1 <- combined_data %>%
  select(SampleID, all_of(group_col), all_of(tfs)) %>%
  rename(Activity = all_of(tfs)) %>%
  filter(!is.na(Activity)) %>%
  mutate(group = factor(.data[[group_col]], levels = c("Responder", "NonResponder")))

pval <- viper_stats %>%
  dplyr::filter(GeneSet == tfs) %>%
  dplyr::transmute(
    group1 = "Responder",
    group2 = "NonResponder",
    p.adj = p.adj,
    p.adj.signif = p.adj.signif,
    label = paste0("p = ", formatC(p.adj, digits = 3, format = "g"),
                   " (", p.adj.signif, ")") )

# add y position 
y_pos <- max(df1$Activity, na.rm = TRUE) * 1.05
pval$y.position <- y_pos

plot <- ggplot(df1, aes(x = group, y = Activity, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.8, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 1.4, alpha = 0.5, colour = "black") +
  scale_fill_manual(values = resp_cols) +
  scale_colour_manual(values = resp_cols) +
  stat_pvalue_manual(pval,
                     label = "label",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01, 
                     bracket.size = 0.5, 
                     size = 4) +
  labs(title = paste0(title, ": ", tfs),
       x = NULL, y = "VIPER activity score") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x  = element_text(size = 14, face = "bold"),
        axis.text.y  = element_text(size = 14),
        legend.position = "none")
ggsave(file = paste0(output, "VIPER_Scores_GLI2.pdf"), plot, width = 7, height = 9)
