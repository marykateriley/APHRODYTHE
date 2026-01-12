##########################################################
# Gene Level PDX, human and TCGA                         #
# Filtered out HPV/Lymphoma                              #
#                                                        #
# Gains / Losses CNA Gene Level Pearsons Corr            #
# PN0039A - Basal HNSCC                                  #
#                                                        #
# Mary-Kate Riley 20-10-25                               #
##########################################################

setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

#BiocManager::install(update = TRUE, ask = FALSE)
# Load packages
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
#BiocManager::install("maftools")
library(maftools)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
#BiocManager::install("NMF")
library(NMF)
library(ComplexHeatmap)
library(pheatmap)
library(gplots)
library(ggrepel) # geom_text_repel()
library(ggrastr)   # for rasterise(geom_point(...), dpi=300)
library(circlize) #colourRamp2()
#library(RColorBrewer)
library(gtools) # orders chromosomes using mixedsort()

#BiocManager::install("karyoploteR")
library(karyoploteR)
library(GenVisR)
library(GenomicRanges)
library(data.table)
library(gwascat)
library(rtracklayer)
library(AnnotationHub)
library(corrr)
library(cowplot)

########################################################## READ DATA ##########################################################
ClinicalInfo <- read.csv("./0_data/Complete_Curated_Clinical_Information.csv", header = T)
ClinicalInfo <- ClinicalInfo %>% 
  filter(WES == "Yes") %>% 
  filter(Remove != "Yes")

tcga_meta <- read.csv("./0_data/hnsc_tcga_pan_can_atlas_2018/HPVRem/clinicalData_mutation_filtered.csv")


#  Load thresholded GISTIC matrix
pdx_thr <- read.delim("./0_data/GISTIC/PDX_HPVLymp_Rem/all_thresholded.by_genes.txt", check.names = FALSE)
human_thr <- read.delim("./0_data/GISTIC/Human_HPVLymp_Rem/all_thresholded.by_genes.txt", check.names = FALSE)
tcga_thr <- read.delim("./0_data/hnsc_tcga_pan_can_atlas_2018/data_cna.txt", check.names = FALSE)


output <- "./3_output/GISTIC/NEW_FilteredStats/GeneLevel/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(output, "CNA_freq_PDX_vs_Human.pdf")
out_csv <- file.path(output, "CNA_freq_PDX_vs_Human_source.csv")
out_tsv <- file.path(output, "CNA_freq_PDX_vs_Human_correlations.tsv")  # new

compute_freq <- function(df, gene_col, direction) {
  gene_col <- rlang::ensym(gene_col)
  
  # Annotation
  annot_cols <- c("Gene Symbol","Gene.Symbol","Hugo_Symbol","Hugo Symbol",
                  "Gene","Gene_Symbol",
                  "Locus ID","Locus.ID",
                  "Cytoband",
                  "Entrez_Gene_Id","Entrez Gene Id")
  
  # Keep only sample columns
  mat <- df %>%
    select(-any_of(annot_cols)) %>%
    select(where(is.numeric)) %>%
    as.matrix()
  
  
  # Per-row frequency of gain or loss
  row_freq <- if (direction == 1) {
    rowMeans(mat >=  1, na.rm = TRUE)   # gains (any-level)
  } else {
    rowMeans(mat <= -1, na.rm = TRUE)   # losses (any-level)
  }
  
  # collapse to one row per gene (duplicates)
  tibble(
    Gene = df %>% pull(!!gene_col),
    freq = row_freq) %>%
    group_by(Gene) %>%
    summarise(freq = mean(freq, na.rm = TRUE), .groups = "drop")
}


##########################################################
## PDX vs HUMAN – gene-level CNA frequencies
##########################################################

## Amp and Del frequencies per gene, per cohort

# Gains (>= +1)
pdx_amp_PH <- compute_freq(pdx_thr, "Gene Symbol", direction =  1) %>%
  transmute(Gene, PDX   = freq)
human_amp_PH <- compute_freq(human_thr, "Gene Symbol", direction =  1) %>%
  transmute(Gene, Human = freq)

# Intersection of genes only
freq_PH_amp <- inner_join(pdx_amp_PH, human_amp_PH, by = "Gene") %>%
  mutate(event = "Amp")

# Losses (<= -1)
pdx_del_PH   <- compute_freq(pdx_thr,   "Gene Symbol", direction = -1) %>%
  transmute(Gene, PDX   = freq)
human_del_PH <- compute_freq(human_thr, "Gene Symbol", direction = -1) %>%
  transmute(Gene, Human = freq)

freq_PH_del <- inner_join(pdx_del_PH, human_del_PH, by = "Gene") %>%
  mutate(event = "Del")

# Amp + Del
freq_PH <- bind_rows(freq_PH_amp, freq_PH_del)

## Pearson Correlations
cor_PH_amp <- cor.test(freq_PH_amp$PDX, freq_PH_amp$Human)
cor_PH_del <- cor.test(freq_PH_del$PDX, freq_PH_del$Human)

fmt_r <- function(x) sprintf("%.2f", unname(x))

title_PH <- sprintf("Gene-level CNA frequencies: PDX vs Human\nAmp: r=%s  Del: r=%s",
                    fmt_r(cor_PH_amp$estimate),
                    fmt_r(cor_PH_del$estimate))

## lot PDX vs Human (Amp + Del)
# Small jitter tuned to sample sizes
n_human <- ncol(human_thr)
n_pdx <- ncol(pdx_thr)
jit_x <- min(0.35 / max(1, n_human), 0.003)
jit_y <- min(0.35 / max(1, n_pdx),   0.003)

p_PDX_Human <- ggplot(freq_PH, aes(x = Human, y = PDX, colour = event)) +
  ggrastr::rasterise(geom_point(alpha = 0.55, size = 0.35, shape = 16, position = position_jitter(jit_x, jit_y)), dpi = 900) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  coord_equal(xlim = c(0, 0.8), ylim = c(0, 1), expand = 0) +
  scale_colour_manual(values = c(Amp = "red", Del = "blue")) +
  theme_classic(base_size = 12) +
  labs(title = title_PH,
       x = "Human (gene-level CNA frequency)",
       y = "PDX (gene-level CNA frequency)",
       colour = NULL)

print(p_PDX_Human)
ggsave(file.path(output, "CNA_freq_PDX_vs_Human.pdf"), p_PDX_Human, width = 8, height = 8, units = "in", device = cairo_pdf)


##########################################################
## PDX vs TCGA – gene-level CNA frequencies
## First subset TCGA to the samples present in tcga_meta$Sample.ID
##########################################################

##Subset TCGA to matched samples

# AnnotationS
gene_id_cols <- intersect(colnames(tcga_thr),
                          c("Hugo_Symbol","Hugo Symbol","Gene Symbol","Gene.Symbol","Gene","Gene_Symbol",
                            "Entrez_Gene_Id","Entrez Gene Id","Cytoband","Locus ID","Locus.ID"))
sample_cols_all <- setdiff(colnames(tcga_thr), gene_id_cols)

tcga_keep <- unique(tcga_meta$Sample.ID)
matched_samples <- intersect(sample_cols_all, tcga_keep)

if (length(matched_samples) == 0) {
  stop("No TCGA columns matched tcga_meta$Sample.ID. Check ID formatting.")
}

tcga_thr_sub <- tcga_thr[, c(gene_id_cols, matched_samples), drop = FALSE]

message(sprintf("TCGA samples: %d -> kept %d", length(sample_cols_all), length(matched_samples)))

## Amp and Del frequencies per gene, per cohort

# Gains (>= +1)
pdx_amp_PT <- compute_freq(pdx_thr, "Gene Symbol", direction =  1) %>%
  transmute(Gene, PDX  = freq)
tcga_amp_PT  <- compute_freq(tcga_thr_sub, "Hugo_Symbol",  direction =  1) %>%
  transmute(Gene, TCGA = freq)

freq_PT_amp <- inner_join(pdx_amp_PT, tcga_amp_PT, by = "Gene") %>%
  mutate(event = "Amp")

# Losses (<= -1)
pdx_del_PT   <- compute_freq(pdx_thr, "Gene Symbol",  direction = -1) %>%
  transmute(Gene, PDX  = freq)
tcga_del_PT  <- compute_freq(tcga_thr_sub, "Hugo_Symbol",  direction = -1) %>%
  transmute(Gene, TCGA = freq)

freq_PT_del <- inner_join(pdx_del_PT, tcga_del_PT, by = "Gene") %>%
  mutate(event = "Del")

# Amp + Del
freq_PT <- bind_rows(freq_PT_amp, freq_PT_del)

## Pearson Correlations
cor_PT_amp <- cor.test(freq_PT_amp$PDX, freq_PT_amp$TCGA)
cor_PT_del <- cor.test(freq_PT_del$PDX, freq_PT_del$TCGA)

title_PT <- sprintf("Gene-level CNA frequencies: PDX vs TCGA\nAmp: r=%s Del: r=%s",
                    fmt_r(cor_PT_amp$estimate),
                    fmt_r(cor_PT_del$estimate))

## Plot PDX vs TCGA (Amp + Del)
n_tcga <- ncol(tcga_thr)
n_pdx  <- ncol(pdx_thr)
jit_x2 <- min(0.35 / max(1, n_tcga), 0.003)
jit_y2 <- min(0.35 / max(1, n_pdx),  0.003)

p_PDX_TCGA <- ggplot(freq_PT, aes(x = TCGA, y = PDX, colour = event)) +
  ggrastr::rasterise(geom_point(alpha = 0.55, size = 0.35, shape = 16, position = position_jitter(width = jit_x2, height = jit_y2)), dpi = 900) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.35) +
  coord_equal(xlim = c(0, 0.8), ylim = c(0, 1), expand = 0) +
  scale_colour_manual(values = c(Amp = "red", Del = "blue")) +
  theme_classic(base_size = 12) +
  labs(title = title_PT,
       x = "TCGA (gene-level CNA frequency)",
       y = "PDX (gene-level CNA frequency)",
       colour = NULL)

print(p_PDX_TCGA)
ggsave(file.path(output, "CNA_freq_PDX_vs_TCGA.pdf"), p_PDX_TCGA, width = 8, height = 8, units = "in", device = cairo_pdf)





########################################################### TCGA PIE CHART ##########################################################
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

# Access the clinical data
clinical_data <- tcga_meta

# Prepare data
pie_data <- clinical_data %>%
  dplyr::count(cases.primary_site) %>%
  dplyr::mutate(percentage = n / sum(n) * 100,
                label = paste0(round(percentage, 1), "%"),
                legend_label = paste0(cases.primary_site, " (n=", n,")"))

# Vector for custom legend labels
names(pie_data$legend_label) <- pie_data$cases.primary_site
pie_data$legend_label_wrapped <- str_wrap(pie_data$legend_label, width = 35) # Wrap after 40 cahracters, for long labels
# Get total sample size
total_samples <- sum(pie_data$n)

# Pie chart
p <- ggplot(pie_data, aes(x = "", y = n, fill = cases.primary_site)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = primary_site_colors, labels = pie_data$legend_label_wrapped) +
  geom_text(aes(label = ifelse(percentage >= 5, label, "")),  # only label >=5% - skip smaller slices
            position = position_stack(vjust = 0.5), size = 7) +
  labs(title = paste0("Copy Number Distribution of Primary Site (n=", total_samples,")"),
       subtitle = "TCGA PanCancer Atlas: HPV+, Lip and Bone/Cartilage Samples Removed",
       fill = "Primary Site") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
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


# Create a function to map values to your custom groups
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
  labs(title = paste0("Copy Number Grouped Primary Site (n=", total_samples,")"),
       subtitle = "TCGA PanCancer Atlas: HPV+, Lip and Bone/Cartilage Samples Removed",
       fill = "Primary Site") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"), # larger color boxes
        legend.spacing.y = unit(0.4, "cm"), # extra vertical spacing between legend rows
        legend.box.spacing = unit(0.5, "cm")) # space between legend and plot
ggsave(paste0(output, "PieChart_PanCancerAtlas_Grouped.pdf"), plot = p, width = 10.5, height = 10)
