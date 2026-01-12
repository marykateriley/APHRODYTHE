##########################################################
# Responders v Non-Responders clinical data              #
# 3, 6 and Avg 3-6 week response                         #
# LAB + SEQUENCING                                       #
# PR, SD, PD                                             #
#                                                        #
# PN0039 - Basal HNSCC                                   #
#                                                        #
# Mary-Kate Riley 02-12-25                               #
##########################################################

setwd("")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load packages
library(tidyverse)
library(ggplot2)
library(rstatix) # wilcox_test
library(ggpubr) # stat_pvalue_manual
library(patchwork)


########################################################## READ DATA ##########################################################
metadata <- read.csv("./0_data/SequencingAndLab_ClinicalInformation.csv", header = T)

metadata <- metadata %>% 
  filter(Include != "N")

########################################################## SET VARIABLES ##########################################################

output <- "./3_output/Response_Clinical_Plots/SeqAndLab_Cases_Response_PRSDPD/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

########################################################## CURATION ##########################################################

# Set up response categories
metadata <- metadata %>%
  dplyr::rename(Case_ID = Case) %>% 
  dplyr::mutate(Age_Group = cut(Age,breaks = c(30, 40, 50, 60, 70, 80, 90, 100),right  = FALSE,
                                labels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")))
# Factor clinical information
colnames(metadata)

metadata$Sex <- factor(metadata$Sex)
levels(metadata$Sex)

metadata$Age_Group <- factor(metadata$Age_Group)
levels(metadata$Age_Group)

metadata$T_Stage <- factor(metadata$T_Stage)
levels(metadata$T_Stage)
levels(metadata$T_Stage)[levels(metadata$T_Stage) == ""] <- "Unknown"
levels(metadata$T_Stage)

metadata$N_Stage <- factor(metadata$N_Stage)
levels(metadata$N_Stage)
levels(metadata$N_Stage)[levels(metadata$N_Stage) == ""] <- "Unknown"
levels(metadata$N_Stage)

metadata$HISTOLOGICAL.GRADE <- factor(metadata$HISTOLOGICAL.GRADE)
levels(metadata$HISTOLOGICAL.GRADE)
levels(metadata$HISTOLOGICAL.GRADE)[levels(metadata$HISTOLOGICAL.GRADE) == ""] <- "Unknown"
levels(metadata$HISTOLOGICAL.GRADE)

# Remove rows with NA in the response columns (Def_Response_3wk, Def_Response_6wk, Def_Response_3_6wk)
metadata_cleaned_avg <- metadata %>%
  dplyr::filter(!is.na(Response_3_6wk), Response_3_6wk != "")
metadata_avg <- metadata %>%
  dplyr::mutate(Response_3_6wk = dplyr::if_else(
    is.na(Response_3_6wk) | Response_3_6wk == "","No_data", Response_3_6wk))

metadata_cleaned_avg$Response_3_6wk <- factor(metadata_cleaned_avg$Response_3_6wk, levels = c("PD", "SD", "PR"))
metadata_avg$Response_3_6wk <- factor(metadata_avg$Response_3_6wk, levels = c("No_data", "PR", "SD", "PD"))

########################################################## DEFINE FUNCTION  ##########################################################

# Define color scheme for the plots
color_scheme <-  c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF") 

# Function to create count bar plots with customized labels
plot_categorical_timepoint <- function(df, variable, response_column, plot_label) {
  # Convert variables to symbols
  variable_sym <- sym(variable)
  response_column_sym <- sym(response_column)
  
  # Filter out NA from response column
  df <- df %>% filter(!is.na(!!response_column_sym))
  
  # Fisher's Exact Test
  tbl <- table(df[[as.character(variable)]], df[[as.character(response_column)]])
  fisher_result <- fisher.test(tbl)
  p_value <- fisher_result$p.value
  
  # Bar plot with counts
  p <- ggplot(df, aes(x = !!variable_sym, fill = !!response_column_sym)) +
    geom_bar(position = "stack", alpha = 0.85) +
    scale_fill_manual(values = color_scheme) +
    theme_minimal(base_family = "Arial") +
    labs(title = paste(plot_label, "\n(p =", format(p_value, digits = 3), ")"),
         x = NULL, y = "Number of Patients") +
    theme(legend.position = "right",
          legend.direction = "vertical",
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())
  
  return(p)
}


########################################################## CREATE BAR PLOTS SUBSET ##########################################################

# Define the variables to loop
variables <- c("Age_Group", "Sex", "Site.of.Primary", "T_Stage", "N_Stage")

# Rename
rename_map <- c("Age_Group" = "Age Group",
                "Sex" = "Sex",
                "Site.of.Primary" = "Primary Site",
                "T_Stage" = "T-Stage",
                "N_Stage" = "N-Stage")

# Create empty lists to store plots for each timepoint
plots_avg <- list()

# Function to generate plots with renamed variables
generate_plots <- function(df, response_col, variables, rename_map, timepoint) {
  plot_list <- list()
  for (var in variables) {
    renamed_var <- rename_map[[var]]
    plot_list[[var]] <- plot_categorical_timepoint(df, var, response_col, renamed_var)
  }
  return(plot_list)
}

# Generate plots for 3-week, 6-week, and average responses
plots_avg <- generate_plots(metadata_cleaned_avg, "Response_3_6wk", variables, rename_map, "average")

# Combine average plots
average_grid <- (plots_avg$Age_Group | plots_avg$Sex | plots_avg$Site.of.Primary) /
  (plots_avg$N_Stage | plots_avg$T_Stage ) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = paste0("Clinical Characteristics by Average 3-6 Weeks Response (n = ",nrow(metadata_cleaned_avg),")"),
                  subtitle = "Fisher's Exact Test",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                                plot.subtitle = element_text(size = 14, hjust = 0.5)))

# Print
average_grid

# Save
ggsave(paste0(output, "/seubset_Average_Grid.pdf"), plot = average_grid, width = 12, height = 10, device = grDevices::cairo_pdf)













########################################################## DEFINE FUNCTION ALL DATA ##########################################################

# Define color scheme for the plots
color_scheme <-  c("PR" = "#5AAE61FF", "SD" = "#F7F7F7FF", "PD" = "#762A83FF", "No_data" = "black")

# Function to create count bar plots with customized labels
plot_categorical_timepoint <- function(df, variable, response_column, plot_label) {
  variable_sym <- sym(variable)
  response_column_sym <- sym(response_column)
  df <- df %>% filter(!is.na(!!response_column_sym))
  
  # Fisher's Exact Test
  tbl <- table(df[[as.character(variable)]], df[[as.character(response_column)]])
  fisher_result <- fisher.test(tbl, workspace = 2e8)
  p_value <- fisher_result$p.value
  
  # Bar plot with counts
  p <- ggplot(df, aes(x = !!variable_sym, fill = !!response_column_sym)) +
    geom_bar(position = "stack", alpha = 0.85) +
    scale_fill_manual(values = color_scheme) +
    theme_minimal(base_family = "Arial") +
    labs(title = paste(plot_label, "\n(p =", format(p_value, digits = 3), ")"),
         x = NULL, y = "Number of Patients") +
    theme(legend.position = "right",
          legend.direction = "vertical",
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())
  return(p)
}

########################################################## CREATE BAR PLOTS ##########################################################

# Define the variables to loop 
variables <- c("Age_Group", "Sex", "Site.of.Primary", "T_Stage", "N_Stage")


# Rename map for the variables
rename_map <- c("Age_Group" = "Age Group",
                "Sex" = "Sex",
                "Site.of.Primary" = "Primary Site",
                "T_Stage" = "T-Stage",
                "N_Stage" = "N-Stage")

# Create empty lists to store plots for each timepoint
plots_avg <- list()

# Function to generate plots with renamed variables
generate_plots <- function(df, response_col, variables, rename_map, timepoint) {
  plot_list <- list()
  for (var in variables) {
    renamed_var <- rename_map[[var]]
    plot_list[[var]] <- plot_categorical_timepoint(df, var, response_col, renamed_var)
  }
  return(plot_list)
}

# Generate plots for 3-week, 6-week, and average responses
plots_avg <- generate_plots(metadata_avg, "Response_3_6wk", variables, rename_map, "average")

# Combine average plots
average_grid <-  (plots_avg$Age_Group | plots_avg$Sex | plots_avg$Site.of.Primary) /
  (plots_avg$N_Stage | plots_avg$T_Stage ) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = paste0("Clinical Characteristics by Average 3-6 Weeks Response (n = ",nrow(metadata_avg),")"),
                  subtitle = "Fisher's Exact Test",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                                plot.subtitle = element_text(size = 14, hjust = 0.5)))

# Print
average_grid

# Save
ggsave(paste0(output, "/Average_Grid_all.pdf"), plot = average_grid, width = 12, height = 10, device = grDevices::cairo_pdf)

