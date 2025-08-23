# The quality control (QC) and normalization script for chord
# Part of the chord pipeline

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)

# Source helper functions
source("R/functions.R")

# Load the raw data (example using BaronPancreasData)
if (!require("scRNAseq")) BiocManager::install("scRNAseq")
library(scRNAseq)

# Load pancreas data
pancreas_data <- BaronPancreasData(which = "human")

# Convert to Seurat object
seu_obj <- CreateSeuratObject(counts = counts(pancreas_data),
                              meta.data = as.data.frame(colData(pancreas_data)))

# Process the data using our custom function
processed_data <- process_sc_data(seu_obj,
                                  min_features = 100,
                                  max_mito = 15,
                                  n_features = 2000)

# Create QC visualizations
create_qc_plots(processed_data, output_dir = "figures")

# Save the processed object
saveRDS(processed_data, "data/processed_data/pancreas_processed.rds")

# Print session info for reproducibility
sessionInfo()
