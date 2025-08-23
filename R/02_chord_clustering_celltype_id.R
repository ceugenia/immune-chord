# The clustering and cell type identification script for chord
# Part of the chord pipeline

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load processed data
processed_data <- readRDS("data/processed_data/pancreas_processed.rds")

# Perform linear dimensional reduction (PCA)
processed_data <- RunPCA(processed_data, verbose = FALSE)

# Examine PCA results
pca_elbow <- ElbowPlot(processed_data, ndims = 50)
ggsave("figures/02_pca_elbow_plot.png", plot = pca_elbow, width = 8, height = 6)

# Cluster the cells
processed_data <- FindNeighbors(processed_data, dims = 1:20)
processed_data <- FindClusters(processed_data, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP)
processed_data <- RunUMAP(processed_data, dims = 1:20)

# Visualize clusters
umap_clusters <- DimPlot(processed_data, reduction = "umap", label = TRUE)
ggsave("figures/02_umap_initial_clusters.png",
       plot = umap_clusters, width = 8, height = 6)

# Find markers for each cluster
cluster_markers <- FindAllMarkers(processed_data,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

# Save marker results
write.csv(cluster_markers, "data/processed_data/initial_cluster_markers.csv")

# Attempt basic annotation using known marker genes
# Example: Endocrine cell markers
marker_genes <- c("INS", "GCG", "SST", "PPY", "GHRL") # Beta, Alpha, Delta, PP, Epsilon

# Visualize marker expression
feature_plots <- FeaturePlot(processed_data, features = marker_genes, ncol = 3)
ggsave("figures/02_feature_plots_markers.png",
       plot = feature_plots, width = 12, height = 8)

# Save object with clusters
saveRDS(processed_data, "data/processed_data/pancreas_clustered.rds")

# Print session info for reproducibility
sessionInfo()
