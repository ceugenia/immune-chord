# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load clustered data
clustered_data <- readRDS("data/processed_data/pancreas_clustered.rds")

# Filter out genes with zero counts across all cells
count_matrix <- GetAssayData(clustered_data, assay = "RNA", slot = "counts")
zero_count_genes <- which(Matrix::rowSums(count_matrix) == 0)
if(length(zero_count_genes) > 0) {
  count_matrix <- count_matrix[-zero_count_genes, ]
  filtered_data <- CreateSeuratObject(counts = count_matrix)
  filtered_data@meta.data <- clustered_data@meta.data
  clustered_data <- filtered_data
}

# Run BigSur to identify variable features
bigsur_result <- BigSur(
  seurat.obj = clustered_data,
  assay = "RNA",
  counts.slot = "counts",
  variable.features = TRUE,
  correlations = FALSE,
  fano.alpha = 0.05,
  min.fano = 1.5
)

# Get variable features identified by BigSur
variable_genes <- VariableFeatures(bigsur_result)

# Scale the data using the variable features
bigsur_result <- ScaleData(bigsur_result, features = variable_genes)

# Run PCA with the scaled data
bigsur_result <- RunPCA(bigsur_result, features = variable_genes)

# Find neighbors and clusters with higher resolution to get 15 communities
bigsur_result <- FindNeighbors(bigsur_result, dims = 1:20)
bigsur_result <- FindClusters(bigsur_result, resolution = 1.2)  # Increased resolution

# Check the number of clusters
cluster_count <- length(levels(bigsur_result@meta.data$seurat_clusters))
print(paste("Number of communities:", cluster_count))

# If we don't have exactly 15 clusters, adjust resolution
if (cluster_count != 15) {
  # Try to find the right resolution for 15 clusters
  resolutions <- seq(0.5, 2.0, by = 0.1)
  for (res in resolutions) {
    bigsur_result <- FindClusters(bigsur_result, resolution = res)
    cluster_count <- length(levels(bigsur_result@meta.data$seurat_clusters))
    if (cluster_count == 15) break
  }
}

# Run UMAP
bigsur_result <- RunUMAP(bigsur_result, dims = 1:20)

# Identify small clusters as potential rare cell populations
cluster_sizes <- table(bigsur_result$seurat_clusters)
rare_clusters <- names(cluster_sizes[cluster_sizes < 5])

# Mark rare cells
bigsur_result$is_rare <- ifelse(
  bigsur_result$seurat_clusters %in% rare_clusters,
  "Rare",
  "Common"
)

# Create UMAP visualization with all 15 communities
umap_plot <- DimPlot(bigsur_result, 
                     reduction = "umap", 
                     label = TRUE, 
                     label.size = 4,
                     pt.size = 1) +
  ggtitle("UMAP with 27 Communities") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(umap_plot)

# Save the plot
ggsave("figures/umap_27_communities.png", 
       plot = umap_plot, 
       width = 10, 
       height = 8)

# Create a separate plot highlighting rare cells
rare_plot <- DimPlot(bigsur_result, 
                     group.by = "is_rare", 
                     cols = c("grey", "red"),
                     pt.size = 1) +
  ggtitle("Rare Cell Populations") +
  theme(plot.title = element_text(hjust = 0.5))

print(rare_plot)
ggsave("figures/umap_rare_cells.png", 
       plot = rare_plot, 
       width = 10, 
       height = 8)

# Save the results with 15 communities
saveRDS(bigsur_result, "data/processed_data/pancreas_27_communities.rds")

# Print cluster sizes
print("Cluster sizes:")
print(table(bigsur_result$seurat_clusters))

# Print rare cluster information
print("Rare clusters (size < 5):")
print(rare_clusters)
