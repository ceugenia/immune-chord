# The differential expression analysis and visualization script for chord
# Part of the chord pipeline

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(viridis)

# Load data with rare cells identified
final_data <- readRDS("data/processed_data/pancreas_27_communities.rds")

# Check if rare cells were identified
if (!"is_rare" %in% colnames(final_data@meta.data) || 
    all(final_data$is_rare == "Common")) {
  
  message("No rare cells identified with current threshold. Analyzing smallest clusters instead.")
  
  # Identify the smallest clusters (bottom 5% of cluster sizes)
  cluster_sizes <- table(final_data$seurat_clusters)
  size_threshold <- quantile(cluster_sizes, probs = 0.05)
  small_clusters <- names(cluster_sizes[cluster_sizes <= size_threshold])
  
  # Mark cells in smallest clusters as "Potential Rare"
  final_data$cell_status <- ifelse(
    final_data$seurat_clusters %in% small_clusters,
    "Potential Rare",
    "Common"
  )
  
  print(paste("Analyzing", length(small_clusters), 
              "small clusters as potential rare populations"))
  print("Cluster sizes:")
  print(cluster_sizes[small_clusters])
  
} else {
  final_data$cell_status <- final_data$is_rare
}

# Perform differential expression analysis
de_results <- FindMarkers(
  final_data,
  ident.1 = ifelse("Potential Rare" %in% final_data$cell_status, 
                   "Potential Rare", "Rare"),
  group.by = "cell_status",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Save the results
write.csv(de_results, "data/processed_data/differential_expression_results.csv")

# Get top differentially expressed genes
top_genes <- de_results %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

print("Top differentially expressed genes:")
print(top_genes)

# Create visualizations

# 1. UMAP colored by cell status
status_plot <- DimPlot(
  final_data,
  group.by = "cell_status",
  cols = c("Common" = "grey", "Potential Rare" = "red", "Rare" = "red"),
  pt.size = 1
) + ggtitle("Cell Status on UMAP")

ggsave("figures/04_umap_cell_status.png", 
       plot = status_plot, 
       width = 10, 
       height = 8)

# 2. Feature plots of top genes
if (nrow(top_genes) > 0) {
  feature_plots <- FeaturePlot(
    final_data,
    features = rownames(top_genes)[1:4],
    reduction = "umap",
    ncol = 2,
    order = TRUE
  ) + plot_annotation(title = "Expression of Top Differential Genes")
  
  ggsave("figures/04_feature_plots_de_genes.png", 
         plot = feature_plots, 
         width = 10, 
         height = 8)
}

# 3. Violin plots of top genes
if (nrow(top_genes) > 0) {
  violin_plots <- VlnPlot(
    final_data,
    features = rownames(top_genes)[1:4],
    group.by = "cell_status",
    pt.size = 0,
    ncol = 2
  ) + plot_annotation(title = "Expression Differences by Cell Status")
  
  ggsave("figures/04_violin_de_genes.png", 
         plot = violin_plots, 
         width = 10, 
         height = 8)
}

# Biological interpretation
if (nrow(top_genes) > 0) {
  print("Biological Interpretation:")
  print(paste("The analysis identified", nrow(de_results), 
              "differentially expressed genes."))
  print(paste("Top genes include:", 
              paste(rownames(top_genes)[1:3], collapse = ", ")))
  print("These genes may represent markers of rare or specialized cell populations.")
} else {
  print("No significant differentially expressed genes found.")
  print("This could indicate that:")
  print("1. The dataset doesn't contain distinct rare cell populations")
  print("2. The clustering resolution needs adjustment")
  print("3. Different parameters should be used for differential expression")
}

# Save final results
saveRDS(final_data, "data/processed_data/pancreas_final_analysis.rds")

# Print session info for reproducibility
sessionInfo()
