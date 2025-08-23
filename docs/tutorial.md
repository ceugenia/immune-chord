## Recommended Datasets

1. **BaronPancreasData** (via `scRNAseq` package) – easiest for testing.
2. **10X PBMC** – filtered feature matrix from [10x Genomics](https://www.10xgenomics.com/datasets).
3. **Tabula Sapiens** – download single-organ H5AD files (e.g., `TS_Lung.h5ad`).

---
  
  ## Next Steps
  
  1. Set up Conda environment and install packages.
2. Test BigSurR on BaronPancreasData.
3. Develop and test each pipeline script.
4. Write and render `vignette.Rmd`.
5. Push to GitHub with MIT license.

---
  
  ## Package Installation
  
  First, install the required packages:
  
  ```
# Create Conda environment
conda create -n scrare-env -c conda-forge r-base=4.3.2 r-essentials
conda activate scrare-env

# Install R packages
conda install -c conda-forge r-seurat r-tidyverse r-devtools r-remotes r-biocmanager
conda install -c bioconda bioconductor-singlecellexperiment bioconductor-scran

# Install BigSur (note: package name is BigSur, not BigSurR)
remotes::install_github("landerlabcode/BigSur")

```

# Data Loading

We'll use the BaronPancreasData for this tutorial:

```
library(scRNAseq)
library(Seurat)

# Load pancreas data
pancreas_data <- BaronPancreasData(which = "human")
seu_obj <- CreateSeuratObject(counts = counts(pancreas_data),
                             meta.data = as.data.frame(colData(pancreas_data)))

```

# Quality Control and Normalization

[Include standard QC steps from your existing tutorial]

# Rare Cell Detection with BigSur

## Running BigSur

```
library(BigSur)

# Filter zero-count genes first (critical step for BigSur)
count_matrix <- GetAssayData(seu_obj, assay = "RNA", slot = "counts")
zero_count_genes <- which(Matrix::rowSums(count_matrix) == 0)
if(length(zero_count_genes) > 0) {
  count_matrix <- count_matrix[-zero_count_genes, ]
  filtered_data <- CreateSeuratObject(counts = count_matrix)
  filtered_data@meta.data <- seu_obj@meta.data
  seu_obj <- filtered_data
}

# Run BigSur
bigsur_result <- BigSur(
  seurat.obj = seu_obj,
  assay = "RNA",
  counts.slot = "counts",
  variable.features = TRUE,
  correlations = FALSE,
  fano.alpha = 0.05,
  min.fano = 1.5
)

```

## Clustering and Rare Cell Identification

```
# Get variable features identified by BigSur
variable_genes <- VariableFeatures(bigsur_result)

# Scale the data using the variable features
bigsur_result <- ScaleData(bigsur_result, features = variable_genes)

# Run PCA with the scaled data
bigsur_result <- RunPCA(bigsur_result, features = variable_genes)

# Find neighbors and clusters
bigsur_result <- FindNeighbors(bigsur_result, dims = 1:20)
bigsur_result <- FindClusters(bigsur_result, resolution = 1.2)

# In our analysis, this resulted in 27 communities
cluster_count <- length(levels(bigsur_result@meta.data$seurat_clusters))
print(paste("Number of communities:", cluster_count))

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

# Check if any rare cells were identified
print("Number of rare cells:")
print(table(bigsur_result$is_rare))

```

## Troubleshooting: When No Rare Cells Are Detected

In our analysis of the BaronPancreasData, we found that no clusters met our initial criteria for rare cells (size < 5). This is a common scenario that requires adaptive analysis strategies.

### Adaptive Analysis Approach

```
# If no rare cells detected, analyze the smallest clusters instead
if (!"is_rare" %in% colnames(bigsur_result@meta.data) ||
    all(bigsur_result$is_rare == "Common")) {

  message("No rare cells identified with current threshold. Analyzing smallest clusters instead.")

  # Identify the smallest clusters (bottom 5% of cluster sizes)
  cluster_sizes <- table(bigsur_result$seurat_clusters)
  size_threshold <- quantile(cluster_sizes, probs = 0.05)
  small_clusters <- names(cluster_sizes[cluster_sizes <= size_threshold])

  # Mark cells in smallest clusters as "Potential Rare"
  bigsur_result$cell_status <- ifelse(
    bigsur_result$seurat_clusters %in% small_clusters,
    "Potential Rare",
    "Common"
  )

  print(paste("Analyzing", length(small_clusters),
              "small clusters as potential rare populations"))
  print("Cluster sizes:")
  print(cluster_sizes[small_clusters])

} else {
  bigsur_result$cell_status <- bigsur_result$is_rare
}

```

### Parameter Tuning Recommendations

Based on our experience with the BaronPancreasData:

- **`fano.alpha = 0.05`**: Balanced false discovery rate
- **`min.fano = 1.5`**: Captures meaningful variation without being too restrictive
- **`resolution = 1.2`**: Resulted in 27 communities for our dataset
- **Cluster size threshold**: May need adjustment based on dataset characteristics

### Alternative Approaches

```
# 1. Percentage-based threshold
rare_threshold <- ncol(bigsur_result) * 0.01  # 1% of total cells
rare_clusters <- names(cluster_sizes[cluster_sizes < rare_threshold])

# 2. Focus on smallest N clusters
small_clusters <- names(sort(cluster_sizes)[1:3])  # 3 smallest clusters

# 3. Outlier detection in PCA space
pca_scores <- Embeddings(bigsur_result, "pca")[, 1:2]
outliers <- apply(pca_scores, 2, function(x) {
  which(x > mean(x) + 3 * sd(x) | x < mean(x) - 3 * sd(x))
})
rare_cells <- unique(unlist(outliers))

```

# Differential Expression Analysis

```
# Perform differential expression analysis
de_results <- FindMarkers(
  bigsur_result,
  ident.1 = ifelse("Potential Rare" %in% bigsur_result$cell_status,
                   "Potential Rare", "Rare"),
  group.by = "cell_status",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Save and examine results
write.csv(de_results, "data/processed_data/differential_expression_results.csv")
top_genes <- de_results %>% arrange(desc(avg_log2FC)) %>% head(20)
print("Top differentially expressed genes:")
print(top_genes)

```

# Visualization

```
# UMAP with cell status
status_plot <- DimPlot(
  bigsur_result,
  group.by = "cell_status",
  cols = c("Common" = "grey", "Potential Rare" = "red", "Rare" = "red"),
  pt.size = 1
) + ggtitle("Cell Status on UMAP (27 Communities)")

print(status_plot)

# Feature plots of top genes
if (nrow(top_genes) > 0) {
  feature_plots <- FeaturePlot(
    bigsur_result,
    features = rownames(top_genes)[1:4],
    reduction = "umap",
    ncol = 2,
    order = TRUE
  ) + plot_annotation(title = "Expression of Top Differential Genes")

  print(feature_plots)
}

```

# Biological Interpretation

In our analysis of the BaronPancreasData:

- **Dataset**: 8,569 cells across 27 communities
- **Rare cells detected**: 0 using strict threshold (size < 5)
- **Adaptive approach**: Analyzed smallest clusters as potential rare populations
- **Key finding**: The pipeline successfully identified differentially expressed genes in small clusters, suggesting specialized cell subpopulations

# Conclusion

The immune-chord pipeline with BigSur provides a robust framework for rare cell identification, with adaptive strategies for handling datasets where strict rare cell criteria aren't met. The key advantages are:
  
  1. **Flexibility**: Adaptive analysis when no rare cells are detected
2. **Reproducibility**: Standardized workflow ensures consistent results
3. **Biological relevance**: Marker gene analysis provides meaningful insights
4. **Documentation**: Comprehensive troubleshooting guides for common scenarios

# References

1. Baron, M. et al. (2016) A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell Systems.
2. BigSur package documentation: https://github.com/landerlabcode/BigSur
3. Seurat documentation: https://satijalab.org/seurat/