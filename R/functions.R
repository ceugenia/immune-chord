#' Functions related to the pre-process of single-cell sequencing of data

#' This performs the quality control (QC), then normalization, and finally scaling of the single cell RNA sequencing data

#' The following parameters are used:

#' @param seurat_object which is a raw Seurat object
#' @param min_features which is the minimum number of features per cell
#' @param max_mito_genes_genes which is the maximum percentage of mitochondrial genes
#' @param n_features which is the number of variable features to identify
#' 
#' @report returns a normalized and scaled Seurat object
#' @export
#'
#'

process_sc_data <- function(seurat_object, min_features = 200, max_mito_genes = 10, n_features = 2000) {
  # Calculate mitochondrial percentage
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  # Filter cells based on QC metrics
  seurat_object <- subset(seurat_object,
                          subset = nFeature_RNA > min_features &
                            percent.mt < max_mito_genes)
  
  # Normalize data
  seurat_object <- NormalizeData(seurat_object)
  
  # Find variable features
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures = n_features)
  
  # Scale data
  seurat_object <- ScaleData(seurat_object)
  
  return(seurat_object)
}

#' Create QC Plots
#'
#' Generates and saves quality control visualizations
#'
#' @param seurat_object A Seurat object
#' @param output_dir Directory to save plots
#'
#' @export
#'
create_qc_plots <- function(seurat_object, output_dir = "figures") {
  # Create Violin plots of QC metrics
  qc_violin <- VlnPlot(seurat_object,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       ncol = 3)
  
  # Saving the violin plots as png files for the analysis of the single cell sequencing data
  ggsave(file.path(output_dir, "01_qc_violin_plots.png"),
         plot = qc_violin, width = 12, height = 5)
  
  # Create feature scatter plots
  feature_scatter <- FeatureScatter(seurat_object,
                                    feature1 = "nCount_RNA",
                                    feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm")
  
  ggsave(file.path(output_dir, "01_qc_feature_scatter.png"),
         plot = feature_scatter, width = 8, height = 6)
  
  # Create mitochondrial percentage plot
  mt_violin <- VlnPlot(seurat_object, features = "percent.mt") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red")
  
  ggsave(file.path(output_dir, "01_qc_mitochondrial.png"),
         plot = mt_violin, width = 6, height = 6)
}