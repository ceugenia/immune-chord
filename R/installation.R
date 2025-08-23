# 1. Install CRAN packages first (most stable)
install.packages(c(
  "BiocManager",
  "remotes",
  "numDeriv",
  "ggrepel",
  "msigdbr",
  "clusterProfiler",
  "dplyr",
  "patchwork",
  "viridis",
  "hdf5r"
))

# Enter commands in R (or R studio, if installed)
install.packages('Seurat')

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE) #1

# 2. Install specific Seurat versions (as requested)
# Note: This might cause dependency conflicts with newer Bioconductor packages
#remotes::install_version("SeuratObject", "4.1.4", 
                         #repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", 
                         #repos = c("https://satijalab.r-universe.dev", getOption("repos"))) #choose option 3 - None

# 3. Install Bioconductor packages
library(BiocManager)
# Install a compatible Bioconductor version
#BiocManager::install(version = "3.18")  # Adjust based on your R version

BiocManager::install(c(
  "SingleCellExperiment",
  "scran",
  "scater",
  "ComplexHeatmap",
  "limma",
  "edgeR",
  "DESeq2",
  "zellkonverter"
))

# 4. Install GitHub packages
# For SeuratDisk (if not available on CRAN)
remotes::install_github("mojaveazure/seurat-disk")

# For CIARA
install.packages("CIARA")

# First, install devtools if you haven't already
install.packages("devtools")

# Install GiniClust2 from your local directory

# Load the package
#library(GiniClust2)

# 5. Load all packages to verify installation
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(scran)
library(scater)
library(ComplexHeatmap)
library(zellkonverter)
library(SeuratDisk)

# 6. Check versions
packageVersion("Seurat")
packageVersion("SingleCellExperiment")

# Test basic functionality
library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

# Test data loading
data("pbmc_small")
pbmc_small

# Test H5AD reading capability
library(SeuratDisk)

# Install potential dependencies
install.packages(c("Rcpp", "RcppEigen", "Matrix", "ggplot2", "gplots", "stats", "utils", "graphics", "grDevices", "methods"))

# Also install Bioconductor dependencies if needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scran", "scater"))

# Install with force = TRUE to ensure fresh installation
devtools::install_github("landerlabcode/BigSurR", force = TRUE)
