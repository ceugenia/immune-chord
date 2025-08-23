# Using crazyhottommy's ATAC-seq tools (Python focus)
import pyatac
from pyatac import ChIPSeq

# Download melanoma scRNA-seq + scATAC-seq data (GSE123139)
!prefetch SRX1234567 && fasterq-dump SRX1234567

# Process ATAC-seq data with ATAC-seq pipeline
atac = ChIPSeq("SRX1234567.fastq", genome="hg38")
atac.run_align() \
   .call_peaks(method="macs2") \
   .make_bigwig()

# Process scRNA-seq data with Scanpy
import scanpy as sc
adata = sc.read_10x_mtx("filtered_gene_bc_matrices")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.normalize_total(adata, target_sum=1e4)