# Diaz Lab Mouse scRNA Pipeline

## Overview
This repository contains R scripts for performing single-cell RNA sequencing (scRNA-seq) analysis of mouse data. The pipeline processes raw gene expression data, identifies cell populations, annotates cell types, and generates multiple visualizations to support biological insights. After the initial automatic annotation, a second semi-automatic annotation is performed using insights from the first, providing higher accuracy and greater flexibility. The scripts are suitable for use on small-scale datasets locally or on an HPC cluster for larger datasets.

## Features
### Data Preprocessing
- Imports raw scRNA-seq count matrices from the 10x Genomics Cellranger pipeline.
- Filters low-confidence barcodes, ambient RNA, and duplicate UMIs to retain high-quality data.
- Imports filtered UMI count matrices into Seurat objects for additional quality control and filtering.
- Applies mitochondrial and ribosomal RNA content thresholds, gene count, and UMI count thresholds for droplet filtering.
- Implements Principal Component Analysis (PCA) for dimensionality reduction.
- Visualizes clusters using Uniform Manifold Approximation and Projection (UMAP).
- Utilizes DoubletFinder to identify and exclude droplets containing more than one cell.

### Clustering and Cell Type Annotation
- Constructs k-nearest neighbor (kNN) graphs and clusters cells using the Louvain algorithm.
- Annotates cell types using SingleR and the MouseRNAseqData reference dataset.
- Generates CSV files showing top markers per cluster for refining annotations.
- Performs semi-automatic annotation with scType for greater accuracy and flexibility.
- Fully annotated Seurat objects are saved as outputs for downstream analysis.

### Differential Gene Expression (DGE)
- Integrates datasets across samples or batches using Seurat's anchor-based workflow.
- Identifies upregulated and downregulated genes under different conditions.
- Supports multiple comparisons in a single run.

### Gene Ontology and Pathway Enrichment Analysis
- Conducts Gene Set Enrichment Analysis (GSEA) using the clusterProfiler package.
- Analyzes Gene Ontology (GO) terms and KEGG pathways for biological insights.
- Maps gene identifiers to ENSEMBL and ENTREZID formats and ranks genes by log-fold change.

## Software
- **R (version 4.0 or higher)**
- **RStudio (recommended)**

### Required R Packages:
- Seurat: 5.1.0
- ggplot2: 3.5.1
- dplyr: 1.1.4
- DoubletFinder: 2.0.4
- openxlsx: 4.2.7.1
- limma: 3.60.6
- HGNChelper: 0.8.14
- cowplot: 1.1.3
- tidyr: 1.3.1
- pheatmap: 1.0.12
- RColorBrewer: 1.1.3
- SingleR: 2.6.0
- SingleCellExperiment: 1.26.0

## Usage
### Order of Execution:
1. `mouse_first_annotation.R`
2. `mouse_gen_seurat_object.R`
3. `mouse_DEGs_GSEA.R`

### Input Files
- **Gene Expression Matrix:** Place the following files (output of Cellranger) in a directory named after the sample:
  - `barcodes.tsv.gz`
  - `features.tsv.gz`
  - `matrix.mtx.gz`

### Configuration
Edit the script to adjust parameters:
- Cell and gene filtering thresholds.
- UMAP and clustering parameters.
- Paths for input and output files.
- Specific genes of interest.

## Customization
- **Cell Filtering:** Adjust `min.cells` and `min.features` for filtering low-quality cells and genes.
- **Clustering Resolution:** Modify the resolution parameter in `Seurat::FindClusters()`.
- **Cell Type Annotation:** Update ScType references as needed.

## Troubleshooting
- Use a high-performance compute cluster for large datasets.
- Ensure all required R packages are installed.

## Contributions
Feel free to fork this repository and submit pull requests for bug fixes or new features.

## License
This project is licensed under the Creative Commons Attribution 4.0 International.

## Contact
**Name:** Georgios Vlachos  
**Email:** vlachog1@mskcc.com

