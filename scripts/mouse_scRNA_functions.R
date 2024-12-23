# ------------------------------------------------------------------------------
# R Script for scRNA-Seq Data Preprocessing Functions
# ------------------------------------------------------------------------------
#
# Description:
# This script provides a set of reusable functions for preprocessing and 
# analyzing single-cell RNA sequencing (scRNA-seq) data. These functions are 
# designed to handle quality control, normalization, dimensionality reduction, 
# clustering, annotation, and visualization.
#
# ------------------------------------------------------------------------------
# Requirements:
# 
# 1. Install required libraries:
#    - `Seurat`, `ggplot2`, `dplyr`, `DoubletFinder`, `openxlsx`, `limma`,
#      `HGNChelper`, `cowplot`, `tidyr`, `pheatmap`, `RColorBrewer`, 
#      `SingleR`, `SingleCellExperiment`.
# 
# 2. Set the working environment:
#    - Load the required libraries at the beginning of the script.
#    - Ensure all input data and reference databases (e.g., ScType) are properly 
#      specified in your workflow scripts.
# 
# ------------------------------------------------------------------------------
# Provided Functions:
# 
# 1. **Preprocessing Functions**:
#    - `set_default_assay(seurat_obj)`: Sets the default assay to "RNA".
#    - `normalize_data(seurat_obj)`: Normalizes the data using `NormalizeData`.
#    - `find_variable_features(seurat_obj)`: Identifies variable features.
#    - `scale_data(seurat_obj)`: Scales the data for downstream analysis.
#    - `run_pca(seurat_obj, npcs)`: Performs PCA for dimensionality reduction.
#    - `find_neighbors(seurat_obj, dims)`: Identifies neighbors for clustering.
#    - `find_clusters(seurat_obj, resolution)`: Clusters the data.
#    - `run_umap(seurat_obj, dims)`: Runs UMAP for visualization.
#    - `apply_doubletfinder(seurat_obj)`: Identifies and removes doublets.
# 
# 2. **Quality Control Functions**:
#    - `plot_violin_qc_before(seurat_obj, base_folder)`: Generates QC violin 
#      plots before filtering.
#    - `plot_violin_qc_after(seurat_obj, base_folder)`: Generates QC violin 
#      plots after filtering.
#    - `filter_qc(seurat_obj)`: Filters low-quality cells based on metrics.
# 
# 3. **Annotation Functions**:
#    - Functions for ScType annotation and integration with metadata.
# 
# 4. **Visualization Functions**:
#    - `generate_gene_expression_umaps(seurat_objects, genes_of_interest, base_folder)`:
#      Generates UMAPs for gene expression.
#    - `generate_sctype_umap(seurat_objects, base_folder)`: Generates UMAPs with 
#      ScType annotations.
#    - `generate_combined_ifng_pdcd1_umap(seurat_objects, base_folder)`: Creates 
#      UMAPs for combined IFNG and PDCD1 expression.
#    - `generate_upset_plot_from_dataframe(data, genes, output_dir, dataset_name)`:
#      Generates UpSet plots for gene expression.
# 
# 5. **Export Functions**:
#    - `extract_gene_expression_and_save(seurat_obj, genes, base_folder, dataset_name)`:
#      Extracts gene expression data and annotations to save as CSV.
#    - `extract_and_save_top_genes_sctype(combined_seurat_object, top_n, output_dir)`:
#      Extracts top-expressed genes for each cell type.
# 
# 6. **Heatmap and Scatterplot Functions**:
#    - `generate_and_save_heatmaps(combined_data, genes_list, output_dir)`: 
#      Creates heatmaps for selected genes.
#    - `generate_and_save_scatterplots(seurat_obj, base_folder, genes_of_interest)`:
#      Generates scatterplots for pairs of genes of interest.
# 
# ------------------------------------------------------------------------------
# Usage:
# 
# - Source this script to include the functions in your analysis workflow.
# - Call the desired functions with appropriate parameters in your main script.
# 
# ------------------------------------------------------------------------------
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(openxlsx)
library(limma)
library(HGNChelper)
library(cowplot)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(cowplot) 
library(SingleR)
library(SingleCellExperiment)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")


### Preprocess Functions ###
# Set Default Assay
set_default_assay <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  return(seurat_obj)
}

# Normalize Data
normalize_data <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = TRUE)
  
  # Check if normalization was successful using slot presence checking in a more robust manner
  data_layer <- tryCatch({
    GetAssayData(seurat_obj, slot = "data")
  }, error = function(e) {
    stop("Normalization failed: 'data' slot not found. Please check if NormalizeData() ran successfully.")
  })
  
  return(seurat_obj)
}


# Find Variable Features
find_variable_features <- function(seurat_obj) {
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  return(seurat_obj)
}

# Scale Data
scale_data <- function(seurat_obj) {
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)
  return(seurat_obj)
}

# Perform PCA
run_pca <- function(seurat_obj, npcs = 50) {
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs, verbose = FALSE)
  return(seurat_obj)
}

# Find Neighbors
find_neighbors <- function(seurat_obj, dims = 1:10) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, verbose = FALSE)
  return(seurat_obj)
}

# Find Clusters
find_clusters <- function(seurat_obj, resolution = 0.5) {
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  return(seurat_obj)
}

# Run UMAP
run_umap <- function(seurat_obj, dims = 1:10) {
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, verbose = FALSE)
  return(seurat_obj)
}



### Apply DoubletFinder Function ###
apply_doubletfinder <- function(seurat_obj) {
  pcs <- 1:20
  pN <- 0.25
  pK <- 0.005
  nExp_poi <- round(0.075 * ncol(seurat_obj))
  
  seurat_obj <- DoubletFinder::doubletFinder(seurat_obj, PCs = pcs, pN = pN, pK = pK, nExp = nExp_poi)
  seurat_obj$DF.classifications <- seurat_obj@meta.data[,paste0("DF.classifications_", pN, "_", pK, "_", nExp_poi)]
  seurat_obj <- subset(seurat_obj, subset = DF.classifications == "Singlet")
  
  return(seurat_obj)
}

### ScType Annotation Functions ###

### Violin Plot Functions ###
plot_violin_qc_before <- function(seurat_obj, base_folder) {
  qc_folder <- file.path(base_folder, "QC")
  if (!dir.exists(qc_folder)) {
    dir.create(qc_folder, recursive = TRUE)
  }
  
  if (anyNA(seurat_obj@meta.data$nFeature_RNA) || anyNA(seurat_obj@meta.data$percent.mt) || anyNA(seurat_obj@meta.data$nCount_RNA)) {
    cat("NA values detected in QC metrics before filtering, removing NAs...\n")
    seurat_obj <- seurat_obj %>% subset(subset = !is.na(nFeature_RNA) & !is.na(percent.mt) & !is.na(nCount_RNA))
  }
  
  if (all(seurat_obj@meta.data$percent.mt == seurat_obj@meta.data$percent.mt[1])) {
    cat("Uniform percent.mt values detected. Skipping violin plot...\n")
  } else {
    p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt", "nCount_RNA"), group.by = "orig.ident", pt.size = 0.1) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 7, face = "bold")
      )
    
    ggsave(filename = file.path(qc_folder, paste0(seurat_obj@project.name, "_violin_plot_before_filter.png")), plot = p, dpi = 600, width = 12, height = 6, units = "in")
  }
}

plot_violin_qc_after <- function(seurat_obj, base_folder) {
  qc_folder <- file.path(base_folder, "QC")
  if (!dir.exists(qc_folder)) {
    dir.create(qc_folder, recursive = TRUE)
  }
  
  if (anyNA(seurat_obj@meta.data$nFeature_RNA) || anyNA(seurat_obj@meta.data$percent.mt) || anyNA(seurat_obj@meta.data$nCount_RNA)) {
    cat("NA values detected in QC metrics after filtering, removing NAs...\n")
    seurat_obj <- seurat_obj %>% subset(subset = !is.na(nFeature_RNA) & !is.na(percent.mt) & !is.na(nCount_RNA))
  }
  
  if (all(seurat_obj@meta.data$percent.mt == seurat_obj@meta.data$percent.mt[1])) {
    cat("Uniform percent.mt values detected. Skipping violin plot...\n")
  } else {
    p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt", "nCount_RNA"), group.by = "orig.ident", pt.size = 0.1) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 7, face = "bold")
      )
    
    ggsave(filename = file.path(qc_folder, paste0(seurat_obj@project.name, "_violin_plot_after_filter.png")), plot = p, dpi = 600, width = 12, height = 6, units = "in")
  }
}

### Filtering Function ###
filter_qc <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15 & nCount_RNA > 500)
  return(seurat_obj)
}

# Extract gene expression and annotations function
extract_gene_expression_and_annotation <- function(seurat_obj, genes, output_dir, dataset_name) {
  available_genes <- genes[genes %in% rownames(seurat_obj)]
  if (length(available_genes) == 0) {
    cat("No genes of interest found in the dataset:", dataset_name, "\n")
    return(NULL)
  }
  
  expression_data <- GetAssayData(seurat_obj, slot = "data")[available_genes, , drop = FALSE]
  annotations <- seurat_obj@meta.data$sctype_classification
  
  if (length(annotations) != ncol(expression_data)) {
    stop("Mismatch between number of annotations and number of cells in the Seurat object for dataset: ", dataset_name)
  }
  
  data_to_save <- data.frame(
    Barcode = colnames(expression_data),
    Annotation = annotations,
    t(as.matrix(expression_data))
  )
  
  output_file <- file.path(output_dir, paste0(dataset_name, "_IFNg-PD1_gene_expression.csv"))
  write.csv(data_to_save, file = output_file, row.names = FALSE)
  
  cat("Saved gene expression and annotation data to", output_file, "\n")
}

# Find top expressing genes for each ScType cluster
find_top_genes_per_sctype <- function(seurat_obj, top_n = 5) {
  DefaultAssay(seurat_obj) <- "RNA"
  avg_exp <- AverageExpression(seurat_obj, assays = "RNA", group.by = "sctype_classification", return.seurat = FALSE)$RNA
  
  top_genes <- as.data.frame(avg_exp) %>%
    rownames_to_column(var = "Gene") %>%
    gather(key = "Cluster", value = "Expression", -Gene) %>%
    group_by(Cluster) %>%
    top_n(n = top_n, wt = Expression) %>%
    arrange(Cluster, desc(Expression))
  
  return(top_genes)
}

# Extract and save top expressing genes for each dataset
extract_and_save_top_genes_sctype <- function(combined_seurat_object, top_n = 5, output_dir) {
  for (dataset_name in names(combined_seurat_object)) {
    cat("Processing dataset:", dataset_name, "\n")
    
    seurat_obj <- combined_seurat_object[[dataset_name]]
    top_genes <- find_top_genes_per_sctype(seurat_obj, top_n = top_n)
    
    output_file <- file.path(output_dir, paste0(dataset_name, "_top_", top_n, "_genes_per_sctype.csv"))
    write.csv(top_genes, file = output_file, row.names = FALSE)
    
    cat("Saved top expressing genes for", dataset_name, "to", output_file, "\n")
  }
}

# Function to extract gene expression and annotations and save them to CSV
extract_gene_expression_and_save <- function(seurat_obj, genes, base_folder, dataset_name) {
  # Check if the genes are present in the dataset
  available_genes <- genes[genes %in% rownames(seurat_obj)]
  if (length(available_genes) == 0) {
    cat("No genes of interest found in the dataset:", dataset_name, "\n")
    return(NULL)
  }
  
  # Extract the expression data for the specified genes
  expression_data <- GetAssayData(seurat_obj, slot = "data")[available_genes, , drop = FALSE]
  
  # Get annotations from metadata
  annotations <- seurat_obj@meta.data$sctype_classification
  
  # Create a data frame with barcodes, annotations, and expression levels
  data_to_save <- data.frame(
    Barcode = colnames(expression_data),
    Annotation = annotations,
    t(as.matrix(expression_data))
  )
  
  # Define the output file path
  output_file <- file.path(base_folder, paste0(dataset_name, "_IFNg_PD1_gene_expression.csv"))
  
  # Save the data frame as a CSV
  write.csv(data_to_save, file = output_file, row.names = FALSE)
  
  cat("Saved gene expression and annotation data to", output_file, "\n")
}

# Function to generate and save UpSet plots from gene expression data
generate_upset_plot_from_dataframe <- function(data, genes, output_dir, dataset_name) {
  # Ensure the genes of interest are present in the data
  available_genes <- genes[genes %in% colnames(data)]
  if (length(available_genes) == 0) {
    stop("None of the genes of interest are found in the dataset: ", dataset_name)
  }
  
  # Filter the data to include only the columns of interest (genes and barcodes)
  upset_data <- data %>%
    select(Barcode, all_of(available_genes)) %>%
    mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))  # Convert expression data to binary (1 = expressed, 0 = not expressed)
  
  # Create the UpSet plot
  upset_plot <- upset(upset_data, intersect = available_genes, base_annotations = list(
    'Intersection size' = intersection_size(counts = TRUE)
  )) + ggtitle(paste(dataset_name, "- UpSet Plot for Selected Genes"))
  
  # Define the output file path based on the genes
  output_file <- file.path(output_dir, paste0(dataset_name, "_upset_plot_", paste(genes, collapse = "_"), ".png"))
  
  # Save the plot
  ggsave(filename = output_file, plot = upset_plot, dpi = 600, width = 10, height = 8, units = "in")
  
  cat("UpSet plot saved for", dataset_name, "with genes", paste(genes, collapse = ", "), "\n")
}


# Function to generate and save heatmaps with specific genes, clustered by cell type (Annotation) and in specified gene order
generate_and_save_heatmaps <- function(combined_data, genes_list, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each dataset and generate heatmaps for each gene list
  for (genes in genes_list) {
    # Filter the data to include only the specified genes and ensure the correct column order
    filtered_data <- combined_data %>% filter(Annotation %in% unique(combined_data$Annotation))
    
    # Ensure the columns are in the specified order
    filtered_data <- filtered_data %>% select(Barcode, Annotation, Dataset, all_of(genes))
    
    # Get dataset names
    dataset_names <- unique(filtered_data$Dataset)
    
    for (dataset_name in dataset_names) {
      # Subset data for the current dataset
      dataset_data <- filtered_data %>% filter(Dataset == dataset_name)
      
      # Handle duplicate barcodes by adding a suffix to make them unique
      dataset_data$Barcode <- make.unique(as.character(dataset_data$Barcode))
      
      # Create a matrix for the heatmap
      heatmap_data <- dataset_data %>% select(all_of(genes))  # Explicitly ensure column order
      rownames(heatmap_data) <- dataset_data$Barcode
      
      # Create a row annotation for the heatmap
      annotation_row <- data.frame(Annotation = dataset_data$Annotation)
      rownames(annotation_row) <- dataset_data$Barcode
      
      # Generate annotation colors for the levels in `Annotation`
      unique_annotations <- unique(dataset_data$Annotation)
      annotation_colors <- setNames(brewer.pal(n = length(unique_annotations), "Set3"), unique_annotations)
      
      # Ensure the levels of `Annotation` match the keys in `annotation_colors`
      annotation_row$Annotation <- factor(annotation_row$Annotation, levels = unique_annotations)
      
      # Sort cells within each annotation by the expression of the first gene (e.g., IFNG)
      ordered_barcodes <- unlist(lapply(unique_annotations, function(cell_type) {
        cells_in_group <- dataset_data %>% filter(Annotation == cell_type)
        cells_in_group <- cells_in_group %>% arrange(desc(!!sym(genes[1])))  # Sort by the expression of the first gene
        return(cells_in_group$Barcode)
      }))
      
      # Reorder the heatmap data and annotations based on the sorted barcodes
      heatmap_data <- heatmap_data[ordered_barcodes, ]
      annotation_row <- annotation_row[ordered_barcodes, , drop = FALSE]
      
      # Generate the heatmap
      heatmap_plot <- pheatmap(
        heatmap_data,
        cluster_rows = FALSE, # Rows are manually ordered by annotation and expression
        cluster_cols = FALSE,  # Disable clustering for columns (preserve gene order)
        show_rownames = FALSE,
        annotation_row = annotation_row,
        main = paste(dataset_name, "- Heatmap for Selected Genes"),
        annotation_colors = list(Annotation = annotation_colors)
      )
      
      # Save the heatmap
      ggsave(
        filename = file.path(output_dir, paste0(dataset_name, "_heatmap_", paste(genes, collapse = "_"), ".png")),
        plot = heatmap_plot,
        dpi = 600,
        width = 10,
        height = 8,
        units = "in"
      )
      
      cat("Saved heatmap for", dataset_name, "with genes", paste(genes, collapse = ", "), "\n")
    }
  }
}


# Function to generate and save scatterplots correlating genes of interest
generate_and_save_scatterplots <- function(seurat_obj, base_folder, genes_of_interest) {
  # Define the output directory for the scatterplots
  output_dir <- file.path(base_folder, "scatterplots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if at least two genes of interest are present in the dataset
  available_genes <- genes_of_interest[genes_of_interest %in% rownames(seurat_obj)]
  if (length(available_genes) < 2) {
    cat("Not enough genes of interest found in the dataset:", seurat_obj@project.name, "\n")
    return(NULL)
  }
  
  # Loop through each unique cell type annotation
  unique_annotations <- unique(seurat_obj@meta.data$sctype_classification)
  
  for (annotation in unique_annotations) {
    # Subset the data for the current annotation
    subset_data <- subset(seurat_obj, subset = sctype_classification == annotation)
    
    # Extract the expression data for the available genes of interest
    expression_data <- GetAssayData(subset_data, slot = "data")[available_genes, , drop = FALSE]
    
    # Generate scatterplots for all possible pairs of the available genes of interest
    for (i in 1:(length(available_genes) - 1)) {
      for (j in (i + 1):length(available_genes)) {
        gene_x <- available_genes[i]
        gene_y <- available_genes[j]
        
        # Create scatterplot
        scatter_plot <- ggplot(as.data.frame(t(expression_data)), aes_string(x = gene_x, y = gene_y)) +
          geom_point(color = "blue", alpha = 0.5) +
          geom_smooth(method = "lm", color = "red", se = FALSE) +
          ggtitle(paste(gene_x, "vs", gene_y, "-", annotation)) +
          xlab(paste(gene_x, "Expression")) +
          ylab(paste(gene_y, "Expression")) +
          theme_minimal()
        
        # Save the scatterplot
        ggsave(
          filename = file.path(output_dir, paste0(seurat_obj@project.name, "_", annotation, "_", gene_x, "_vs_", gene_y, ".png")),
          plot = scatter_plot, dpi = 600, width = 8, height = 6, units = "in"
        )
        
        cat("Saved scatterplot for", annotation, "with genes", gene_x, "and", gene_y, "\n")
      }
    }
  }
}


# Function to generate and save UpSet plots for all cell types
generate_and_save_upset_plots_cell_types <- function(seurat_obj, genes_of_interest, base_folder) {
  # Define the output directory for the UpSet plots
  output_dir <- file.path(base_folder, "upset_plots_cell_types")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each unique cell type annotation
  unique_annotations <- unique(seurat_obj@meta.data$sctype_classification)
  
  for (annotation in unique_annotations) {
    # Subset the data for the current annotation
    subset_data <- subset(seurat_obj, subset = sctype_classification == annotation)
    
    # Check if all genes of interest are present in the dataset
    available_genes <- genes_of_interest[genes_of_interest %in% rownames(seurat_obj)]
    if (length(available_genes) < 2) {
      cat("Not enough genes of interest found for UpSet plot in annotation:", annotation, "\n")
      next
    }
    
    # Extract the expression data for the genes of interest
    expression_data <- GetAssayData(subset_data, slot = "data")[available_genes, , drop = FALSE]
    
    # Convert expression data to a binary matrix (1 = expressed, 0 = not expressed)
    binary_data <- as.data.frame(t(expression_data)) %>%
      mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))
    
    # Generate the UpSet plot
    upset_plot <- upset(binary_data, intersect = available_genes, base_annotations = list(
      'Intersection size' = intersection_size(counts = TRUE)
    )) + ggtitle(paste(annotation, "- UpSet Plot"))
    
    # Save the UpSet plot
    ggsave(filename = file.path(output_dir, paste0(seurat_obj@project.name, "_", annotation, "_upset_plot.png")), plot = upset_plot, dpi = 600, width = 10, height = 8, units = "in")
    
    cat("Saved UpSet plot for", annotation, "to", output_dir, "\n")
  }
}


# Function to generate UMAP plots for each gene of interest, showing expression levels
generate_gene_expression_umaps <- function(seurat_objects, genes_of_interest, base_folder) {
  umap_folder <- file.path(base_folder, "UMAPs")
  if (!dir.exists(umap_folder)) {
    dir.create(umap_folder, recursive = TRUE)
  }
  
  for (sample_name in names(seurat_objects)) {
    # Set the default assay to "RNA"
    DefaultAssay(seurat_objects[[sample_name]]) <- "RNA"
    
    # Get the UMAP coordinates for the cells
    umap_coords <- Embeddings(seurat_objects[[sample_name]], "umap")
    
    # Get the gene names from the RNA assay
    gene_names <- rownames(seurat_objects[[sample_name]][["RNA"]])
    
    # Loop through each gene of interest and create a UMAP plot
    for (gene in genes_of_interest) {
      if (gene %in% gene_names) {
        # Get expression data for the current gene
        gene_expression <- GetAssayData(seurat_objects[[sample_name]], slot = "data")[gene, ]
        
        # Create a data frame for UMAP coordinates and gene expression
        umap_gene <- data.frame(UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2], Expression = gene_expression)
        
        # Plot non-expressing cells in light gray and expressing cells on top with color gradient
        p_gene <- ggplot() +
          geom_point(data = umap_gene %>% filter(Expression == 0), 
                     aes(x = UMAP_1, y = UMAP_2), color = "lightgray", size = 0.5) +
          geom_point(data = umap_gene %>% filter(Expression > 0), 
                     aes(x = UMAP_1, y = UMAP_2, color = Expression), size = 0.5) +
          scale_color_gradientn(colors = c("blue", "yellow", "red"),
                                name = "Expression Level") +
          theme_minimal() +
          ggtitle(paste(sample_name, "-", gene, "Expression")) +
          theme(plot.title = element_text(hjust = 0.5))
        
        # Save the plot as PNG
        ggsave(filename = file.path(umap_folder, paste0(sample_name, "_", gene, "_expression_umap.png")),
               plot = p_gene, dpi = 600, width = 8, height = 6, units = "in")
        
        # Save the plot as PDF
        ggsave(filename = file.path(umap_folder, paste0(sample_name, "_", gene, "_expression_umap.pdf")),
               plot = p_gene, width = 8, height = 6, units = "in")
        
      } else {
        cat("Gene", gene, "not found in dataset", sample_name, "\n")
      }
    }
  }
}



# Function to generate UMAP plots with ScType Annotations and consistent axis limits
generate_sctype_umap <- function(seurat_objects, base_folder) {
  umap_folder <- file.path(base_folder, "UMAPs")
  if (!dir.exists(umap_folder)) {
    dir.create(umap_folder, recursive = TRUE)
  }
  
  for (sample_name in names(seurat_objects)) {
    # Set the default assay to "RNA"
    DefaultAssay(seurat_objects[[sample_name]]) <- "RNA"
    
    # Get the UMAP coordinates for the cells
    umap_coords <- Embeddings(seurat_objects[[sample_name]], "umap")
    
    # Calculate axis limits for consistent scaling
    x_limits <- range(umap_coords[, 1])
    y_limits <- range(umap_coords[, 2])
    
    # Generate UMAP plot with ScType annotations
    umap_plot <- DimPlot(seurat_objects[[sample_name]], reduction = "umap", label = TRUE, group.by = 'sctype_classification') +
      xlim(x_limits) +  # Apply consistent x axis limits
      ylim(y_limits) +  # Apply consistent y axis limits
      ggtitle(paste(sample_name, "UMAP with ScType Annotations")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Correct output path to save in UMAPs folder
    output_path <- file.path(umap_folder, paste0(sample_name, "_UMAP_with_ScType_Annotations.png"))
    
    # Save the plot
    ggsave(filename = output_path, plot = umap_plot, dpi = 600, width = 8, height = 6, units = "in")
    
    cat("Finished processing ScType UMAP for dataset:", sample_name, "\n")
  }
}

generate_combined_ifng_pdcd1_umap <- function(seurat_objects, base_folder) {
  umap_folder <- file.path(base_folder, "UMAPs")
  if (!dir.exists(umap_folder)) {
    dir.create(umap_folder, recursive = TRUE)
  }
  
  for (sample_name in names(seurat_objects)) {
    # Set the default assay to "RNA"
    DefaultAssay(seurat_objects[[sample_name]]) <- "RNA"
    
    # Get the UMAP coordinates for the cells
    umap_coords <- Embeddings(seurat_objects[[sample_name]], "umap")
    
    # Get expression data for IFNG and PDCD1
    ifng_expression <- GetAssayData(seurat_objects[[sample_name]], slot = "data")["IFNG", ]
    pdcd1_expression <- GetAssayData(seurat_objects[[sample_name]], slot = "data")["PDCD1", ]
    
    # Create a data frame for plotting
    umap_data <- data.frame(
      UMAP_1 = umap_coords[, 1],
      UMAP_2 = umap_coords[, 2],
      IFNG_Expression = ifng_expression,
      PDCD1_Expression = pdcd1_expression
    )
    
    # Determine the color for each cell
    umap_data$Expression <- ifelse(umap_data$IFNG_Expression > 0 & umap_data$PDCD1_Expression > 0, "Both",
                                   ifelse(umap_data$IFNG_Expression > 0, "IFNG",
                                          ifelse(umap_data$PDCD1_Expression > 0, "PDCD1", "None")))
    
    # Calculate percentages
    total_cells <- nrow(umap_data)
    percentage_ifng <- round(sum(umap_data$Expression == "IFNG") / total_cells * 100, 2)
    percentage_pdcd1 <- round(sum(umap_data$Expression == "PDCD1") / total_cells * 100, 2)
    percentage_both <- round(sum(umap_data$Expression == "Both") / total_cells * 100, 2)
    percentage_none <- round(sum(umap_data$Expression == "None") / total_cells * 100, 2)
    
    # Create custom legend labels including the percentages
    custom_labels <- c(paste0("IFNG (", percentage_ifng, "%)"),
                       paste0("PDCD1 (", percentage_pdcd1, "%)"),
                       paste0("Both (", percentage_both, "%)"),
                       paste0("None (", percentage_none, "%)"))
    
    # Generate the plot with the non-expressing cells (base layer) and add the expressing ones on top
    p_layered <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
      # Base layer: non-expressing cells
      geom_point(data = umap_data %>% filter(Expression == "None"), color = "white", size = 0.5) +
      # Expressing cells: IFNG, PDCD1, and both
      geom_point(data = umap_data %>% filter(Expression != "None"), aes(color = Expression), size = 1) +
      # Set custom colors for IFNG, PDCD1, both, and none
      scale_color_manual(values = c("IFNG" = "red", "PDCD1" = "blue", "Both" = "green", "None" = "white"),
                         breaks = c("IFNG", "PDCD1", "Both", "None"),
                         labels = custom_labels) +
      # Add title and minimal theme
      theme_minimal() +
      ggtitle(paste(sample_name, "UMAP - IFNG and PDCD1 Expression")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      # Force the legend to show and adjust legend
      guides(color = guide_legend(override.aes = list(size = 5))) +
      theme(legend.position = "right", legend.text = element_text(size = 10))
    
    # Save the plot
    ggsave(
      filename = file.path(umap_folder, paste0(sample_name, "_combined_IFNG_PDCD1_umap_with_legend.png")),
      plot = p_layered, dpi = 600, width = 8, height = 6, units = "in"
    )
    
    cat("Combined UMAP with layered expression and percentages for IFNG and PDCD1 saved for dataset:", sample_name, "\n")
  }
}
