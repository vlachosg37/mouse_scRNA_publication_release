# ------------------------------------------------------------------------------
# R Script for Preliminary Cell Type Annotation of Mouse scRNA Data
# ------------------------------------------------------------------------------
#
# Description:
# This script performs preliminary cell type annotation for single-cell RNA 
# sequencing (scRNA-seq) data from mouse samples. It integrates quality control, 
# dimensionality reduction, clustering, and ScType-based cell type annotation. 
# The annotated data can be refined in subsequent analyses.
#
# ------------------------------------------------------------------------------
# Requirements:
# 
# 1. Set the path to required scripts:
#    - Specify the path to the `mouse_scRNA_functions.R` script containing 
#      helper functions for preprocessing and visualization.
# 
# 2. Input data:
#    - Provide input scRNA-seq data directories with filtered feature-barcode 
#      matrices (output of Cell Ranger or similar tools).
# 
# 3. Set the ScType database:
#    - Specify the path to the ScType database for initial cell type annotation.
#    - Use the results of this annotation (_top_genes_per_cluster.csv) to refine 
#      and adjust the database for higher accuracy in downstream analyses.
# 
# 4. Specify the output folder:
#    - Define the base folder where all results (UMAPs, spreadsheets, heatmaps) 
#      will be saved.
# 
# 5. Genes of interest:
#    - Set a list of genes of interest for generating gene expression UMAPs, 
#      spreadsheets, and heatmaps.
# 
# ------------------------------------------------------------------------------
# Provided Features:
#
# 1. **Preprocessing**:
#    - Reads filtered feature-barcode matrices and creates Seurat objects.
#    - Performs quality control to add mitochondrial and ribosomal content 
#      metrics.
# 
# 2. **Dimensionality Reduction and Clustering**:
#    - Runs PCA, UMAP, and clustering algorithms to identify cell populations.
# 
# 3. **Cell Type Annotation**:
#    - Utilizes the ScType and SingleR framework for preliminary cell type annotation.
#    - Saves annotation results in Seurat objects and spreadsheets.
# 
# 4. **Gene Expression Visualization**:
#    - Generates UMAPs for genes of interest with expression-level coloring.
#    - Creates heatmaps to display the expression of selected genes.
#    - Exports gene expression data and annotations as spreadsheets.
# 
# ------------------------------------------------------------------------------
# Outputs:
#
# - Annotated Seurat objects saved in the specified output folder.
# - Gene expression UMAPs saved as PNG and PDF files.
# - Heatmaps of selected genes saved as PNG files.
# - Spreadsheets containing annotations and gene expression data for selected 
#   genes and clusters.
# 
# ------------------------------------------------------------------------------
# Usage:
# 
# - Adjust input file paths and parameters at the beginning of the script.
# - Source the required helper functions from `mouse_bothannot_functions.R`.
# - Run the script to generate preliminary annotations and visualizations.
#
# ------------------------------------------------------------------------------

# Load required libraries and functions
source("/Users/vlachog1/Desktop/scRNA_mice/241022_analysis/Scripts/mouse_bothannot_functions.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")

### Prepare gene sets from the ScType database ###
db_path <- "/Users/vlachog1/Desktop/scRNA_mice/ScType_databases/ScTypeDB_full_mouse.xlsx"
gs_list <- gene_sets_prepare(db_path, "Mouse")

### Define Genes of Interest ###
genes_of_interest <- c("Msh2", "Ifng", "Cd274", "Pdcd1", "Gzma", "Prf1")

### Set Base Folder ###
base_folder <- "/Users/vlachog1/Desktop/scRNA_mice/output"

# Check if the base folder exists, and create it if it doesn't
if (!dir.exists(base_folder)) {
  dir.create(base_folder, recursive = TRUE)  # recursive = TRUE creates any necessary parent directories as well
  cat("Base folder created:", base_folder, "\n")
} else {
  cat("Base folder already exists:", base_folder, "\n")
}
setwd(base_folder)

### Load and preprocess the data ###
# Specify directories for one or multiple datasets
data_dirs <- list(
  "/Users/vlachog1/Desktop/scRNA_mice/data/filtered/CDDP",
  "/Users/vlachog1/Desktop/scRNA_mice/data/filtered/TMZ",
  "/Users/vlachog1/Desktop/scRNA_mice/data/filtered/Parental",
  "/Users/vlachog1/Desktop/scRNA_mice/data/filtered/Combo"
)

# Initialize lists to store Seurat objects
seurat_objects <- list()

# Load data, create Seurat objects
for (data_dir in data_dirs) {
  sample_name <- basename(data_dir)
  data <- Read10X(data.dir = data_dir)
  
  if (is.list(data) && "Gene Expression" %in% names(data)) {
    rna_counts <- data$`Gene Expression`
  } else if (inherits(data, "dgCMatrix")) {
    rna_counts <- data
  } else {
    stop("Unexpected data type from Read10X")
  }
  
  # Create the Seurat object
  seurat_obj <- CreateSeuratObject(counts = rna_counts, project = sample_name)
  
  # Store the created Seurat object in the list
  seurat_objects[[sample_name]] <- seurat_obj
  
  # Print basic information
  cat("Created Seurat object for:", sample_name, "\n")
}

# Print summary of created Seurat objects
print(lapply(seurat_objects, summary))

# Add mitochondrial and ribosomal RNA content
for (sample_name in names(seurat_objects)) {
  DefaultAssay(seurat_objects[[sample_name]]) <- "RNA"
  seurat_objects[[sample_name]][["percent.mt"]] <- PercentageFeatureSet(seurat_objects[[sample_name]], pattern = "^MT-")
  seurat_objects[[sample_name]][["percent.rb"]] <- PercentageFeatureSet(seurat_objects[[sample_name]], pattern = "^RP[SL]")
  
  cat("Sample:", sample_name, "- percent.mt and percent.rb added\n")
}

###Preprocess, DoubletFinder, and Filtering ###
for (sample_name in names(seurat_objects)) {
  cat("Processing sample:", sample_name, "\n")
  
  # Step 1: Set Default Assay
  seurat_objects[[sample_name]] <- set_default_assay(seurat_objects[[sample_name]])
  
  # Step 2: Normalize Data
  seurat_objects[[sample_name]] <- normalize_data(seurat_objects[[sample_name]])
  
  # Step 3: Find Variable Features
  seurat_objects[[sample_name]] <- find_variable_features(seurat_objects[[sample_name]])
  
  # Step 4: Scale Data
  seurat_objects[[sample_name]] <- scale_data(seurat_objects[[sample_name]])
  
  # Step 5: Run PCA
  seurat_objects[[sample_name]] <- run_pca(seurat_objects[[sample_name]], npcs = 50)
  
  # Step 6: Find Neighbors
  seurat_objects[[sample_name]] <- find_neighbors(seurat_objects[[sample_name]], dims = 1:10)
  
  # Step 7: Find Clusters
  seurat_objects[[sample_name]] <- find_clusters(seurat_objects[[sample_name]], resolution = 0.5)
  
  # Step 8: Run UMAP
  seurat_objects[[sample_name]] <- run_umap(seurat_objects[[sample_name]], dims = 1:10)
  
  # Save filtered Seurat object
  saveRDS(seurat_objects[[sample_name]], file = file.path(base_folder, paste0(sample_name, "_processed_seurat.rds")))
}

### Visualize UMAPs and gene expression UMAPs ###
umap_folder <- file.path(base_folder, "UMAPs")
if (!dir.exists(umap_folder)) {
  dir.create(umap_folder, recursive = TRUE)
}

for (sample_name in names(seurat_objects)) {
  # Plot and save the RNA UMAP
  p1 <- DimPlot(seurat_objects[[sample_name]], reduction = "umap", label = TRUE, pt.size = 0.5) + 
    NoLegend() + 
    ggtitle(paste(sample_name, "- RNA"))
  
  ggsave(filename = file.path(umap_folder, paste0(sample_name, "_umap_rna.png")), plot = p1, dpi = 600, width = 8, height = 6, units = "in")
}

#Generate UMAP Plots for Gene Expression
generate_gene_expression_umaps(seurat_objects, genes_of_interest, base_folder)

# Loop through each Seurat object in the list
for (dataset_name in names(seurat_objects)) {
  cat("Processing dataset:", dataset_name, "\n")
  
  # Get the current Seurat object
  seurat_obj <- seurat_objects[[dataset_name]]
  
  # Set the default assay to "RNA"
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
  seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[["RNA"]])))
  print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))
  
  # Extract scaled scRNA-seq matrix
  scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(seurat_obj[["RNA"]]$scale.data) else as.matrix(seurat_obj[["RNA"]]@scale.data)
  
  # Run ScType
  es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # Merge by cluster
  cL_resutls <- do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters == cl)), 10)
  }))
  
  sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # Set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
  print(sctype_scores[,1:3])
  
  # Annotate the Seurat object (ensure 'sctype_classification' column is properly set in metadata)
  seurat_obj@meta.data$sctype_classification <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    seurat_obj@meta.data$sctype_classification[seurat_obj@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
  }
  
  # Save the annotated Seurat object back into the list
  seurat_objects[[dataset_name]] <- seurat_obj
  
  # Generate UMAP plot with labels
  umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification') +
    ggtitle(paste(dataset_name, "UMAP with ScType Annotations"))
  
  # Generate UMAP plot without labels
  umap_plot_no_labels <- DimPlot(seurat_obj, reduction = "umap", label = FALSE, group.by = 'sctype_classification') +
    ggtitle(paste(dataset_name, "UMAP with ScType Annotations"))
  
  # Define output paths for labeled and no-labels versions
  output_path_labeled_png <- file.path(umap_folder, paste0(dataset_name, "_UMAP_with_ScType_Annotations.png"))
  output_path_labeled_pdf <- file.path(umap_folder, paste0(dataset_name, "_UMAP_with_ScType_Annotations.pdf"))
  output_path_no_labels_png <- file.path(umap_folder, paste0(dataset_name, "_UMAP_with_ScType_Annotations_no_labels.png"))
  output_path_no_labels_pdf <- file.path(umap_folder, paste0(dataset_name, "_UMAP_with_ScType_Annotations_no_labels.pdf"))
  
  # Save labeled version
  ggsave(filename = output_path_labeled_png, plot = umap_plot, dpi = 600, width = 8, height = 6, units = "in")
  ggsave(filename = output_path_labeled_pdf, plot = umap_plot, width = 8, height = 6, units = "in")
  
  # Save no-labels version
  ggsave(filename = output_path_no_labels_png, plot = umap_plot_no_labels, dpi = 600, width = 8, height = 6, units = "in")
  ggsave(filename = output_path_no_labels_pdf, plot = umap_plot_no_labels, width = 8, height = 6, units = "in")
  
  cat("Finished processing dataset:", dataset_name, "\n")
}

cat("All datasets processed and UMAP plots saved.\n")

### Save the gene expression CSV in the "Spreadsheets" folder ###
# Create the "Spreadsheets" folder inside the base folder
spreadsheets_folder <- file.path(base_folder, "Spreadsheets")
if (!dir.exists(spreadsheets_folder)) {
  dir.create(spreadsheets_folder, recursive = TRUE)
}

# Initialize an empty data frame to store the combined data
combined_gene_expression_data <- data.frame()

# Loop through each dataset and save gene expression data to CSV
for (dataset_name in names(seurat_objects)) {
  cat("Processing dataset:", dataset_name, "\n")
  
  # Extract the Seurat object for the current dataset
  seurat_obj <- seurat_objects[[dataset_name]]
  
  # Check if the genes are present in the dataset
  available_genes <- genes_of_interest[genes_of_interest %in% rownames(seurat_obj)]
  if (length(available_genes) == 0) {
    cat("No genes of interest found in the dataset:", dataset_name, "\n")
    next
  }
  
  # Extract the expression data for the specified genes
  expression_data <- GetAssayData(seurat_obj, slot = "data")[available_genes, , drop = FALSE]
  
  # Get annotations from metadata
  annotations <- seurat_obj@meta.data$sctype_classification
  
  # Create a data frame with barcodes, annotations, expression levels, and dataset name
  data_to_save <- data.frame(
    Barcode = colnames(expression_data),
    Annotation = annotations,
    t(as.matrix(expression_data)),
    Dataset = dataset_name  # Add the dataset name
  )
  
  # Append the current dataset's data to the combined data frame
  combined_gene_expression_data <- rbind(combined_gene_expression_data, data_to_save)
  
  # Define the output file path for the individual dataset
  output_file <- file.path(spreadsheets_folder, paste0(dataset_name, "_genes_of_interest_expression.csv"))
  
  # Save the data frame as a CSV for the individual dataset
  write.csv(data_to_save, file = output_file, row.names = FALSE)
  
  cat("Saved gene expression and annotation data to", output_file, "\n")
}

# Define the output file path for the combined data
combined_output_file <- file.path(spreadsheets_folder, paste0("combined_genes_of_interest_expression.csv"))

# Save the combined data frame as a CSV
write.csv(combined_gene_expression_data, file = combined_output_file, row.names = FALSE)

cat("Combined gene expression and annotation data saved to", combined_output_file, "\n")



# Define the output directory for the heatmaps
output_dir <- file.path(base_folder, "Heatmaps")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set the correct path for the combined CSV file
combined_output_file <- file.path(spreadsheets_folder, "combined_genes_of_interest_expression.csv")

# Load the combined data from the CSV file
combined_data <- read.csv(combined_output_file)

# Generate and save heatmaps using the provided genes of interest
generate_and_save_heatmaps(combined_data, list(genes_of_interest), output_dir)

# Convert Seurat objects to SingleCellExperiment and perform SingleR annotation with MouseRNAseqData
# Load reference dataset
mouse_rnaseq_ref <- celldex::MouseRNAseqData()

# Annotate all Seurat objects using MouseRNAseqData and save annotations to the metadata
for (sample_name in names(seurat_objects)) {
  # Set default assay to RNA
  DefaultAssay(seurat_objects[[sample_name]]) <- "RNA"
  
  # Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(DietSeurat(seurat_objects[[sample_name]]))
  
  # Perform SingleR annotation using MouseRNAseqData reference
  mouse_main <- SingleR(test = sce, assay.type.test = 1, ref = mouse_rnaseq_ref, labels = mouse_rnaseq_ref$label.main)
  mouse_fine <- SingleR(test = sce, assay.type.test = 1, ref = mouse_rnaseq_ref, labels = mouse_rnaseq_ref$label.fine)
  
  # Add MouseRNAseqData annotations to the Seurat object metadata in separate columns
  seurat_objects[[sample_name]]@meta.data$mouse_main <- mouse_main$pruned.labels
  seurat_objects[[sample_name]]@meta.data$mouse_fine <- mouse_fine$pruned.labels
  
  # Save the updated Seurat object back to the list
  seurat_objects[[sample_name]] <- seurat_objects[[sample_name]]
  
  # Plot UMAP with MouseRNAseqData fine annotations and save it
  umap_plot <- DimPlot(seurat_objects[[sample_name]], reduction = "umap", group.by = "mouse_fine", label = TRUE, repel = TRUE) +
    ggtitle(paste(sample_name, "UMAP with MouseRNAseqData Fine Annotations"))
  
  ggsave(filename = file.path(umap_folder, paste0(sample_name, "_UMAP_with_MouseRNAseq_Fine_Annotations.png")), plot = umap_plot, dpi = 600, width = 8, height = 6, units = "in")
  
  # Optionally, plot UMAP with MouseRNAseqData main annotations and save it
  umap_plot_main <- DimPlot(seurat_objects[[sample_name]], reduction = "umap", group.by = "mouse_main", label = TRUE, repel = TRUE) +
    ggtitle(paste(sample_name, "UMAP with MouseRNAseqData Main Annotations"))
  
  ggsave(filename = file.path(umap_folder, paste0(sample_name, "_UMAP_with_MouseRNAseq_Main_Annotations.png")), plot = umap_plot_main, dpi = 600, width = 8, height = 6, units = "in")
  
  cat("MouseRNAseqData annotations added and UMAP plots saved for dataset:", sample_name, "\n")
}

cat("All datasets processed and MouseRNAseqData annotations saved to the Seurat objects.\n")

# Loop through each Seurat object and save the combined data to CSV
for (dataset_name in names(seurat_objects)) {
  cat("Processing dataset:", dataset_name, "\n")
  
  # Extract the current Seurat object
  seurat_obj <- seurat_objects[[dataset_name]]
  
  # Step 1: Ensure ScType confidence scores are in the metadata
  cat("Checking and adding ScType confidence scores...\n")
  if (!("sctype_confidence" %in% colnames(seurat_obj@meta.data))) {
    cat("ScType confidence score not found in metadata. Calculating confidence scores...\n")
    
    # Check if the Seurat object uses the Assay5 structure (indicative of Seurat v5)
    if (inherits(seurat_obj[["RNA"]], "Assay5")) {
      # For Seurat v5, extract scaled data using the new method
      scRNAseqData_scaled <- as.matrix(GetAssayData(seurat_obj, slot = "scale.data", assay = "RNA"))
    } else {
      # For Seurat v4, extract scaled data using the old method
      scRNAseqData_scaled <- as.matrix(seurat_obj[["RNA"]]@scale.data)
    }
    
    # Calculate ScType scores using the prepared gene sets
    es.max <- sctype_score(scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    
    # Calculate cluster-wise results for ScType
    cL_results <- do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl) {
      es.max.cl = sort(rowSums(es.max[, rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters == cl)), 10)
    }))
    
    # Get the top score for each cluster
    sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    
    # Add a new column for ScType confidence score
    seurat_obj@meta.data$sctype_confidence <- NA
    
    # Assign confidence scores for each cluster
    for (cluster in unique(sctype_scores$cluster)) {
      confidence_score <- sctype_scores %>% filter(cluster == !!cluster) %>% pull(scores)
      
      # Apply the confidence score to all cells in the cluster
      seurat_obj@meta.data$sctype_confidence[seurat_obj@meta.data$seurat_clusters == cluster] <- confidence_score
    }
    
    # Update the Seurat object in the list
    seurat_objects[[dataset_name]] <- seurat_obj
    
    cat("Confidence scores added to metadata for dataset:", dataset_name, "\n")
  } else {
    cat("ScType confidence score already exists in metadata for dataset:", dataset_name, "\n")
  }
  
  # Step 2: Extract gene expression data and save to CSV
  
  # Check if the genes are present in the dataset
  available_genes <- genes_of_interest[genes_of_interest %in% rownames(seurat_obj)]
  if (length(available_genes) == 0) {
    cat("No genes of interest found in the dataset:", dataset_name, "\n")
    next
  }
  
  # Extract the gene expression data for the specified genes
  expression_data <- GetAssayData(seurat_obj, slot = "data")[available_genes, , drop = FALSE]
  cat("Dimensions of expression data (genes, cells):", dim(expression_data), "\n")
  
  # Ensure we only take the cells that are in the expression data
  matching_barcodes <- colnames(expression_data)
  cat("Number of matching barcodes:", length(matching_barcodes), "\n")
  
  # Subset the metadata to match the cells in expression_data
  metadata <- seurat_obj@meta.data[matching_barcodes, , drop = FALSE]
  cat("Dimensions of metadata after subsetting (cells, columns):", dim(metadata), "\n")
  
  # Check if the required columns exist in the metadata and handle missing columns
  required_columns <- c("sctype_classification", "mouse_main", "mouse_fine", "sctype_confidence")
  for (col in required_columns) {
    if (!(col %in% colnames(metadata))) {
      cat("Warning: Column", col, "is missing from metadata for dataset:", dataset_name, "\n")
      metadata[[col]] <- NA  # Add an NA-filled column if it doesn't exist
    }
  }
  
  # Combine ScType, MouseRNAseqData annotations, gene expressions, and other metadata
  combined_annotations <- data.frame(
    Barcode = matching_barcodes,
    ScType_annotation = metadata$sctype_classification,   
    MouseRNAseq_main_annotation = metadata$mouse_main,   
    MouseRNAseq_fine_annotation = metadata$mouse_fine,  
    sctype_annotation_confidence_score = metadata$sctype_confidence,  
    t(as.matrix(expression_data)),   
    Dataset = metadata$orig.ident   
  )
  
  cat("Combined data dimensions (rows, columns):", dim(combined_annotations), "\n")
  
  # Define the output folder for saving the data
  spreadsheets_folder <- file.path(base_folder, "Spreadsheets")
  if (!dir.exists(spreadsheets_folder)) {
    dir.create(spreadsheets_folder, recursive = TRUE)
  }
  
  # Define the output path for saving the data
  output_file <- file.path(spreadsheets_folder, paste0(dataset_name, "_cell_annotations_and_expression.csv"))
  
  # Save the combined annotations and expression data to CSV
  write.csv(combined_annotations, file = output_file, row.names = FALSE)
  
  cat("Saved annotations and expression data to", output_file, "\n")
}

cat("All datasets processed and CSV files saved.\n")


### Save markers in a CSV in a new folder named "Spreadsheets" ### 
# Loop through each Seurat object in the list
for (dataset_name in names(seurat_objects)) {
  cat("Processing dataset:", dataset_name, "\n")
  
  # Get the current Seurat object
  seurat_obj <- seurat_objects[[dataset_name]]
  
  # Set the default assay to "RNA"
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Find all markers (differentially expressed genes) for each cluster
  all_markers <- FindAllMarkers(seurat_obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Initialize a list to store the results for each cluster
  cluster_results <- list()
  
  # Loop through each cluster to identify top 30 overexpressed and top 10 underexpressed genes
  for (cluster in unique(all_markers$cluster)) {
    cat("Processing cluster:", cluster, "\n")
    
    cluster_markers <- all_markers[all_markers$cluster == cluster, ]
    
    # Get top 30 overexpressed genes (sorted by avg_log2FC in descending order)
    top_overexpressed <- cluster_markers %>% 
      filter(avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      head(30)
    
    # Get top 10 underexpressed genes (sorted by avg_log2FC in ascending order)
    top_underexpressed <- cluster_markers %>% 
      filter(avg_log2FC < 0) %>% 
      arrange(avg_log2FC) %>% 
      head(10)
    
    # Combine results
    combined_results <- rbind(top_overexpressed, top_underexpressed)
    combined_results$expression_type <- ifelse(combined_results$avg_log2FC > 0, "Overexpressed", "Underexpressed")
    
    # Add ScType annotation and MouseRNAseq annotations to the results
    combined_results$scType_annotation <- seurat_obj@meta.data$sctype_classification[seurat_obj@meta.data$seurat_clusters == cluster][1]
    combined_results$MouseRNAseq_main_annotation <- seurat_obj@meta.data$mouse_main[seurat_obj@meta.data$seurat_clusters == cluster][1]
    combined_results$MouseRNAseq_fine_annotation <- seurat_obj@meta.data$mouse_fine[seurat_obj@meta.data$seurat_clusters == cluster][1]
    
    # Extract the confidence score for the current cluster
    confidence_score <- sctype_scores %>% filter(cluster == !!cluster) %>% pull(scores) %>% unique()
    
    # Handle case where confidence_score is empty
    if (length(confidence_score) == 0) {
      confidence_score <- NA  # or some default value
    }
    
    # Ensure the confidence score is repeated for all rows
    combined_results$confidence_score <- rep(confidence_score, nrow(combined_results))
    
    # Add the results to the list
    cluster_results[[cluster]] <- combined_results
  }
  
  # Combine all clusters' results into a single data frame
  final_results <- do.call(rbind, cluster_results)
  
  # Save the results to a CSV file in the "Spreadsheets" folder
  output_path <- file.path(spreadsheets_folder, paste0(dataset_name, "_top_genes_per_cluster.csv"))
  write.csv(final_results, file = output_path, row.names = FALSE)
  
  cat("Finished processing dataset:", dataset_name, "Results saved to:", output_path, "\n")
}


### Save Statistically Significant Overexpressed Genes in a CSV ### 
# Define adjusted p-value threshold for statistical significance
p_val_threshold <- 0.05

# Loop through each Seurat object in the list
for (dataset_name in names(seurat_objects)) {
  cat("Processing significant overexpressed genes for dataset:", dataset_name, "\n")
  
  # Get the current Seurat object
  seurat_obj <- seurat_objects[[dataset_name]]
  
  # Set the default assay to "RNA"
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Find all markers (differentially expressed genes) for each cluster
  all_markers <- FindAllMarkers(seurat_obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Initialize a list to store the results for each cluster
  significant_results <- list()
  
  # Loop through each cluster to identify statistically significant overexpressed genes
  for (cluster in unique(all_markers$cluster)) {
    cat("Processing significant genes for cluster:", cluster, "\n")
    
    cluster_markers <- all_markers[all_markers$cluster == cluster, ]
    
    # Get statistically significant overexpressed genes (sorted by avg_log2FC in descending order)
    significant_overexpressed <- cluster_markers %>% 
      filter(avg_log2FC > 0, p_val_adj < p_val_threshold) %>% 
      arrange(desc(avg_log2FC))
    
    if (nrow(significant_overexpressed) == 0) {
      next
    }
    
    # Add ScType annotation and MouseRNAseq annotations to the results
    significant_overexpressed$scType_annotation <- seurat_obj@meta.data$sctype_classification[seurat_obj@meta.data$seurat_clusters == cluster][1]
    significant_overexpressed$MouseRNAseq_main_annotation <- seurat_obj@meta.data$mouse_main[seurat_obj@meta.data$seurat_clusters == cluster][1]
    significant_overexpressed$MouseRNAseq_fine_annotation <- seurat_obj@meta.data$mouse_fine[seurat_obj@meta.data$seurat_clusters == cluster][1]
    
    # Extract the confidence score for the current cluster
    confidence_score <- sctype_scores %>% filter(cluster == !!cluster) %>% pull(scores) %>% unique()
    
    # Handle case where confidence_score is empty
    if (length(confidence_score) == 0) {
      confidence_score <- NA  # or some default value
    }
    
    # Ensure the confidence score is repeated for all rows
    significant_overexpressed$confidence_score <- rep(confidence_score, nrow(significant_overexpressed))
    
    # Add the results to the list
    significant_results[[cluster]] <- significant_overexpressed
  }
  
  # Combine all clusters' results into a single data frame
  final_significant_results <- do.call(rbind, significant_results)
  
  # Save the results to a CSV file in the "Spreadsheets" folder
  output_path <- file.path(spreadsheets_folder, paste0(dataset_name, "_overexpressed_significant_genes.csv"))
  write.csv(final_significant_results, file = output_path, row.names = FALSE)
  
  cat("Finished processing significant genes for dataset:", dataset_name, "Results saved to:", output_path, "\n")
}

cat("All datasets processed and CSV files saved.\n")

### Save Complete Annotations and Gene Expression Data at Cluster Level ###
# Define the folder to save the cluster-level annotations and expression data in the existing Spreadsheets folder
annotations_folder <- file.path(base_folder, "Spreadsheets")

# Initialize an empty data frame to store the combined cluster-level data
combined_cluster_annotations_data <- data.frame()

# Loop through each dataset
for (dataset_name in names(seurat_objects)) {
  cat("Processing dataset for cluster-level annotations:", dataset_name, "\n")
  
  # Extract the Seurat object
  seurat_obj <- seurat_objects[[dataset_name]]
  
  # Get the metadata with cluster information and annotations
  metadata <- seurat_obj@meta.data
  
  # Ensure all required annotation columns exist
  required_columns <- c("seurat_clusters", "sctype_classification", "mouse_main", "mouse_fine")
  for (col in required_columns) {
    if (!(col %in% colnames(metadata))) {
      metadata[[col]] <- NA  # Add NA-filled column if it doesn’t exist
    }
  }
  
  # Calculate the number of cells per cluster and the percentage
  cluster_stats <- metadata %>%
    group_by(seurat_clusters) %>%
    summarize(
      num_cells = n(),
      pct_cells = (n() / nrow(metadata)) * 100,
      ScType_annotation = dplyr::first(sctype_classification),
      MouseRNAseq_main_annotation = dplyr::first(mouse_main),
      MouseRNAseq_fine_annotation = dplyr::first(mouse_fine),
      .groups = 'drop'
    )
  
  # Calculate the average gene expression per cluster and remove unwanted columns
  avg_expression_data <- as.data.frame(AverageExpression(seurat_obj, assays = "RNA", slot = "data")$RNA)
  avg_expression_data <- avg_expression_data[ , !(colnames(avg_expression_data) %in% c("g0", "g1", "g2", "g3", "g4", "g5", 
                                                                                       "g6", "g7", "g8", "g9", "g10", "g11", "g12"))]
  
  # Combine cluster annotations, cell count, cell percentage, and average gene expression
  cluster_level_data <- cbind(
    Cluster = cluster_stats$seurat_clusters, 
    cluster_stats[, -1],  # Remove the cluster column from stats to avoid duplication
    avg_expression_data[match(cluster_stats$seurat_clusters, rownames(avg_expression_data)), ]
  )
  
  # Add dataset name for reference
  cluster_level_data$Dataset <- dataset_name
  
  # Append current dataset’s cluster-level data to the combined data frame
  combined_cluster_annotations_data <- rbind(combined_cluster_annotations_data, cluster_level_data)
  
  # Define the output path for individual dataset
  output_file <- file.path(annotations_folder, paste0(dataset_name, "_cluster_level_annotations.csv"))
  
  # Save the cluster-level annotations and expression data to CSV
  write.csv(cluster_level_data, file = output_file, row.names = FALSE)
  
  cat("Saved cluster-level annotations and expression data to", output_file, "\n")
}

# Define the output file path for the combined data in the Spreadsheets folder
combined_output_file <- file.path(annotations_folder, "combined_cluster_level_annotations.csv")

# Save the combined data frame as a CSV in the Spreadsheets folder
write.csv(combined_cluster_annotations_data, file = combined_output_file, row.names = FALSE)

cat("Combined cluster-level annotations and expression data saved to", combined_output_file, "\n")
