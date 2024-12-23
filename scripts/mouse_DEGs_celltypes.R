# ------------------------------------------------------------------------------------------------
# R Script for Differential Gene Expression and Gene Ontology Pathway Analysis of mouse scRNA data
# ------------------------------------------------------------------------------------------------
# 
# Description:
# This script performs the following analyses:
# 
# 1. Differential Gene Expression (DGE) Analysis:
#    - Compares differentially expressed genes across datasets and cell types.
#    - Supports the identification of biologically significant genes using customizable thresholds.
# 
# 2. Gene Ontology (GO) and Pathway Analysis:
#    - Conducts Gene Set Enrichment Analysis (GSEA) to identify enriched biological processes 
#      and pathways.
#    - Allows comparisons across the entire dataset or between specific cell types.
# 
# 3. Visualization:
#    - Generates plots and spreadsheets for genes of interest to visualize expression patterns.
# 
# ------------------------------------------------------------------------------------------------
# Requirements:
# 
# 1. Set the path to the `mouse_bothannot_functions.R` script:
#    - This script contains the helper functions required for preprocessing, analysis, and visualization.
# 
# 2. Specify the folders of the input files:
#    - Input data (Seurat objects with appropriate metadata) should be provided.
# 
# 3. Set the base folder for output:
#    - Define the directory where all output files, including UMAP plots, heatmaps, and 
#      annotation files, will be saved.
# 
# -------------------------------------------------------------------------------------------------
# Libraries and Functions: 
# - Functions and libraries are loaded from:
#   `mouse_scRNA_functions.R` for essential preprocessing steps and custom utilities.
# -------------------------------------------------------------------------------------------------

# Load required functions and libraries
source("/data/ldiaz/vlachog1/scRNA_mice/R/mouse_scRNA_functions.R")

# Set working directory
setwd("/data/ldiaz/vlachog1/scRNA_mice")

# Set base folder for output
base_folder <- file.path(getwd(), "output/DEGs_Cell_Types/test")
if (!dir.exists(base_folder)) {
  dir.create(base_folder, recursive = TRUE)
  cat("Base folder created:", base_folder, "\n")
} else {
  cat("Base folder already exists:", base_folder, "\n")
}

# Define Seurat object file paths
seurat_paths <- list(
  "Parental" = "/data/ldiaz/vlachog1/scRNA_mice/output/Seurat_objects/Parental_final_seurat.rds",
  "CDDP" = "/data/ldiaz/vlachog1/scRNA_mice/output/Seurat_objects/CDDP_final_seurat.rds",
  "TMZ" = "/data/ldiaz/vlachog1/scRNA_mice/output/Seurat_objects/TMZ_final_seurat.rds",
  "Combo" = "/data/ldiaz/vlachog1/scRNA_mice/output/Seurat_objects/Combo_final_seurat.rds"
)

# Define cell subtypes of interest
cell_subtypes_of_interest <- c(
  "Endothelial cells",
  "Cancer cells"
)


# Define comparisons
comparisons <- list(
  list(group1 = c("Parental"), group2 = c("CDDP"), comparison_name = "Parental_CDDP")
)

# Loop through each defined comparison and cell subtype
for (cell_subtype in cell_subtypes_of_interest) {
  # Reload Seurat objects for each cell subtype
  seurat_objects <- list()
  for (sample_name in names(seurat_paths)) {
    seurat_objects[[sample_name]] <- readRDS(seurat_paths[[sample_name]])
    cat("Reloaded Seurat object:", sample_name, "\n")
  }
  
  for (comparison in comparisons) {
    group1_samples <- comparison$group1
    group2_samples <- comparison$group2
    comparison_name <- paste0(comparison$comparison_name, "_", cell_subtype)
    
    cat("Starting analysis for comparison:", comparison_name, "in cell subtype:", cell_subtype, "\n")
    
    # Create output folder for comparison
    comparison_folder <- file.path(base_folder, comparison_name)
    if (!dir.exists(comparison_folder)) {
      dir.create(comparison_folder, recursive = TRUE)
    }
    cat("Directory created for comparison:", comparison_folder, "\n")
    
    # Check if samples exist in loaded Seurat objects
    group1_samples <- group1_samples[group1_samples %in% names(seurat_objects)]
    group2_samples <- group2_samples[group2_samples %in% names(seurat_objects)]
    if (length(group1_samples) == 0 || length(group2_samples) == 0) {
      cat("Skipping comparison due to missing samples in", comparison_name, "\n")
      next
    }
    
    # Subset Seurat objects by cell subtype and check if the cell subtype exists
    valid_comparison <- TRUE
    for (sample in c(group1_samples, group2_samples)) {
      subset_seurat <- subset(seurat_objects[[sample]], subset = sctype_classification == cell_subtype)
      
      if (nrow(subset_seurat@meta.data) == 0) {
        cat("No cells found for cell subtype", cell_subtype, "in sample", sample, "\n")
        valid_comparison <- FALSE
        break
      } else {
        cat("Number of", cell_subtype, "cells in", sample, ":", nrow(subset_seurat@meta.data), "\n")
        seurat_objects[[sample]] <- subset_seurat
      }
    }
    
    if (!valid_comparison) {
      cat("Skipping comparison due to missing cell subtype", cell_subtype, "\n")
      next
    }
    
    # Proceed with normalization, integration, and analysis
    for (sample in c(group1_samples, group2_samples)) {
      DefaultAssay(seurat_objects[[sample]]) <- "RNA"
      seurat_objects[[sample]] <- NormalizeData(seurat_objects[[sample]])
    }
    
    # Select integration features and find anchors
    cat("Selecting integration features for", comparison_name, "\n")
    features <- SelectIntegrationFeatures(object.list = lapply(c(group1_samples, group2_samples), function(sample) seurat_objects[[sample]]))
    anchors <- FindIntegrationAnchors(object.list = lapply(c(group1_samples, group2_samples), function(sample) seurat_objects[[sample]]), anchor.features = features)
    
    # Integrate data
    cat("Integrating data for comparison:", comparison_name, "\n")
    integrated_data <- IntegrateData(anchorset = anchors, k.weight = 30)
    
    # Ensure data layers are joined in the RNA assay
    integrated_data <- JoinLayers(integrated_data, assay = "RNA")
    
    # Assign unique cell names and dataset identity
    integrated_data <- RenameCells(integrated_data, add.cell.id = comparison_name)
    integrated_data$dataset <- factor(integrated_data$orig.ident, levels = c(group1_samples, group2_samples))
    Idents(integrated_data) <- "dataset"
    
    # Debugging step: Check the unique identities in the integrated data
    cat("Unique identities in integrated data:\n")
    print(unique(integrated_data@meta.data$orig.ident))
    
    # Perform scaling, PCA, and UMAP
    DefaultAssay(integrated_data) <- "RNA"
    integrated_data <- ScaleData(integrated_data)
    integrated_data <- RunPCA(integrated_data)
    integrated_data <- RunUMAP(integrated_data, dims = 1:30)
    
    # Differential Expression Analysis
    de_results <- FindMarkers(integrated_data, ident.1 = group1_samples[1], ident.2 = group2_samples[1], assay = "RNA", logfc.threshold = 0.25, min.pct = 0.1)
    de_results <- de_results %>%
      rownames_to_column(var = "Gene") %>%
      mutate(logFC = avg_log2FC, p_value = p_val, p_adjusted = p_val_adj) %>%
      dplyr::select(Gene, logFC, p_value, p_adjusted, pct.1, pct.2)
    
    # Save DE results
    de_output_file <- file.path(comparison_folder, paste0("DE_results_RNA_", comparison_name, ".csv"))
    write.csv(de_results, file = de_output_file, row.names = FALSE)
    
    # Define cell counts for each group
    cell_counts <- list()
    for (sample in c(group1_samples, group2_samples)) {
      # Subset Seurat object for the current cell subtype
      subset_seurat <- subset(seurat_objects[[sample]], subset = sctype_classification == cell_subtype)
      
      # Check if any cells were found, count cells, and update seurat object if valid
      if (nrow(subset_seurat@meta.data) > 0) {
        cell_counts[[sample]] <- nrow(subset_seurat@meta.data)
        seurat_objects[[sample]] <- subset_seurat
        cat("Number of", cell_subtype, "cells in", sample, ":", cell_counts[[sample]], "\n")
      } else {
        cat("No cells found for cell subtype", cell_subtype, "in sample", sample, "\n")
        cell_counts[[sample]] <- 0
      }
    }
    
    # Update volcano plot labels to include a line break before the cell counts
    group1_label <- paste("Upregulated in", group1_samples[1], "\n(", cell_counts[[group1_samples[1]]], "cells)")
    group2_label <- paste("Upregulated in", group2_samples[1], "\n(", cell_counts[[group2_samples[1]]], "cells)")
    
    # Define columns for volcano plot with dynamic group labels
    de_results <- de_results %>%
      mutate(
        expression = case_when(
          logFC < 0 ~ group2_label,
          logFC > 0 ~ group1_label,
          TRUE ~ "NS"
        )
      )
    
    # Color mapping for samples
    color_mapping <- setNames(c("blue", "red", "grey"), c(group1_label, group2_label, "NS"))
    
    # Create volcano plot with updated title and legend labels
    volcano_plot <- ggplot(de_results, aes(x = logFC, y = -log10(p_value))) +
      geom_point(aes(color = expression), size = 2) +
      scale_color_manual(values = color_mapping) +
      theme_minimal() +
      labs(
        title = paste("Volcano Plot of DEGs:", comparison_name),
        x = "Log2 Fold Change",
        y = "-Log10 P-Value",
        color = paste(cell_subtype)
      ) +
      geom_text_repel(
        data = de_results %>%
          filter(p_adjusted < 0.05) %>%
          arrange(p_value, desc(abs(logFC))) %>%
          slice(1:20),
        aes(label = Gene),
        max.overlaps = 50,
        box.padding = 0.5,
        point.padding = 0.3,
        size = 3.5
      )
    
    # Save volcano plot
    volcano_output_path <- file.path(comparison_folder, paste0("Volcano_Plot_DEGs_", comparison_name, ".png"))
    ggsave(volcano_output_path, plot = volcano_plot, dpi = 600, width = 8, height = 12)
    print(volcano_plot)
    
    
    
    # GSEA for GO and KEGG pathways
    cat("Starting GSEA analysis for RNA DEGs\n")
    output_dir <- file.path(comparison_folder, "GSEA_results")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Map gene identifiers from DE results
    mapped_genes <- bitr(de_results$Gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
    de_results_mapped <- merge(de_results, mapped_genes, by.x = "Gene", by.y = "SYMBOL")
    
    ranked_gene_list_go <- de_results_mapped$logFC
    names(ranked_gene_list_go) <- de_results_mapped$ENSEMBL
    ranked_gene_list_go <- sort(ranked_gene_list_go, decreasing = TRUE)
    
    # GSEA for GO Terms
    gsea_go <- gseGO(
      geneList = ranked_gene_list_go,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      keyType = "ENSEMBL",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05
    )
    go_output_csv <- file.path(output_dir, paste0("gsea_go_results_", comparison_name, ".csv"))
    write.csv(as.data.frame(gsea_go), file = go_output_csv)
    cat("GO GSEA results saved to:", go_output_csv, "\n")
    
    # Dot Plot of Top GO Terms if available
    if (!is.null(gsea_go) && length(gsea_go$ID) > 0) {
      dotplot_go <- dotplot(gsea_go, showCategory = 10) + 
        ggtitle(paste("Top 10 GO Biological Processes -", comparison_name))
      ggsave(filename = file.path(output_dir, paste0("go_dotplot_", comparison_name, ".png")), plot = dotplot_go, width = 10, height = 7)
    }
  }
}


