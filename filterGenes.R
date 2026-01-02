# quality_filter_module.R
# Additional quality filtering: genes and cells

library(Matrix)
library(Seurat)

#' Apply additional quality filters to the matrix
#'
#' @param matrix_final Final filtered matrix from Scrublet
#' @param barcodes_final Final barcodes
#' @param features Final features/genes
#' @param min_genes Minimum number of genes per cell (default: 60)
#' @param min_cells Minimum number of cells per gene (default: 2)
#' @param min_umis_per_gene Minimum UMIs per gene in those cells (default: 2)
#' @param project_name Project name for Seurat object (default: "scRNAseq")
#' @param output_dir Directory to save final matrix (optional)
#'
#' @return List containing Seurat object and statistics
apply_quality_filters <- function(matrix_final,
                                  barcodes_final,
                                  features,
                                  min_genes = 60,
                                  min_cells = 2,
                                  min_umis_per_gene = 2,
                                  project_name = "scRNAseq",
                                  output_dir = NULL) {
  
  cat("=== Quality Filtering ===\n\n")
  
  # Create initial Seurat object with all data
  cat("Creating initial Seurat object...\n")
  
  # Prepare data for Seurat
  rownames(matrix_final) <- features[, 1]
  if (is.data.frame(barcodes_final)) {
    colnames(matrix_final) <- barcodes_final[, 1]
  } else {
    colnames(matrix_final) <- barcodes_final
  }
  
  seurat_obj <- CreateSeuratObject(
    counts = matrix_final,
    project = project_name
  )
  
  cat("Initial Seurat object created:\n")
  cat("  Cells:", ncol(seurat_obj), "\n")
  cat("  Genes:", nrow(seurat_obj), "\n\n")
  
  # Calculate metrics before filtering
  n_genes_per_cell <- Matrix::colSums(matrix_final > 0)
  n_cells_per_gene <- Matrix::rowSums(matrix_final > 0)
  max_umi_per_gene <- apply(matrix_final, 1, max)
  
  cat("Pre-filtering statistics:\n")
  cat("  Median genes per cell:", median(n_genes_per_cell), "\n")
  cat("  Mean genes per cell:", round(mean(n_genes_per_cell), 1), "\n")
  cat("  Cells with <", min_genes, "genes:", sum(n_genes_per_cell < min_genes), "\n\n")
  
  cat("  Median cells per gene:", median(n_cells_per_gene), "\n")
  cat("  Genes in <", min_cells, "cells:", sum(n_cells_per_gene < min_cells), "\n\n")
  
  # Filter 1: Remove cells with < min_genes
  cat("Applying cell filter (min genes =", min_genes, ")...\n")
  cells_to_keep <- n_genes_per_cell >= min_genes
  
  cat("  Cells removed:", sum(!cells_to_keep), "\n")
  cat("  Cells remaining:", sum(cells_to_keep), "\n\n")
  
  # Subset Seurat object
  seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
  
  # Filter 2: Remove genes detected in < min_cells with < min_umis_per_gene UMIs
  cat("Applying gene filter (min cells =", min_cells, 
      ", min UMIs per gene =", min_umis_per_gene, ")...\n")
  
  # Get the updated matrix after cell filtering
  current_matrix <- GetAssayData(seurat_obj, layer = "counts")
  
  # For each gene, count cells with >= min_umis_per_gene
  cells_with_sufficient_umis <- Matrix::rowSums(current_matrix >= min_umis_per_gene)
  genes_to_keep <- cells_with_sufficient_umis >= min_cells
  
  cat("  Genes removed:", sum(!genes_to_keep), "\n")
  cat("  Genes remaining:", sum(genes_to_keep), "\n\n")
  
  # Subset Seurat object by genes
  seurat_obj <- subset(seurat_obj, features = rownames(seurat_obj)[genes_to_keep])
  
  # Final statistics
  final_matrix <- GetAssayData(seurat_obj, layer = "counts")
  
  cat("Final filtered Seurat object:\n")
  cat("  Cells:", ncol(seurat_obj), "\n")
  cat("  Genes:", nrow(seurat_obj), "\n")
  cat("  Total UMIs:", sum(final_matrix), "\n")
  cat("  Mean UMIs per cell:", round(mean(Matrix::colSums(final_matrix)), 1), "\n")
  cat("  Median UMIs per cell:", median(Matrix::colSums(final_matrix)), "\n")
  cat("  Mean genes per cell:", round(mean(Matrix::colSums(final_matrix > 0)), 1), "\n")
  cat("  Median genes per cell:", median(Matrix::colSums(final_matrix > 0)), "\n\n")
  
  # Save if requested
  if (!is.null(output_dir)) {
    cat("Saving final filtered matrix...\n")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    writeMM(final_matrix, file = file.path(output_dir, "matrix.mtx"))
    write.table(colnames(final_matrix),
                file = file.path(output_dir, "barcodes.tsv"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(data.frame(gene = rownames(final_matrix)),
                file = file.path(output_dir, "features.tsv"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cat("  Saved to:", output_dir, "\n\n")
  }
  
  # Compile statistics
  stats <- list(
    cells_before = ncol(matrix_final),
    cells_removed = sum(!cells_to_keep),
    cells_after = ncol(seurat_obj),
    genes_before = nrow(matrix_final),
    genes_removed = sum(!genes_to_keep),
    genes_after = nrow(seurat_obj),
    total_umis = sum(final_matrix),
    mean_umis_per_cell = mean(Matrix::colSums(final_matrix)),
    median_umis_per_cell = median(Matrix::colSums(final_matrix)),
    mean_genes_per_cell = mean(Matrix::colSums(final_matrix > 0)),
    median_genes_per_cell = median(Matrix::colSums(final_matrix > 0))
  )
  
  cat("=== Quality Filtering Complete ===\n")
  
  return(list(
    seurat_obj = seurat_obj,
    stats = stats
  ))
}


#' Print complete filtering pipeline summary
#'
#' @param emptydrops_stats Statistics from EmptyDrops
#' @param rrna_stats Statistics from rRNA filtering
#' @param scrublet_stats Statistics from Scrublet
#' @param quality_stats Statistics from quality filtering
print_complete_filtering_summary <- function(emptydrops_stats, 
                                             rrna_stats,
                                             scrublet_stats, 
                                             quality_stats) {
  
  cat("\n", rep("=", 50), "\n", sep="")
  cat("COMPLETE FILTERING PIPELINE SUMMARY\n")
  cat(rep("=", 50), "\n\n", sep="")
  
  cat("Step 1: EmptyDrops - Remove empty droplets\n")
  cat("  Starting barcodes: ", format(emptydrops_stats$total_barcodes_seen, big.mark = ","), "\n")
  cat("  Cells kept: ", format(emptydrops_stats$total_barcodes, big.mark = ","), "\n")
  cat("  Empty droplets removed:", 
      format(emptydrops_stats$total_barcodes_seen - emptydrops_stats$total_barcodes, big.mark = ","), 
      sprintf(" (%.1f%%)\n\n", 
              100 * (emptydrops_stats$total_barcodes_seen - emptydrops_stats$total_barcodes) / emptydrops_stats$total_barcodes_seen))
  
  cat("Step 2: rRNA Filtering - Remove ribosomal RNA contamination\n")
  cat("  Genes before: ", format(rrna_stats$genes_before, big.mark = ","), "\n")
  cat("  Genes kept: ", format(rrna_stats$genes_after, big.mark = ","), "\n")
  cat("  rRNA genes removed: ", format(rrna_stats$rrna_genes_removed, big.mark = ","))
  
  if (rrna_stats$rrna_genes_removed > 0) {
    cat(" (", paste(rrna_stats$rrna_gene_names, collapse = ", "), ")", sep="")
  }
  cat("\n\n")
  
  cat("Step 3: Scrublet - Remove doublets\n")
  cat("  Cells before: ", format(scrublet_stats$cells_before_scrublet, big.mark = ","), "\n")
  cat("  Cells kept: ", format(scrublet_stats$cells_after_scrublet, big.mark = ","), "\n")
  cat("  Doublets removed: ", format(scrublet_stats$doublets_removed, big.mark = ","),
      sprintf(" (%.1f%%)\n\n", 100 * scrublet_stats$doublet_rate))
  
  cat("Step 4: Quality Filtering - Remove low-quality cells and genes\n")
  cat("  Cells before: ", format(quality_stats$cells_before, big.mark = ","), "\n")
  cat("  Cells kept: ", format(quality_stats$cells_after, big.mark = ","), "\n")
  cat("  Low-quality cells: ", format(quality_stats$cells_removed, big.mark = ","),
      sprintf(" (%.1f%%)\n", 100 * quality_stats$cells_removed / quality_stats$cells_before))
  
  cat("  Genes before: ", format(quality_stats$genes_before, big.mark = ","), "\n")
  cat("  Genes kept: ", format(quality_stats$genes_after, big.mark = ","), "\n")
  cat("  Low-detection genes: ", format(quality_stats$genes_removed, big.mark = ","),
      sprintf(" (%.1f%%)\n\n", 100 * quality_stats$genes_removed / quality_stats$genes_before))
  
  cat("FINAL CLEAN DATASET\n")
  cat(rep("-", 70), "\n", sep="")
  cat("  Total cells: ", format(quality_stats$cells_after, big.mark = ","), "\n")
  cat("  Total genes: ", format(quality_stats$genes_after, big.mark = ","), "\n")
  cat("  Total UMIs: ", format(quality_stats$total_umis, big.mark = ","), "\n")
  cat("  Median UMIs per cell: ", quality_stats$median_umis_per_cell, "\n")
  cat("  Median genes per cell: ", quality_stats$median_genes_per_cell, "\n")
  
  # Calculate overall filtering stats
  total_barcodes_removed <- emptydrops_stats$total_barcodes_seen - quality_stats$cells_after
  total_genes_removed <- emptydrops_stats$total_genes - quality_stats$genes_after
  
  cat("\n  Overall filtering:\n")
  cat("    Barcodes removed: ", format(total_barcodes_removed, big.mark = ","),
      sprintf(" (%.1f%%)\n", 100 * total_barcodes_removed / emptydrops_stats$total_barcodes_seen))
  cat("    Genes removed: ", format(total_genes_removed, big.mark = ","),
      sprintf(" (%.1f%%)\n", 100 * total_genes_removed / emptydrops_stats$total_genes))
  
  cat("\n", rep("=", 50), "\n", sep="")
}