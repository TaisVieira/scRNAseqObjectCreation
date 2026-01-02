# emptydrops_module.R
# Module for running EmptyDrops filtering on scRNA-seq data

library(Matrix)
library(DropletUtils)
library(ggplot2)

#' Run EmptyDrops filtering on a sparse matrix
#'
#' @param matrix_path Path to the matrix.mtx file
#' @param barcodes_path Path to the barcodes.tsv file
#' @param features_path Path to the features.tsv file
#' @param lower Lower threshold for UMI counts (default: 40 for sporozoites)
#' @param fdr_threshold FDR threshold for calling cells (default: 0.01)
#' @param output_dir Directory to save filtered matrix (default: "matrix_filtered")
#' @param seed Random seed for reproducibility (default: 100)
#'
#' @return List containing filtered matrix, barcodes, features, and EmptyDrops results
run_emptydrops_filtering <- function(matrix_path, 
                                     barcodes_path, 
                                     features_path,
                                     lower = 40,
                                     fdr_threshold = 0.01,
                                     output_dir = "matrix_filtered",
                                     seed = 100) {
  
  cat("=== EmptyDrops Filtering Pipeline ===\n\n")
  
  # Read the sparse matrix files
  cat("Loading data...\n")
  matrix <- readMM(matrix_path)
  barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  features <- read.table(features_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  cat("Matrix dimensions:", nrow(matrix), "genes x", ncol(matrix), "barcodes\n\n")
  
  # Calculate UMI statistics
  umi_per_barcode <- Matrix::colSums(matrix)
  
  cat("Pre-filtering statistics:\n")
  cat("  Total UMIs:", sum(umi_per_barcode), "\n")
  cat("  Barcodes with >0 UMIs:", sum(umi_per_barcode > 0), "\n")
  cat("  Median UMIs (>0):", median(umi_per_barcode[umi_per_barcode > 0]), "\n\n")
  
  # Run EmptyDrops
  cat("Running EmptyDrops with lower =", lower, "...\n")
  set.seed(seed)
  #matrix_filtered_0 <- matrix[, colSums(matrix) > 0]
  e.out <- emptyDrops(matrix, lower = lower)
  
  cat("EmptyDrops completed!\n\n")
  
  # Summary of results
  cat("EmptyDrops Results:\n")
  cat("  Barcodes tested:", sum(!is.na(e.out$FDR)), "\n")
  cat("  Cells called (FDR <", fdr_threshold, "):", sum(e.out$FDR < fdr_threshold, na.rm = TRUE), "\n")
  
  result_table <- table(Sig = e.out$FDR < fdr_threshold, 
                        Limited = e.out$Limited, 
                        useNA = "always")
  print(result_table)
  cat("\n")
  
  # Filter matrix
  cells_to_keep <- which(e.out$FDR < fdr_threshold)
  matrix_filtered <- matrix[, cells_to_keep]
  barcodes_filtered <- barcodes[cells_to_keep, ]
  
  cat("Filtered matrix:", nrow(matrix_filtered), "genes x", ncol(matrix_filtered), "cells\n")
  cat("Total UMIs in filtered matrix:", sum(matrix_filtered), "\n")
  cat("Mean UMIs per cell:", round(mean(Matrix::colSums(matrix_filtered)), 1), "\n")
  cat("Median UMIs per cell:", median(Matrix::colSums(matrix_filtered)), "\n\n")
  
  # Save filtered matrix
  cat("Saving filtered matrix to", output_dir, "/\n")
  dir.create(output_dir, showWarnings = FALSE)
  
  writeMM(matrix_filtered, file = file.path(output_dir, "matrix.mtx"))
  write.table(barcodes_filtered, 
              file = file.path(output_dir, "barcodes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(features,
              file = file.path(output_dir, "features.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  cat("Filtered matrix saved successfully!\n\n")
  
  # Create diagnostic plot
  cat("Creating diagnostic plot...\n")
  plot_obj <- create_emptydrops_diagnostic_plot(e.out, fdr_threshold)
  
  # Collect statistics 
  stats <- list(
    total_barcodes_seen = ncol(matrix),
    total_barcodes = ncol(matrix_filtered),
    total_genes = nrow(matrix_filtered),
    total_umis = sum(matrix_filtered),
    mean_umis_per_cell = mean(Matrix::colSums(matrix_filtered)),
    median_umis_per_cell = median(Matrix::colSums(matrix_filtered))
  )
  
  return(list(
    matrix_filtered = matrix_filtered,
    barcodes_filtered = barcodes_filtered,
    features = features,
    emptydrops_results = e.out,
    cells_kept = cells_to_keep,
    diagnostic_plot = plot_obj,
    stats = stats
  ))
}

#' Create EmptyDrops diagnostic plot
#'
#' @param e.out EmptyDrops output object
#' @param fdr_threshold FDR threshold used for calling cells
#'
#' @return ggplot object
create_emptydrops_diagnostic_plot <- function(e.out, fdr_threshold = 0.01) {
  
  emptydrops_df <- data.frame(
    Total = e.out$Total,
    LogProb = e.out$LogProb,
    FDR = e.out$FDR,
    Call = ifelse(is.na(e.out$FDR), "Not tested",
                  ifelse(e.out$FDR < fdr_threshold, "Cell", "Empty"))
  )
  
  p <- ggplot(emptydrops_df[!is.na(e.out$FDR), ], 
              aes(x = Total, y = LogProb, color = Call)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_x_log10() +
    scale_color_manual(values = c("Cell" = "blue", "Empty" = "red")) +
    labs(
      title = "EmptyDrops: Log-Probability vs Total UMIs",
      x = "Total UMIs (log scale)",
      y = "Log-Probability under ambient profile",
      color = "Classification"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(p)
  ggsave("emptydrops_diagnostic.png", p, width = 10, height = 6, dpi = 300)
  
  return(p)
}