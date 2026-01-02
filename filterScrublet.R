# scrublet_filter_module.R
# Module for applying Scrublet doublet filtering in R

library(Matrix)

#' Load Scrublet results and filter matrix to remove doublets
#'
#' @param matrix_filtered Filtered sparse matrix from EmptyDrops
#' @param barcodes_filtered Barcodes from EmptyDrops filtering
#' @param features Features/genes
#' @param scrublet_csv_path Path to Scrublet results CSV
#' @param output_dir Directory to save final filtered matrix (optional)
#'
#' @return List containing final filtered matrix, barcodes, features, and statistics
apply_scrublet_filtering <- function(matrix_filtered,
                                     barcodes_filtered,
                                     features,
                                     scrublet_csv_path,
                                     output_dir = NULL) {
  
  cat("=== Scrublet Doublet Filtering ===\n\n")
  
  # Load Scrublet results
  cat("Loading Scrublet results...\n")
  scrublet_results <- read.csv(scrublet_csv_path, stringsAsFactors = FALSE)
  
  # Convert predicted_doublet to logical if it's character
  if (is.character(scrublet_results$predicted_doublet)) {
    scrublet_results$predicted_doublet <- as.logical(scrublet_results$predicted_doublet)
  }
  
  cat("Scrublet results loaded:\n")
  cat("  Total barcodes in results:", nrow(scrublet_results), "\n")
  cat("  Doublets detected:", sum(scrublet_results$predicted_doublet), "\n")
  cat("  Doublet rate:", 
      sprintf("%.1f%%", 100 * mean(scrublet_results$predicted_doublet)), "\n\n")
  
  # Match barcodes between EmptyDrops output and Scrublet results
  cat("Matching barcodes...\n")
  
  # Get singlet barcodes from Scrublet
  singlet_barcodes <- scrublet_results$barcode[scrublet_results$predicted_doublet == FALSE]
  
  # Find indices of singlets in the current matrix
  barcode_matches <- match(singlet_barcodes, barcodes_filtered)
  
  # Remove NAs (barcodes in Scrublet but not in current matrix)
  valid_matches <- barcode_matches[!is.na(barcode_matches)]
  
  cat("  Barcodes in EmptyDrops matrix:", length(barcodes_filtered), "\n")
  cat("  Singlets from Scrublet:", length(singlet_barcodes), "\n")
  cat("  Matched singlets to keep:", length(valid_matches), "\n")
  
  # Filter matrix
  cat("\nFiltering matrix to remove doublets...\n")
  matrix_final <- matrix_filtered[, valid_matches]
  barcodes_final <- barcodes_filtered[valid_matches, drop = FALSE]
  
  cat("Final matrix:", nrow(matrix_final), "genes x", ncol(matrix_final), "cells\n")
  
  # Calculate statistics
  total_umis <- sum(matrix_final)
  mean_umis <- mean(Matrix::colSums(matrix_final))
  median_umis <- median(Matrix::colSums(matrix_final))
  
  cat("\nFinal matrix statistics:\n")
  cat("  Total UMIs:", format(total_umis, big.mark = ","), "\n")
  cat("  Mean UMIs per cell:", round(mean_umis, 1), "\n")
  cat("  Median UMIs per cell:", median_umis, "\n")
  
  # Save if output directory provided
  if (!is.null(output_dir)) {
    cat("\nSaving final filtered matrix...\n")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    writeMM(matrix_final, file = file.path(output_dir, "matrix.mtx"))
    write.table(barcodes_final,
                file = file.path(output_dir, "barcodes.tsv"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(features,
                file = file.path(output_dir, "features.tsv"),
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    
    cat("  Saved to:", output_dir, "\n")
  }
  
  # Summary statistics
  stats <- list(
    cells_before_scrublet = ncol(matrix_filtered),
    doublets_removed = sum(scrublet_results$predicted_doublet),
    cells_after_scrublet = ncol(matrix_final),
    doublet_rate = mean(scrublet_results$predicted_doublet),
    total_umis = total_umis,
    mean_umis_per_cell = mean_umis,
    median_umis_per_cell = median_umis
  )
  
  cat("\n=== Filtering Complete ===\n")
  
  return(list(
    matrix_final = matrix_final,
    barcodes_final = barcodes_final,
    features = features,
    scrublet_results = scrublet_results,
    stats = stats
  ))
}


#' Print summary of the complete filtering pipeline
#'
#' @param emptydrops_stats Statistics from EmptyDrops filtering
#' @param scrublet_stats Statistics from Scrublet filtering
print_filtering_summary <- function(emptydrops_stats, scrublet_stats) {
  
  cat("\n", rep("=", 60), "\n", sep="")
  cat("COMPLETE FILTERING PIPELINE SUMMARY\n")
  cat(rep("=", 60), "\n\n", sep="")
  
  cat("Starting barcodes:", format(emptydrops_stats$total_barcodes_seen, big.mark = ","), "\n\n")
  
  cat("After EmptyDrops:\n")
  cat("  Cells kept:", format(emptydrops_stats$total_barcodes, big.mark = ","), "\n")
  cat("  Empty droplets removed:", 
      format(emptydrops_stats$total_barcodes_seen - emptydrops_stats$total_barcodes, big.mark = ","), "\n\n")
  
  cat("After Scrublet:\n")
  cat("  Cells kept:", format(scrublet_stats$cells_after_scrublet, big.mark = ","), "\n")
  cat("  Doublets removed:", format(scrublet_stats$doublets_removed, big.mark = ","), "\n")
  cat("  Doublet rate:", sprintf("%.1f%%", 100 * scrublet_stats$doublet_rate), "\n\n")
  
  cat("Final Clean Dataset:\n")
  cat("  Total cells:", format(scrublet_stats$cells_after_scrublet, big.mark = ","), "\n")
  cat("  Total UMIs:", format(scrublet_stats$total_umis, big.mark = ","), "\n")
  cat("  Median UMIs per cell:", scrublet_stats$median_umis_per_cell, "\n")
  
  cat("\n", rep("=", 60), "\n", sep="")
}