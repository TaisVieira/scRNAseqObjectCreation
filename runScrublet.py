import scrublet as scr
import pandas as pd
import scipy.io
import scipy.sparse
import os
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def load_data(data_dir):
    """
    Load 10X-style data (matrix, barcodes, features).
    
    Parameters:
    -----------
    data_dir : str
        Path to directory containing matrix.mtx, barcodes.tsv, features.tsv
    
    Returns:
    --------
    tuple
        (matrix, barcodes, features) - matrix is cells x genes for Scrublet
    """
    data_dir = Path(data_dir)
    print("Loading data...")
    
    # Load matrix (transpose to get cells x genes for Scrublet)
    matrix = scipy.io.mmread(data_dir / "matrix.mtx").T
    matrix = scipy.sparse.csr_matrix(matrix)
    
    # Load barcodes and features
    barcodes = pd.read_csv(data_dir / "barcodes.tsv", header=None)[0].values
    features = pd.read_csv(data_dir / "features.tsv", header=None)[0].values
    
    print(f"Loaded matrix: {matrix.shape} (cells x genes)")
    print(f"Barcodes: {len(barcodes)}")
    print(f"Features: {len(features)}")
    
    return matrix, barcodes, features


def run_scrublet_analysis(matrix, expected_doublet_rate=0.075):
    """
    Run Scrublet doublet detection.
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Expression matrix (cells x genes)
    expected_doublet_rate : float
        Expected doublet rate for the experiment
    
    Returns:
    --------
    tuple
        (scrub, doublet_scores, predicted_doublets)
    """
    print("Running Scrublet doublet detection...")
    
    scrub = scr.Scrublet(matrix, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
    print(f"Doublet detection complete. Auto threshold: {scrub.threshold_:.4f}")
    return scrub, doublet_scores, predicted_doublets


def plot_score_distributions(scrub, percentiles=[80, 95, 99], save_path=None):
    """
    Plot observed and simulated doublet linear score distributions with key thresholds.

    Parameters:
    -----------
    scrub : Scrublet object
        Fitted Scrublet object
    percentiles : list of int, optional
        List of percentiles to plot (default: [80, 95, 99])
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    dict
        Dictionary with percentile values, e.g. {80: p80, 95: p95, 99: p99}
    """
    # Compute requested percentiles from simulated distribution
    perc_values = {p: np.percentile(scrub.doublet_scores_sim_, p) for p in percentiles}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Colors for percentiles
    colors = ['orange', 'blue', 'green', 'purple', 'brown', 'pink']  # extend as needed

    # Observed doublet scores
    ax1.hist(scrub.doublet_scores_obs_, bins=50, alpha=0.7, color='lightblue', edgecolor='black')
    ax1.axvline(scrub.threshold_, color='red', linestyle='--', linewidth=2, label=f'Auto threshold: {scrub.threshold_:.3f}')
    for i, p in enumerate(percentiles):
        ax1.axvline(perc_values[p], color=colors[i % len(colors)], linestyle='--', linewidth=2,
                    label=f'{p}th percentile: {perc_values[p]:.3f}')
    ax1.set_xlabel('Doublet Score')
    ax1.set_ylabel('Count')
    ax1.set_title('Observed Transcriptomes (linear count)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Simulated doublet scores
    ax2.hist(scrub.doublet_scores_sim_, bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
    ax2.axvline(scrub.threshold_, color='red', linestyle='--', linewidth=2, label=f'Auto threshold: {scrub.threshold_:.3f}')
    for i, p in enumerate(percentiles):
        ax2.axvline(perc_values[p], color=colors[i % len(colors)], linestyle='--', linewidth=2,
                    label=f'{p}th percentile: {perc_values[p]:.3f}')
    ax2.set_xlabel('Doublet Score')
    ax2.set_ylabel('Count')
    ax2.set_title('Simulated Doublets (linear count)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()
    return perc_values


def compare_percentiles(scrub, percentiles=[70, 80, 85, 90, 95]):
    """Compare multiple percentile thresholds quickly."""
    print("\nPercentile threshold comparison:")
    print("-" * 50)
    
    for p in percentiles:
        threshold = np.percentile(scrub.doublet_scores_sim_, p)
        pred_mask = scrub.doublet_scores_obs_ >= threshold
        rate = pred_mask.mean()
        n_doublets = pred_mask.sum()
        
        print(f"p{p:2d}: {threshold:.4f} -> {n_doublets:4,} doublets ({rate:.1%})")


def analyze_umi_correlation(matrix, scrub, threshold_percentile=80):
    """
    Analyze correlation between doublet scores and UMI counts.
    
    Parameters:
    -----------
    matrix : scipy.sparse matrix
        Expression matrix (cells x genes)
    scrub : Scrublet object
        Fitted Scrublet object
    threshold_percentile : float
        Percentile of simulated scores to use as threshold
    
    Returns:
    --------
    dict
        Dictionary with correlation statistics
    """
    # Calculate UMI counts per cell
    n_counts = np.array(matrix.sum(axis=1)).flatten()
    
    # Apply threshold
    threshold = float(np.percentile(scrub.doublet_scores_sim_, threshold_percentile))
    doublet_mask = scrub.doublet_scores_obs_ >= threshold
    
    print(f"\nUMI Count Analysis (using {threshold_percentile}th percentile threshold):")
    print(f"Threshold: {threshold:.4f}")
    print(f"Predicted doublet rate: {doublet_mask.mean():.1%}")
    print(f"Median UMI counts - Singlets: {np.median(n_counts[~doublet_mask]):,.0f}")
    print(f"Median UMI counts - Doublets: {np.median(n_counts[doublet_mask]):,.0f}")
    
    # Calculate correlation
    correlation = np.corrcoef(scrub.doublet_scores_obs_, n_counts)[0, 1]
    print(f"Correlation between doublet scores and UMI counts: {correlation:.3f}")
    
    return {
        'threshold': threshold,
        'doublet_rate': doublet_mask.mean(),
        'median_umi_singlets': np.median(n_counts[~doublet_mask]),
        'median_umi_doublets': np.median(n_counts[doublet_mask]),
        'correlation': correlation
    }


def set_manual_threshold(scrub, threshold):
    """
    Manually set doublet detection threshold.
    
    Parameters:
    -----------
    scrub : Scrublet object
        Scrublet object to update
    threshold : float
        Doublet score threshold to use
    """
    scrub.call_doublets(threshold=float(threshold))
    print(f"Manual threshold set: {threshold:.4f}")
    print(f"Predicted doublets: {np.sum(scrub.predicted_doublets_):,} "
          f"({scrub.predicted_doublets_.mean():.1%})")


def save_results(barcodes, doublet_scores, predicted_doublets, output_path):
    """
    Save doublet detection results to CSV.
    
    Parameters:
    -----------
    barcodes : np.ndarray
        Cell barcodes
    doublet_scores : np.ndarray
        Doublet scores
    predicted_doublets : np.ndarray
        Predicted doublet labels
    output_path : str
        Path to save the CSV file
    """
    # Ensure all arrays have the same length
    n_cells = len(doublet_scores)
    
    results_df = pd.DataFrame({
        'barcode': barcodes[:n_cells],
        'doublet_score': doublet_scores,
        'predicted_doublet': predicted_doublets
    })
    
    results_df.to_csv(output_path, index=False)
    print(f"Results saved to: {output_path}")

    
#################### run analysis ####################

data_dir = "./matrix/scrublet_input"
output_dir = "./matrix"
expected_doublet_rate = 0.075

# 1. Load data and run initial Scrublet analysis
matrix, barcodes, features = load_data(data_dir)
scrub, doublet_scores, predicted_doublets = run_scrublet_analysis(matrix, expected_doublet_rate)

# 2. Initial visualization and threshold exploration

# Plot score distributions
scrub.plot_histogram()
thresholds = plot_score_distributions(scrub, percentiles=[90, 95, 99], save_path=Path(output_dir) / "doublet_score_distributions.png")

# Show threshold comparison
threshold_df = compare_percentiles(scrub, percentiles=[90, 95, 99])

# Analyze UMI correlation
umi_stats = analyze_umi_correlation(matrix, scrub, threshold_percentile=95)

# 3. Set your final threshold 
final_threshold = np.percentile(scrub.doublet_scores_sim_, 95)
set_manual_threshold(scrub, final_threshold)

# 4. Save results
save_results(barcodes, doublet_scores, scrub.predicted_doublets_, 
              Path(output_dir) / "scrublet_doublet_scores.csv")

