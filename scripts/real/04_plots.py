#!/usr/bin/env python3
"""
Generate diagnostic plots for QTL predictions.

Creates:
1. Scatter plot: predicted vs observed effect sizes
2. ROC curve: classification performance
3. Precision-Recall curve
4. Calibration plot
5. Distribution comparisons

Input:
- results/predictions/{dataset}/ag_predictions.parquet
- results/tables/{dataset}_metrics.tsv

Output:
- results/plots/{dataset}/{modality}/*.png
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import roc_curve, precision_recall_curve, auc

sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300


def plot_scatter(df: pd.DataFrame, output_dir: Path, dataset: str, modality: str):
    """Scatter plot of predicted vs observed effects."""
    # Aggregate by variant
    if 'quantile_score' in df.columns:
        df['abs_score'] = df['quantile_score'].abs()
        idx = df.groupby('variant_id')['abs_score'].idxmax()
        agg_df = df.loc[idx].copy()
        score_col = 'quantile_score'
    elif 'raw_score' in df.columns:
        df['abs_score'] = df['raw_score'].abs()
        idx = df.groupby('variant_id')['abs_score'].idxmax()
        agg_df = df.loc[idx].copy()
        score_col = 'raw_score'
    else:
        agg_df = df.copy()
        score_col = 'delta'
    
    valid = ~(agg_df[score_col].isna() | agg_df['beta'].isna())
    scores = agg_df.loc[valid, score_col].values
    beta = agg_df.loc[valid, 'beta'].values
    
    if len(scores) < 2:
        return
    
    r, p = spearmanr(scores, beta)
    
    plt.figure(figsize=(8, 8))
    plt.scatter(beta, scores, alpha=0.3, s=20)
    plt.xlabel('Observed Effect Size (β)', fontsize=12)
    plt.ylabel(f'Predicted Score ({score_col})', fontsize=12)
    plt.title(f'{dataset} - {modality}\nSpearman r = {r:.3f}, p = {p:.3e}', fontsize=14)
    plt.axhline(0, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Add regression line
    z = np.polyfit(beta, scores, 1)
    p_reg = np.poly1d(z)
    x_line = np.linspace(beta.min(), beta.max(), 100)
    plt.plot(x_line, p_reg(x_line), "r--", alpha=0.8, linewidth=2)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'{dataset}_{modality}_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_roc(df: pd.DataFrame, output_dir: Path, dataset: str, modality: str):
    """ROC curve for binary classification."""
    # Aggregate by variant
    if 'quantile_score' in df.columns:
        df['abs_score'] = df['quantile_score'].abs()
        idx = df.groupby('variant_id')['abs_score'].idxmax()
        agg_df = df.loc[idx].copy()
        score_col = 'quantile_score'
    elif 'raw_score' in df.columns:
        df['abs_score'] = df['raw_score'].abs()
        idx = df.groupby('variant_id')['abs_score'].idxmax()
        agg_df = df.loc[idx].copy()
        score_col = 'raw_score'
    else:
        agg_df = df.copy()
        score_col = 'delta'
    
    valid = ~(agg_df[score_col].isna() | agg_df['beta'].isna())
    scores = agg_df.loc[valid, score_col].values
    beta = agg_df.loc[valid, 'beta'].values
    
    if len(scores) < 2:
        return
    
    # Binary labels: high vs low effect
    beta_median = np.median(np.abs(beta))
    labels = (np.abs(beta) > beta_median).astype(int)
    
    fpr, tpr, _ = roc_curve(labels, np.abs(scores))
    roc_auc = auc(fpr, tpr)
    
    plt.figure(figsize=(8, 8))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC (AUC = {roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Chance')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title(f'{dataset} - {modality}\nROC Curve', fontsize=14)
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(output_dir / f'{dataset}_{modality}_roc.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_pr(df: pd.DataFrame, output_dir: Path, dataset: str, modality: str):
    """Precision-Recall curve."""
    # Aggregate by variant
    if 'quantile_score' in df.columns:
        df['abs_score'] = df['quantile_score'].abs()
        idx = df.groupby('variant_id')['abs_score'].idxmax()
        agg_df = df.loc[idx].copy()
        score_col = 'quantile_score'
    elif 'raw_score' in df.columns:
        df['abs_score'] = df['raw_score'].abs()
        idx = df.groupby('variant_id')['abs_score'].idxmax()
        agg_df = df.loc[idx].copy()
        score_col = 'raw_score'
    else:
        agg_df = df.copy()
        score_col = 'delta'
    
    valid = ~(agg_df[score_col].isna() | agg_df['beta'].isna())
    scores = agg_df.loc[valid, score_col].values
    beta = agg_df.loc[valid, 'beta'].values
    
    if len(scores) < 2:
        return
    
    # Binary labels
    beta_median = np.median(np.abs(beta))
    labels = (np.abs(beta) > beta_median).astype(int)
    
    precision, recall, _ = precision_recall_curve(labels, np.abs(scores))
    pr_auc = auc(recall, precision)
    
    plt.figure(figsize=(8, 8))
    plt.plot(recall, precision, color='darkorange', lw=2, label=f'PR (AUC = {pr_auc:.3f})')
    plt.axhline(0.5, color='navy', lw=2, linestyle='--', label='Chance')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall', fontsize=12)
    plt.ylabel('Precision', fontsize=12)
    plt.title(f'{dataset} - {modality}\nPrecision-Recall Curve', fontsize=14)
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(output_dir / f'{dataset}_{modality}_pr.png', dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Generate plots")
    parser.add_argument('--datasets', nargs='+',
                       default=['caQTLs', 'dsQTLs', 'eQTLs', 'hQTLs'],
                       help='Datasets to plot')
    parser.add_argument('--base-dir', type=Path,
                       default=Path(__file__).parent.parent.parent,
                       help='Base directory')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing plots')
    args = parser.parse_args()
    
    base_dir = args.base_dir
    
    for dataset in args.datasets:
        pred_path = base_dir / "results/predictions" / dataset / "ag_predictions.parquet"
        
        if not pred_path.exists():
            print(f"Warning: {pred_path} not found, skipping {dataset}")
            continue
        
        print(f"Plotting {dataset}...")
        
        # Load predictions
        df = pd.read_parquet(pred_path)
        
        if len(df) == 0 or 'modality' not in df.columns:
            print(f"  Skipping {dataset} (empty or invalid)")
            continue
        
        # Plot for each modality
        for modality in df['modality'].unique():
            output_dir = base_dir / "results/plots" / dataset / modality
            output_dir.mkdir(parents=True, exist_ok=True)
            
            mod_df = df[df['modality'] == modality]
            
            print(f"  {modality}: generating plots...")
            plot_scatter(mod_df, output_dir, dataset, modality)
            plot_roc(mod_df, output_dir, dataset, modality)
            plot_pr(mod_df, output_dir, dataset, modality)
    
    print("\n✓ Plotting complete")


if __name__ == '__main__':
    main()
