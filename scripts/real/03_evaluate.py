#!/usr/bin/env python3
"""
Evaluate AlphaGenome predictions against observed QTL effects.

This script computes:
- Spearman and Pearson correlations between predicted and observed effect sizes
- AUROC and AUPRC for binary classification (high vs low effect)
- Sign concordance (directional accuracy)

Input:
- results/predictions/{dataset}/ag_predictions.parquet

Output:
- results/tables/{dataset}_metrics.tsv
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score, average_precision_score


def compute_metrics(df: pd.DataFrame, dataset: str = None) -> dict:
    """Compute evaluation metrics.
    
    Args:
        df: DataFrame with predictions
        dataset: Dataset name (e.g., 'eQTLs') for special handling
    """
    # CRITICAL FIX: For eQTLs, filter to target gene matches only
    if dataset == 'eQTLs' and 'target_gene' in df.columns and 'gene_name' in df.columns:
        # Only keep predictions where predicted gene matches target gene
        before = len(df)
        df = df[df['gene_name'] == df['target_gene']].copy()
        after = len(df)
        print(f"    Filtered to target genes: {after}/{before} predictions ({after/before*100:.1f}%)")
    
    # Handle new format with gene-level scores
    # Aggregate by variant: use max absolute quantile_score per variant
    if 'quantile_score' in df.columns:
        # Group by variant and take row with max absolute quantile_score
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
        # Old format with single delta per variant
        agg_df = df.copy()
        score_col = 'delta'
    
    # Filter out NaN values
    valid = ~(agg_df[score_col].isna() | agg_df['beta'].isna())
    scores = agg_df.loc[valid, score_col].values
    beta = agg_df.loc[valid, 'beta'].values
    
    if len(scores) < 2:
        return {
            'n_variants': len(scores),
            'sign_concordance': 0.0,
            'spearman_r': 0.0,
            'spearman_p': 1.0,
            'pearson_r': 0.0,
            'pearson_p': 1.0,
            'auroc': 0.5,
            'auprc': 0.5
        }
    
    # Sign concordance
    sign_conc = np.mean(np.sign(scores) == np.sign(beta))
    
    # Correlations
    spearman_r, spearman_p = spearmanr(scores, beta)
    pearson_r, pearson_p = pearsonr(scores, beta)
    
    # Binary classification: |beta| > median
    beta_median = np.median(np.abs(beta))
    labels = (np.abs(beta) > beta_median).astype(int)
    
    # AUROC and AUPRC
    try:
        auroc = roc_auc_score(labels, np.abs(scores))
        auprc = average_precision_score(labels, np.abs(scores))
    except ValueError:
        auroc = 0.5
        auprc = 0.5
    
    return {
        'n_variants': len(scores),
        'sign_concordance': sign_conc,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'auroc': auroc,
        'auprc': auprc
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate predictions")
    parser.add_argument('--datasets', nargs='+',
                       default=['caQTLs', 'dsQTLs', 'eQTLs', 'hQTLs'],
                       help='Datasets to evaluate')
    parser.add_argument('--base-dir', type=Path,
                       default=Path(__file__).parent.parent.parent,
                       help='Base directory')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing metrics')
    args = parser.parse_args()
    
    base_dir = args.base_dir
    output_dir = base_dir / "results/tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for dataset in args.datasets:
        pred_path = base_dir / "results/predictions" / dataset / "ag_predictions.parquet"
        
        if not pred_path.exists():
            print(f"Warning: {pred_path} not found, skipping {dataset}")
            continue
        
        output_path = output_dir / f"{dataset}_metrics.tsv"
        
        if output_path.exists() and not args.force:
            print(f"Skipping {dataset} (metrics exist, use --force to overwrite)")
            continue
        
        print(f"Evaluating {dataset}...")
        
        # Load predictions
        df = pd.read_parquet(pred_path)
        print(f"  Loaded {len(df)} predictions")
        
        # Skip if empty or invalid
        if len(df) == 0 or 'modality' not in df.columns:
            print(f"  Skipping {dataset} (empty or invalid predictions)")
            continue
        
        # Compute metrics per modality
        results = []
        for modality in df['modality'].unique():
            mod_df = df[df['modality'] == modality]
            print(f"  {modality}:")
            metrics = compute_metrics(mod_df, dataset=dataset)  # Pass dataset for eQTL filtering
            metrics['dataset'] = dataset
            metrics['modality'] = modality
            results.append(metrics)
            
            print(f"    n={metrics['n_variants']}, "
                  f"r={metrics['spearman_r']:.3f}, "
                  f"auroc={metrics['auroc']:.3f}")
        
        # Save metrics
        metrics_df = pd.DataFrame(results)
        metrics_df.to_csv(output_path, sep='\t', index=False)
        print(f"  Saved to {output_path}")
    
    print("\n✓ Evaluation complete")


if __name__ == '__main__':
    main()
