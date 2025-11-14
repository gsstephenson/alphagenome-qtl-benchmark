#!/usr/bin/env python3
"""Run AlphaGenome predictions using OFFICIAL API."""

import argparse
import os
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import warnings


def load_alphagenome_client(api_key: str):
    """Initialize AlphaGenome client and load gene annotations."""
    try:
        from alphagenome.models import dna_client, variant_scorers
        from alphagenome.data import genome
        
        client = dna_client.create(api_key=api_key)
        
        print("Loading GENCODE v46 annotations...")
        gtf_path = Path(__file__).parent.parent.parent / "data/annotations/gencode.v46.annotation.gtf.gz.feather"
        
        if not gtf_path.exists():
            print("  Downloading gene annotations...")
            gtf = pd.read_feather(
                'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
            )
            gtf_path.parent.mkdir(parents=True, exist_ok=True)
            gtf.to_feather(gtf_path)
        else:
            print(f"  Loading from {gtf_path}")
            gtf = pd.read_feather(gtf_path)
        
        print(f"  Loaded {len(gtf)} gene annotations")
        
        return client, variant_scorers, genome, gtf
        
    except Exception as e:
        raise RuntimeError(f"Failed to initialize AlphaGenome: {e}")


def create_variant_scorers(variant_scorers_module, modalities):
    """Create appropriate variant scorers for each modality.
    
    Returns:
        tuple: (scorer_map dict, scorer_to_modality dict)
        - scorer_map: modality -> scorer object
        - scorer_to_modality: scorer object id -> list of modality names
    """
    scorer_map = {}
    scorer_to_modality = {}  # Maps scorer object ID to list of modalities
    recommended = variant_scorers_module.RECOMMENDED_VARIANT_SCORERS
    
    for modality in modalities:
        scorer = None
        
        if modality == 'ATAC':
            scorer = recommended['ATAC']
        elif modality == 'DNase':
            scorer = recommended['DNASE']
        elif modality == 'H3K27ac':
            scorer = recommended['CHIP_HISTONE']
        elif modality == 'H3K4me1':
            scorer = recommended['CHIP_HISTONE']
        elif modality == 'RNA_SEQ':  # NEW: Support for eQTLs
            scorer = recommended['RNA_SEQ']
        else:
            warnings.warn(f"Unknown modality: {modality}")
            continue
        
        scorer_map[modality] = scorer
        
        # Track which modalities use this scorer (for reverse lookup)
        scorer_id = id(scorer)
        if scorer_id not in scorer_to_modality:
            scorer_to_modality[scorer_id] = []
        scorer_to_modality[scorer_id].append(modality)
    
    return scorer_map, scorer_to_modality


def predict_variant_official(client, variant_scorers_module, genome_module,
                            variant_row, scorer_map, scorer_to_modality, gtf):
    """Run AlphaGenome prediction using official API.
    
    Args:
        scorer_to_modality: Dict mapping scorer object id to list of modality names
    """
    from alphagenome.models import dna_client
    
    # Ensure chromosome has 'chr' prefix
    chrom = variant_row['chrom']
    if not str(chrom).startswith('chr'):
        chrom = f"chr{chrom}"
    
    variant = genome_module.Variant(
        chromosome=chrom,
        position=int(variant_row['pos']),
        reference_bases=variant_row['ref'],
        alternate_bases=variant_row['alt']
    )
    
    variant_interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    
    # Deduplicate scorers (multiple modalities can map to same scorer)
    # e.g., H3K27ac and H3K4me1 both use CHIP_HISTONE scorer
    seen_scorer_ids = set()
    scorers_to_use = []
    for scorer in scorer_map.values():
        scorer_id = id(scorer)
        if scorer_id not in seen_scorer_ids:
            scorers_to_use.append(scorer)
            seen_scorer_ids.add(scorer_id)
    
    if len(scorers_to_use) == 0:
        return []
    
    try:
        variant_scores = client.score_variant(
            interval=variant_interval,
            variant=variant,
            variant_scorers=scorers_to_use
        )
    except Exception as e:
        print(f"\n  Warning: score_variant failed for {variant_row['variant_id']}: {e}")
        return []
    
    records = []
    
    for scorer_result in variant_scores:
        try:
            tidy_df = variant_scorers_module.tidy_scores(
                [scorer_result], 
                match_gene_strand=True
            )
        except Exception as e:
            print(f"\n  Warning: tidy_scores failed: {e}")
            continue
        
        # Get scorer string representation which contains the requested_output
        # e.g., "CenterMaskScorer(requested_output=ATAC, ...)" or "GeneMaskLFCScorer(requested_output=RNA_SEQ)"
        scorer_str = str(scorer_result.scorers[0]) if hasattr(scorer_result, 'scorers') else str(scorer_result)
        scorer_id = id(scorer_result)  # Get scorer object ID
        
        # Extract the requested_output from scorer string
        # FIXED: Parse requested_output from scorer string instead of relying on missing uns['scorer_name']
        if 'requested_output=ATAC' in scorer_str:
            modality = 'ATAC'
        elif 'requested_output=DNASE' in scorer_str:
            modality = 'DNase'
        elif 'requested_output=RNA_SEQ' in scorer_str:
            modality = 'RNA_SEQ'
        elif 'requested_output=CHIP_HISTONE' in scorer_str:
            # For histone marks, use all modalities that map to this scorer
            # This fixes the hQTLs bug where H3K27ac and H3K4me1 both use CHIP_HISTONE
            modality = ','.join(variant_row.get('modalities', '').split(','))
        else:
            modality = 'unknown'
        
        # Use all columns from tidy_df, excluding non-serializable ones
        exclude_cols = {'scored_interval'}  # Can't serialize Interval objects
        
        for idx, row in tidy_df.iterrows():
            record = {
                'variant_id': variant_row['variant_id'],
                'chrom': variant_row['chrom'],
                'pos': variant_row['pos'],
                'ref': variant_row['ref'],
                'alt': variant_row['alt'],
                'modality': modality,
                'beta': variant_row['beta']
            }
            
            # Add all columns from tidy_df except excluded ones
            for col in tidy_df.columns:
                if col not in record and col not in exclude_cols:
                    val = row[col]
                    # Convert any remaining complex objects to strings
                    if hasattr(val, '__class__') and not isinstance(val, (str, int, float, bool, type(None))):
                        val = str(val)
                    record[col] = val
            
            # Ensure key fields exist with defaults
            record.setdefault('gene_id', 'NA')
            record.setdefault('gene_name', 'NA')
            record.setdefault('tissue_ontology', 'NA')
            record.setdefault('score', np.nan)
            record.setdefault('quantile_score', np.nan)
            
            # Add target gene name for eQTLs (for filtering in evaluation)
            record['target_gene'] = variant_row.get('gene_name', None)
            
            records.append(record)
    
    return records


def generate_predictions(df, client, variant_scorers_module, genome_module, gtf, use_mock=False):
    """Generate predictions for all variants using official API."""
    records = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Predicting"):
        modalities_str = str(row['modalities'])
        modalities = [m.strip() for m in modalities_str.split(',')]
        
        # Skip scorer creation in mock mode
        if not use_mock:
            scorer_map, scorer_to_modality = create_variant_scorers(variant_scorers_module, modalities)
            if len(scorer_map) == 0:
                continue
        else:
            scorer_map = {}  # Not needed for mock
            scorer_to_modality = {}
        
        try:
            if use_mock:
                import hashlib
                for modality in modalities:
                    gene_key = f"{row['chrom']}-{row['pos']}-{modality}"
                    score = (hashlib.md5(gene_key.encode()).digest()[0] / 255.0) - 0.5
                    records.append({
                        'variant_id': row['variant_id'],
                        'chrom': row['chrom'],
                        'pos': row['pos'],
                        'ref': row['ref'],
                        'alt': row['alt'],
                        'gene_id': 'MOCK_GENE',
                        'gene_name': row.get('gene_name', 'MockGene'),  # Use actual gene name if available
                        'tissue_ontology': row['ontology'],
                        'modality': modality,
                        'score': score,
                        'quantile_score': score,
                        'beta': row['beta'],
                        'target_gene': row.get('gene_name', None)  # For eQTL filtering
                    })
            else:
                variant_records = predict_variant_official(
                    client, variant_scorers_module, genome_module,
                    row, scorer_map, scorer_to_modality, gtf
                )
                records.extend(variant_records)
        except Exception as e:
            print(f"\n  Error predicting {row['variant_id']}: {e}")
            continue
    
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="Run AlphaGenome predictions")
    parser.add_argument('--datasets', nargs='+', default=['caQTLs', 'dsQTLs', 'hQTLs'])
    parser.add_argument('--base-dir', type=Path, default=Path(__file__).parent.parent.parent)
    parser.add_argument('--api-key', type=str)
    parser.add_argument('--mock', action='store_true')
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    
    base_dir = args.base_dir
    api_key = args.api_key or os.getenv('ALPHA_GENOME_KEY')
    
    if not api_key and not args.mock:
        print("ERROR: AlphaGenome API key required!")
        return 1
    
    if not args.mock:
        print("Initializing AlphaGenome client with OFFICIAL API...")
        try:
            client, variant_scorers, genome, gtf = load_alphagenome_client(api_key)
            print("✓ AlphaGenome client ready")
        except Exception as e:
            print(f"ERROR: {e}")
            return 1
    else:
        print("Using MOCK predictions")
        client, variant_scorers, genome, gtf = None, None, None, None
    
    for dataset in args.datasets:
        input_path = base_dir / "results/sequences" / dataset / "sequences.parquet"
        
        if not input_path.exists():
            print(f"Warning: {input_path} not found, skipping {dataset}")
            continue
        
        output_dir = base_dir / "results/predictions" / dataset
        output_path = output_dir / "ag_predictions.parquet"
        
        if output_path.exists() and not args.force:
            print(f"Skipping {dataset} (use --force to overwrite)")
            continue
        
        print(f"\n{'='*60}")
        print(f"Generating predictions for {dataset}")
        print(f"{'='*60}")
        
        df = pd.read_parquet(input_path)
        
        if args.limit:
            df = df.head(args.limit)
            print(f"  Limited to {len(df)} variants")
        else:
            print(f"  Loaded {len(df)} variants")
        
        pred_df = generate_predictions(df, client, variant_scorers, genome, gtf, use_mock=args.mock)
        
        print(f"\n  Generated {len(pred_df)} gene-level predictions")
        
        if len(pred_df) > 0:
            print(f"  Unique genes: {pred_df['gene_id'].nunique()}")
            print(f"  Modalities: {pred_df['modality'].unique()}")
            print(f"  Score range: [{pred_df['score'].min():.4f}, {pred_df['score'].max():.4f}]")
        
        output_dir.mkdir(parents=True, exist_ok=True)
        pred_df.to_parquet(output_path, index=False)
        print(f"  Saved to {output_path}")
    
    print("\n" + "="*60)
    print("✓ Prediction generation complete using OFFICIAL API")
    print("="*60)
    return 0


if __name__ == '__main__':
    exit(main())
