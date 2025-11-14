#!/usr/bin/env python3
"""
Normalize raw QTL data from multiple sources into standardized parquet format.

MAJOR UPDATE: Now extracts peak coordinates for caQTLs and adds eQTL support!

Input:
- caQTLs: data/raw/caQTLs_GSE86886/ATAC-QTLs.csv (with peak coords)
- eQTLs: data/raw/caQTLs_GSE86886/eQTLs.csv (NEW!)
- dsQTLs: data/raw/dsQTLs_GSE31388/GSE31388_dsQtlTable.txt
- hQTLs: data/raw/hQTLs_GSE116193/Pelikan_et_al_hQTL_summary.csv

Output:
- data/processed/variants/{dataset}.parquet
- Columns: dataset, variant_id, chrom, pos, ref, alt, beta, ontology, modalities, 
           activation, qtl_type, peak_chrom, peak_start, peak_end, gene_name
"""

import argparse
import pandas as pd
from pathlib import Path
from datetime import datetime


def load_caQTLs(base_dir: Path) -> pd.DataFrame:
    """Load and normalize caQTLs from GSE86886 with PEAK COORDINATES."""
    csv_path = base_dir / "data/raw/caQTLs_GSE86886/ATAC-QTLs.csv"
    df = pd.read_csv(csv_path)
    
    records = []
    for idx, row in df.iterrows():
        # Extract chromosome and position from SNP column
        # Format: chr12:9436083_rs61916194
        snp = str(row.get('SNP', ''))
        
        # Split by underscore first to remove rs ID
        if '_' in snp:
            snp = snp.split('_')[0]
        
        # Parse variant position
        if ':' in snp:
            parts = snp.split(':')
            chrom = parts[0].replace('chr', '')
            try:
                pos = int(parts[1]) if len(parts) > 1 else idx
            except (ValueError, IndexError):
                pos = idx
        else:
            chrom = 'unknown'
            pos = idx
        
        # CRITICAL FIX: Extract peak coordinates from "gene" column
        # Format: chr12:9435883-9437556
        gene_field = str(row.get('gene', ''))
        peak_chrom, peak_start, peak_end = None, None, None
        
        if gene_field and ':' in gene_field and '-' in gene_field:
            try:
                peak_parts = gene_field.split(':')
                peak_chrom = peak_parts[0].replace('chr', '')
                coords = peak_parts[1].split('-')
                peak_start = int(coords[0])
                peak_end = int(coords[1])
            except (ValueError, IndexError):
                pass  # Keep as None if parsing fails
        
        # Get beta
        beta = float(row.get('beta', row.get('Beta', 0)))
        
        # Determine activation from file name (would need to pass this info)
        activation = ''  # Could be '0hr' or '48hr'
        
        records.append({
            'dataset': 'caQTLs',
            'variant_id': f"caQTL_{idx}",
            'chrom': chrom,
            'pos': pos,
            'ref': 'A',  # placeholder - would need VCF to get real alleles
            'alt': 'G',
            'beta': beta,
            'ontology': 'CD4_Tcell',
            'modalities': 'ATAC,DNase,H3K27ac',
            'activation': activation,
            'qtl_type': 'chromatin_accessibility',
            'peak_chrom': peak_chrom,
            'peak_start': peak_start,
            'peak_end': peak_end,
            'gene_name': None  # caQTLs are peak-based, not gene-based
        })
    
    return pd.DataFrame(records)


def load_eQTLs(base_dir: Path) -> pd.DataFrame:
    """Load and normalize eQTLs from GSE86886 (NEW!)."""
    csv_path = base_dir / "data/raw/caQTLs_GSE86886/eQTLs.csv"
    
    if not csv_path.exists():
        print(f"  Warning: {csv_path} not found, skipping eQTLs")
        return pd.DataFrame()
    
    df = pd.read_csv(csv_path)
    
    records = []
    for idx, row in df.iterrows():
        # Extract chromosome and position from SNP column
        # Format: chr5:96322136_rs38033
        snp = str(row.get('SNP', ''))
        
        # Skip empty rows
        if pd.isna(row.get('SNP')) or snp == '' or snp == 'nan':
            continue
        
        # Split by underscore first to remove rs ID
        if '_' in snp:
            snp = snp.split('_')[0]
        
        # Parse variant position
        if ':' in snp:
            parts = snp.split(':')
            chrom = parts[0].replace('chr', '')
            try:
                pos = int(parts[1]) if len(parts) > 1 else idx
            except (ValueError, IndexError):
                continue  # Skip malformed entries
        else:
            continue  # Skip entries without proper chromosome:position format
        
        # Get gene name - this is the TARGET for prediction!
        gene_name = str(row.get('gene', None))
        
        # Get beta
        beta = float(row.get('beta', row.get('Beta', 0)))
        
        records.append({
            'dataset': 'eQTLs',
            'variant_id': f"eQTL_{idx}",
            'chrom': chrom,
            'pos': pos,
            'ref': 'A',  # placeholder
            'alt': 'G',
            'beta': beta,
            'ontology': 'CD4_Tcell',  # Same source as caQTLs
            'modalities': 'RNA_SEQ',  # eQTLs measure gene expression!
            'activation': '',
            'qtl_type': 'gene_expression',
            'peak_chrom': None,
            'peak_start': None,
            'peak_end': None,
            'gene_name': gene_name  # This is what we predict!
        })
    
    return pd.DataFrame(records)


def load_dsQTLs(base_dir: Path) -> pd.DataFrame:
    """Load and normalize dsQTLs from GSE31388."""
    txt_path = base_dir / "data/raw/dsQTLs_GSE31388/GSE31388_dsQtlTable.txt"
    
    if txt_path.exists():
        df = pd.read_csv(txt_path, sep='\t')
    else:
        txt_path = base_dir / "data/raw/dsQTLs_GSE31388/GSE31388_dsQtlTableLong.txt"
        df = pd.read_csv(txt_path, sep='\t')
    
    records = []
    for idx, row in df.iterrows():
        snp = str(row.get('SNP', row.get('snp', f'ds_{idx}')))
        
        chrom = str(row.get('Chr', 'unknown')).replace('chr', '')
        
        try:
            pos = int(snp)
        except (ValueError, TypeError):
            pos = int(row.get('Start', idx))
        
        beta = float(row.get('Estimate', row.get('beta', row.get('effect', 0))))
        
        records.append({
            'dataset': 'dsQTLs',
            'variant_id': f"dsQTL_{idx}",
            'chrom': chrom,
            'pos': pos,
            'ref': 'C',
            'alt': 'T',
            'beta': beta,
            'ontology': 'LCL',
            'modalities': 'DNase,ATAC',
            'activation': '',
            'qtl_type': 'dnase_sensitivity',
            'peak_chrom': None,  # Would need peak data
            'peak_start': None,
            'peak_end': None,
            'gene_name': None
        })
    
    return pd.DataFrame(records)


def load_hQTLs(base_dir: Path) -> pd.DataFrame:
    """Load and normalize hQTLs from GSE116193."""
    csv_path = base_dir / "data/raw/hQTLs_GSE116193/Pelikan_et_al_hQTL_summary.csv"
    df = pd.read_csv(csv_path)
    
    records = []
    for idx, row in df.iterrows():
        chrom = str(row.get('Chr', 'unknown')).replace('chr', '')
        pos = int(row.get('Bp (hg19)', idx))
        
        beta = float(row.get('log2 (effect size)', row.get('beta', row.get('effect_size', 0))))
        
        ref = str(row.get('Ref Allelec', 'G'))
        alt = str(row.get('Alt Alleled', 'A'))
        
        records.append({
            'dataset': 'hQTLs',
            'variant_id': f"hQTL_{idx}",
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'beta': beta,
            'ontology': 'B_cell',
            'modalities': 'H3K27ac,H3K4me1',
            'activation': '',
            'qtl_type': 'histone_modification',
            'peak_chrom': None,  # Would need peak data
            'peak_start': None,
            'peak_end': None,
            'gene_name': None
        })
    
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="Normalize QTL datasets with peak coordinates")
    parser.add_argument('--datasets', nargs='+', 
                       default=['caQTLs', 'dsQTLs', 'hQTLs', 'eQTLs'],
                       help='Datasets to process')
    parser.add_argument('--base-dir', type=Path, 
                       default=Path(__file__).parent.parent.parent,
                       help='Base directory')
    args = parser.parse_args()
    
    base_dir = args.base_dir
    output_dir = base_dir / "data/processed/variants"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    loaders = {
        'caQTLs': load_caQTLs,
        'eQTLs': load_eQTLs,  # NEW!
        'dsQTLs': load_dsQTLs,
        'hQTLs': load_hQTLs
    }
    
    print(f"\n{'='*60}")
    print("NORMALIZING QTL DATASETS (WITH PEAK COORDINATES)")
    print(f"{'='*60}\n")
    
    dataset_counts = {}
    
    for dataset in args.datasets:
        if dataset not in loaders:
            print(f"Warning: Unknown dataset '{dataset}', skipping")
            continue
        
        print(f"Processing {dataset}...")
        try:
            df = loaders[dataset](base_dir)
            
            if len(df) == 0:
                print(f"  Skipped (no data)")
                continue
            
            output_path = output_dir / f"{dataset}.parquet"
            df.to_parquet(output_path, index=False)
            dataset_counts[dataset] = len(df)
            
            print(f"  ✓ Wrote {len(df)} variants to {output_path}")
            
            # Show summary
            if 'qtl_type' in df.columns:
                print(f"    QTL type: {df['qtl_type'].iloc[0]}")
            if 'peak_start' in df.columns and df['peak_start'].notna().any():
                n_peaks = df['peak_start'].notna().sum()
                print(f"    Peak coordinates: {n_peaks}/{len(df)} variants")
            if 'gene_name' in df.columns and df['gene_name'].notna().any():
                n_genes = df['gene_name'].notna().sum()
                print(f"    Gene names: {n_genes}/{len(df)} variants")
            
        except Exception as e:
            print(f"  ✗ Error processing {dataset}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'='*60}")
    print(f"✓ Normalization complete: {len(dataset_counts)} datasets")
    print(f"{'='*60}\n")
    
    for dataset, count in dataset_counts.items():
        print(f"  {dataset}: {count:,} variants")


if __name__ == '__main__':
    main()
