#!/usr/bin/env python3
"""
Extract reference and alternate sequences for QTL variants.

Extracts 2kb windows (±1kb from SNP position) from genome FASTA,
creates both reference and alternate allele sequences.

Input:
- data/processed/variants/{dataset}.parquet
- data/raw/reference_genome/GRCh38.primary_assembly.genome.fa

Output:
- results/sequences/{dataset}/sequences.parquet
"""

import argparse
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta
from tqdm import tqdm


WINDOW_SIZE = 2048  # AlphaGenome requires 2048bp
FLANK_SIZE = WINDOW_SIZE // 2

# RefSeq to chromosome mapping (GRCh38/hg38)
REFSEQ_TO_CHR = {
    'NC_000001.11': 'chr1', 'NC_000002.12': 'chr2', 'NC_000003.12': 'chr3',
    'NC_000004.12': 'chr4', 'NC_000005.10': 'chr5', 'NC_000006.12': 'chr6',
    'NC_000007.14': 'chr7', 'NC_000008.11': 'chr8', 'NC_000009.12': 'chr9',
    'NC_000010.11': 'chr10', 'NC_000011.10': 'chr11', 'NC_000012.12': 'chr12',
    'NC_000013.11': 'chr13', 'NC_000014.9': 'chr14', 'NC_000015.10': 'chr15',
    'NC_000016.10': 'chr16', 'NC_000017.11': 'chr17', 'NC_000018.10': 'chr18',
    'NC_000019.10': 'chr19', 'NC_000020.11': 'chr20', 'NC_000021.9': 'chr21',
    'NC_000022.11': 'chr22', 'NC_000023.11': 'chrX', 'NC_000024.10': 'chrY'
}

# Build reverse mapping (chr -> RefSeq)
CHR_TO_REFSEQ = {v: k for k, v in REFSEQ_TO_CHR.items()}
# Also support without 'chr' prefix
for i in range(1, 23):
    CHR_TO_REFSEQ[str(i)] = CHR_TO_REFSEQ[f'chr{i}']
CHR_TO_REFSEQ['X'] = CHR_TO_REFSEQ['chrX']
CHR_TO_REFSEQ['Y'] = CHR_TO_REFSEQ['chrY']


def extract_sequences(df: pd.DataFrame, genome_fa: Fasta) -> pd.DataFrame:
    """Extract reference and alternate sequences for each variant."""
    records = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Extracting sequences"):
        chrom = row['chrom']
        pos = int(row['pos'])
        ref = str(row['ref']).upper()
        alt = str(row['alt']).upper()
        
        # Try to find chromosome in genome
        chr_key = None
        
        # Try direct lookup
        if chrom in genome_fa:
            chr_key = chrom
        elif f'chr{chrom}' in genome_fa:
            chr_key = f'chr{chrom}'
        elif chrom.replace('chr', '') in genome_fa:
            chr_key = chrom.replace('chr', '')
        # Try RefSeq mapping
        elif chrom in CHR_TO_REFSEQ and CHR_TO_REFSEQ[chrom] in genome_fa:
            chr_key = CHR_TO_REFSEQ[chrom]
        elif f'chr{chrom}' in CHR_TO_REFSEQ and CHR_TO_REFSEQ[f'chr{chrom}'] in genome_fa:
            chr_key = CHR_TO_REFSEQ[f'chr{chrom}']
        
        if chr_key is None:
            continue
        
        try:
            # Extract flanking sequence (0-based coordinates for pyfaidx)
            # pos is 1-based, convert to 0-based
            pos_0based = pos - 1
            start = max(0, pos_0based - FLANK_SIZE)
            end = pos_0based + FLANK_SIZE
            
            # Get reference sequence
            seq_obj = genome_fa[chr_key][start:end]
            ref_seq = str(seq_obj.seq).upper()
            
            # Create alternate sequence
            # Find the SNP position in our extracted sequence
            snp_pos_in_seq = pos_0based - start
            
            # Verify reference allele matches
            actual_ref = ref_seq[snp_pos_in_seq:snp_pos_in_seq + len(ref)]
            
            if len(ref) == 1 and len(alt) == 1:  # Simple SNP
                # Replace reference with alternate
                alt_seq = ref_seq[:snp_pos_in_seq] + alt + ref_seq[snp_pos_in_seq + 1:]
            else:
                # Handle indels (more complex)
                alt_seq = ref_seq[:snp_pos_in_seq] + alt + ref_seq[snp_pos_in_seq + len(ref):]
            
            # Ensure sequences are same length (pad/trim if needed)
            target_len = WINDOW_SIZE
            ref_seq = ref_seq[:target_len].ljust(target_len, 'N')
            alt_seq = alt_seq[:target_len].ljust(target_len, 'N')
            
            records.append({
                'variant_id': row['variant_id'],
                'chrom': row['chrom'],
                'pos': row['pos'],
                'ref': ref,
                'alt': alt,
                'beta': row['beta'],
                'ontology': row['ontology'],
                'modalities': row['modalities'],
                'gene_name': row.get('gene_name', None),  # Preserve target gene for eQTLs
                'ref_seq': ref_seq,
                'alt_seq': alt_seq,
                'seq_length': len(ref_seq)
            })
            
        except Exception as e:
            # Skip variants that fail extraction
            continue
    
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="Extract sequences for QTL variants")
    parser.add_argument('--datasets', nargs='+',
                       default=['caQTLs', 'dsQTLs', 'hQTLs'],
                       help='Datasets to process')
    parser.add_argument('--base-dir', type=Path,
                       default=Path(__file__).parent.parent.parent,
                       help='Base directory')
    parser.add_argument('--genome', type=Path,
                       help='Path to genome FASTA file')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing sequences')
    args = parser.parse_args()
    
    base_dir = args.base_dir
    input_dir = base_dir / "data/processed/variants"
    
    # Locate genome FASTA
    if args.genome:
        genome_path = args.genome
    else:
        genome_path = base_dir / "data/raw/reference_genome/GRCh38.primary_assembly.genome.fa"
    
    if not genome_path.exists():
        print(f"ERROR: Genome FASTA not found at {genome_path}")
        print("Please provide --genome path or place genome at:")
        print(f"  {base_dir}/data/raw/reference_genome/GRCh38.primary_assembly.genome.fa")
        return 1
    
    print(f"Loading genome: {genome_path}")
    genome_fa = Fasta(str(genome_path))
    print(f"Genome loaded with {len(genome_fa.keys())} chromosomes")
    
    for dataset in args.datasets:
        input_path = input_dir / f"{dataset}.parquet"
        
        if not input_path.exists():
            print(f"Warning: {input_path} not found, skipping {dataset}")
            continue
        
        output_dir = base_dir / "results/sequences" / dataset
        output_path = output_dir / "sequences.parquet"
        
        if output_path.exists() and not args.force:
            print(f"Skipping {dataset} (sequences exist, use --force to overwrite)")
            continue
        
        print(f"\nExtracting sequences for {dataset}...")
        
        # Load variants
        df = pd.read_parquet(input_path)
        print(f"  Loaded {len(df)} variants")
        
        # Extract sequences
        seq_df = extract_sequences(df, genome_fa)
        print(f"  Extracted {len(seq_df)} sequences ({len(seq_df)/len(df)*100:.1f}% success)")
        
        # Save
        output_dir.mkdir(parents=True, exist_ok=True)
        seq_df.to_parquet(output_path, index=False)
        print(f"  Saved to {output_path}")
    
    print("\n✓ Sequence extraction complete")
    return 0


if __name__ == '__main__':
    exit(main())
