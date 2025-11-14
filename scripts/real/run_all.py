#!/usr/bin/env python3
"""
Driver script for REAL AlphaGenome QTL benchmark pipeline.

Steps:
1. Normalize raw QTL data
2. Extract sequences from genome
3. Run AlphaGenome predictions
4. Evaluate predictions
5. Generate diagnostic plots

Usage:
    # With real AlphaGenome API
    export ALPHA_GENOME_KEY="your_key_here"
    python scripts/real/run_all.py --datasets all
    
    # With mock predictions for testing
    python scripts/real/run_all.py --datasets caQTLs --mock --limit 100
"""

import argparse
import subprocess
import sys
from pathlib import Path
import time
import os


DATASETS = ['caQTLs', 'eQTLs', 'dsQTLs', 'hQTLs']  # Added eQTLs!


def run_step(script_name: str, args: argparse.Namespace, step_num: int,
             extra_args: list = None) -> bool:
    """Run a pipeline step and return success status."""
    script_path = Path(__file__).parent / script_name
    
    print(f"\n{'='*70}")
    print(f"STEP {step_num}: {script_name}")
    print(f"{'='*70}")
    
    # Use conda run to ensure proper environment
    # Use -u for unbuffered output so progress bars show in real-time
    cmd = ['conda', 'run', '-n', 'alphagenome-env', 'python', '-u', str(script_path)]
    
    # Add dataset arguments
    if args.datasets != ['all']:
        cmd.extend(['--datasets'] + args.datasets)
    else:
        cmd.extend(['--datasets'] + DATASETS)
    
    # Add base dir
    cmd.extend(['--base-dir', str(args.base_dir)])
    
    # Add force flag (but never for normalization scripts)
    if args.force and not script_name.startswith('00_normalize'):
        cmd.append('--force')
    
    # Add extra arguments
    if extra_args:
        cmd.extend(extra_args)
    
    print(f"Running: {' '.join(cmd)}\n")
    
    start_time = time.time()
    # Stream output in real-time so progress bars are visible
    result = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr)
    elapsed = time.time() - start_time
    
    print(f"\nCompleted in {elapsed:.1f}s")
    
    if result.returncode != 0:
        print(f"ERROR: {script_name} failed with return code {result.returncode}")
        return False
    
    return True


def print_summary(base_dir: Path, datasets: list):
    """Print final summary with output locations."""
    print(f"\n{'='*70}")
    print("FINAL SUMMARY")
    print(f"{'='*70}\n")
    
    for dataset in datasets:
        print(f"Dataset: {dataset}")
        print(f"  Variants: data/processed/variants/{dataset}.parquet")
        print(f"  Sequences: results/sequences/{dataset}/sequences.parquet")
        print(f"  Predictions: results/predictions/{dataset}/ag_predictions.parquet")
        print(f"  Metrics: results/tables/{dataset}_metrics.tsv")
        
        # List plot directories
        plot_dir = base_dir / "results/plots" / dataset
        if plot_dir.exists():
            modalities = [d.name for d in plot_dir.iterdir() if d.is_dir()]
            if modalities:
                print(f"  Plots: results/plots/{dataset}/ ({', '.join(modalities)})")
        print()
    
    print(f"Data audit: data/metadata/README.md")
    print(f"\n{'='*70}")
    print("✓ Pipeline complete!")
    print(f"{'='*70}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Run REAL AlphaGenome QTL benchmark pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with real AlphaGenome (requires API key)
  export ALPHA_GENOME_KEY="your_key"
  python scripts/real/run_all.py --datasets all
  
  # Test with mock predictions and limited variants
  python scripts/real/run_all.py --datasets caQTLs --mock --limit 100
  
  # Skip steps
  python scripts/real/run_all.py --skip-normalize --skip-sequences
        """
    )
    
    parser.add_argument('--datasets', nargs='+', default=['all'],
                       help='Datasets to process: all, caQTLs, dsQTLs, hQTLs')
    parser.add_argument('--base-dir', type=Path,
                       default=Path(__file__).parent.parent.parent,
                       help='Base directory (default: repository root)')
    parser.add_argument('--genome', type=Path,
                       help='Path to genome FASTA file')
    parser.add_argument('--api-key', type=str,
                       help='AlphaGenome API key (or set ALPHA_GENOME_KEY env var)')
    parser.add_argument('--mock', action='store_true',
                       help='Use mock predictions (no API calls)')
    parser.add_argument('--limit', type=int,
                       help='Limit variants per dataset (for testing)')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing outputs')
    parser.add_argument('--skip-normalize', action='store_true',
                       help='Skip normalization step')
    parser.add_argument('--skip-sequences', action='store_true',
                       help='Skip sequence extraction')
    parser.add_argument('--skip-predict', action='store_true',
                       help='Skip prediction step')
    parser.add_argument('--skip-evaluate', action='store_true',
                       help='Skip evaluation step')
    parser.add_argument('--skip-plots', action='store_true',
                       help='Skip plot generation')
    
    args = parser.parse_args()
    
    # Expand 'all' to full dataset list
    if 'all' in args.datasets:
        args.datasets = DATASETS
    
    # Check API key if not using mock
    if not args.mock and not args.skip_predict:
        api_key = args.api_key or os.getenv('ALPHA_GENOME_KEY')
        if not api_key:
            print("\nWARNING: No AlphaGenome API key found!")
            print("Set ALPHA_GENOME_KEY environment variable or use --api-key")
            print("Or use --mock for testing without API calls")
            print()
            response = input("Continue with mock predictions? [y/N]: ")
            if response.lower() != 'y':
                return 1
            args.mock = True
    
    print(f"\n{'='*70}")
    print("REAL ALPHAGENOME QTL BENCHMARK PIPELINE")
    print(f"{'='*70}")
    print(f"Base directory: {args.base_dir}")
    print(f"Datasets: {', '.join(args.datasets)}")
    print(f"Mode: {'MOCK predictions' if args.mock else 'REAL AlphaGenome API'}")
    if args.limit:
        print(f"Limit: {args.limit} variants per dataset")
    print(f"Force overwrite: {args.force}")
    print(f"{'='*70}\n")
    
    # Pipeline steps
    step_num = 1
    
    # Step 1: Normalize (with peak coordinates + eQTLs)
    if not args.skip_normalize:
        success = run_step('00_normalize.py', args, step_num)
        if not success:
            print("\n✗ Pipeline failed!")
            return 1
        step_num += 1
    
    # Step 2: Extract sequences
    if not args.skip_sequences:
        extra_args = []
        if args.genome:
            extra_args.extend(['--genome', str(args.genome)])
        
        success = run_step('01_extract_sequences.py', args, step_num, extra_args)
        if not success:
            print("\n✗ Pipeline failed!")
            return 1
        step_num += 1
    
    # Step 3: Predict
    if not args.skip_predict:
        extra_args = []
        if args.mock:
            extra_args.append('--mock')
        if args.api_key:
            extra_args.extend(['--api-key', args.api_key])
        if args.limit:
            extra_args.extend(['--limit', str(args.limit)])
        
        # Use corrected prediction script with fixed modality detection
        success = run_step('02_predict.py', args, step_num, extra_args)
        if not success:
            print("\n✗ Pipeline failed!")
            return 1
        step_num += 1
    
    # Step 4: Evaluate
    if not args.skip_evaluate:
        success = run_step('03_evaluate.py', args, step_num)
        if not success:
            print("\n✗ Pipeline failed!")
            return 1
        step_num += 1
    
    # Step 5: Plots
    if not args.skip_plots:
        success = run_step('04_plots.py', args, step_num)
        if not success:
            print("\n✗ Pipeline failed!")
            return 1
    
    # Print summary
    print_summary(args.base_dir, args.datasets)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
