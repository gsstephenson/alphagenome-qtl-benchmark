# AlphaGenome QTL Benchmark

**A rigorous evaluation of AlphaGenome's ability to predict quantitative trait loci (QTL) effects using the official AlphaGenome API.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)

---

## Overview

This repository contains a systematic benchmark of Google DeepMind's AlphaGenome model for predicting the functional impact of genetic variants on gene regulation. We test AlphaGenome on **16,072 real quantitative trait loci (QTLs)** across four regulatory modalities:

- **eQTLs** (424 variants): Gene expression QTLs
- **caQTLs** (3,317 variants): Chromatin accessibility QTLs  
- **dsQTLs** (6,070 variants): DNase sensitivity QTLs
- **hQTLs** (6,261 variants): Histone modification QTLs (H3K27ac, H3K4me1)

### Key Question

**Can AlphaGenome, trained on DNA sequence alone, accurately predict the direction and magnitude of regulatory variant effects measured in real human cells?**

---

## Methodology

### AlphaGenome API Usage

We use AlphaGenome's **official high-level API** exactly as recommended in their documentation:

```python
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

# Initialize client
client = dna_client.create(api_key=API_KEY)

# Create variant
variant = genome.Variant(
    chromosome="chr1",
    position=12345,
    reference_bases="A",
    alternate_bases="G"
)

# Score variant with recommended scorers
variant_scores = client.score_variant(
    interval=variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
    variant=variant,
    variant_scorers=[
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'],  # for eQTLs
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['ATAC'],     # for caQTLs
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['DNASE'],    # for dsQTLs
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_HISTONE']  # for hQTLs
    ]
)

# Extract gene-level predictions
tidy_df = variant_scorers.tidy_scores([variant_scores], match_gene_strand=True)
```

**Key Design Choices:**

1. **1MB genomic windows**: Standard AlphaGenome input size (`SEQUENCE_LENGTH_1MB`)
2. **Gene-level aggregation**: Using `tidy_scores(match_gene_strand=True)` for proper gene matching
3. **Recommended scorers**: Using `RECOMMENDED_VARIANT_SCORERS` exactly as specified by AlphaGenome
4. **Target gene filtering** (eQTLs only): Predictions filtered to match the actual gene affected by each eQTL
5. **GENCODE v46 annotations**: For accurate gene identification

### Dataset Selection Rationale

| Dataset | Source | N | Why This Dataset? |
|---------|--------|---|-------------------|
| **eQTLs** | GSE86886 | 424 | Gold standard for gene expression; AlphaGenome paper benchmarked on eQTLs (reported r=0.49) |
| **caQTLs** | GSE86886 | 3,317 | ATAC-seq QTLs from same study; tests chromatin accessibility predictions |
| **dsQTLs** | GSE31388 | 6,070 | DNase sensitivity QTLs; independent dataset for robustness |
| **hQTLs** | GSE116193 | 6,261 | Histone modification QTLs; tests epigenetic mark predictions |

**Cell Type:** All QTLs from CD4+ T cells (primary human immune cells)

---

## Repository Structure

```
alphagenome_phase3_qtl_benchmark/
├── README.md                          # This file
├── IMPLEMENTATION_SUMMARY.md          # Technical implementation details
├── QTL_BENCHMARK_RESULTS.md           # Full results and analysis
├── requirements.txt                   # Python dependencies
│
├── data/
│   ├── raw/                          # Original QTL datasets (not in git)
│   │   ├── caQTLs_GSE86886/
│   │   ├── dsQTLs_GSE31388/
│   │   ├── hQTLs_GSE116193/
│   │   └── reference_genome/
│   ├── processed/                    # Normalized QTL data (not in git)
│   │   └── variants/
│   ├── annotations/                  # GENCODE gene annotations (not in git)
│   └── metadata/
│       ├── README.md                 # Data provenance documentation
│       ├── sample_metadata.tsv
│       └── tissue_mappings.tsv
│
├── scripts/
│   └── real/                         # Production pipeline scripts
│       ├── 00_normalize.py           # Convert raw QTL data to standard format
│       ├── 01_extract_sequences.py   # Extract genomic sequences
│       ├── 02_predict.py             # Run AlphaGenome predictions
│       ├── 03_evaluate.py            # Compute correlation metrics
│       ├── 04_plots.py               # Generate diagnostic plots
│       └── run_all.py                # Master pipeline driver
│
├── results/                          # Pipeline outputs (not in git)
│   ├── predictions/                  # AlphaGenome predictions by dataset
│   ├── sequences/                    # Extracted genomic sequences
│   ├── tables/                       # Evaluation metrics
│   └── plots/                        # Correlation and diagnostic plots
│
└── logs/                             # Execution logs (not in git)
```

---

## Installation

### Requirements

- Python 3.11+
- AlphaGenome API key (request from Google DeepMind)
- ~50GB disk space for reference genome and processed data

### Setup

```bash
# Clone repository
git clone https://github.com/yourusername/alphagenome_phase3_qtl_benchmark.git
cd alphagenome_phase3_qtl_benchmark

# Create conda environment
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env

# Install dependencies
pip install -r requirements.txt

# Set AlphaGenome API key
export ALPHA_GENOME_KEY="your_api_key_here"
```

### Download Reference Data

```bash
# Download GRCh38 reference genome
mkdir -p data/raw/reference_genome
cd data/raw/reference_genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa
```

**Note:** Raw QTL datasets are already included in `data/raw/` directories.

---

## Usage

### Quick Start (Mock Mode for Testing)

Test the pipeline with mock predictions (no API calls):

```bash
python3 scripts/real/run_all.py --datasets caQTLs --mock --limit 100
```

### Full Pipeline (Real AlphaGenome API)

Run the complete benchmark on all datasets:

```bash
# Ensure API key is set
export ALPHA_GENOME_KEY="your_api_key_here"

# Run full pipeline
python3 scripts/real/run_all.py --datasets all
```

**⚠️ Warning:** This will make ~16,000 API calls and may take several hours.

### Step-by-Step Execution

Run individual pipeline steps:

```bash
# Step 1: Normalize QTL data
python3 scripts/real/00_normalize.py --datasets caQTLs eQTLs dsQTLs hQTLs

# Step 2: Extract genomic sequences
python3 scripts/real/01_extract_sequences.py --datasets caQTLs eQTLs dsQTLs hQTLs

# Step 3: Run AlphaGenome predictions
python3 scripts/real/02_predict.py --datasets caQTLs eQTLs dsQTLs hQTLs

# Step 4: Evaluate predictions
python3 scripts/real/03_evaluate.py --datasets caQTLs eQTLs dsQTLs hQTLs

# Step 5: Generate plots
python3 scripts/real/04_plots.py --datasets caQTLs eQTLs dsQTLs hQTLs
```

---

## Expected Outputs

### Predictions

Each variant generates gene-level predictions:

```
results/predictions/{dataset}/ag_predictions.parquet
```

Columns:
- `variant_id`: Unique variant identifier
- `chrom`, `pos`, `ref`, `alt`: Variant coordinates
- `gene_id`, `gene_name`: Target gene annotations
- `quantile_score`: Normalized effect prediction (-1 to +1)
- `score`: Raw log-fold change prediction
- `tissue_ontology`: Cell type
- `beta`: Observed QTL effect size (ground truth)

### Evaluation Metrics

```
results/tables/{dataset}_metrics.tsv
```

Metrics per modality:
- **Spearman r**: Correlation between predicted and observed effect sizes
- **Pearson r**: Linear correlation
- **AUROC**: Binary classification accuracy (high vs. low effect)
- **Sign concordance**: Directional accuracy (% correct sign)

### Plots

```
results/plots/{dataset}/{modality}/
├── scatter_plot.png           # Predicted vs. observed effects
├── roc_curve.png              # ROC curve for classification
└── residuals.png              # Residual analysis
```

---

## Expected Results

Based on the AlphaGenome paper (Nature, 2024):

| QTL Type | Expected Spearman r | Published Result |
|----------|-------------------|------------------|
| **eQTLs** | ~0.49 | Figure 4c,d (paper benchmark) |
| **caQTLs** | Unknown | Not reported in paper |
| **dsQTLs** | Unknown | Not reported in paper |
| **hQTLs** | Unknown | Not reported in paper |

---

## Implementation Details

### Critical Design Decisions

1. **Gene Matching for eQTLs**
   - AlphaGenome returns predictions for *all* genes in 1MB window
   - We filter to only the **target gene** specified in eQTL data
   - This is critical for fair evaluation (matching paper methodology)

2. **Scorer Selection**
   - eQTLs → `GeneMaskLFCScorer(RNA_SEQ)`
   - caQTLs → `GeneMaskLFCScorer(ATAC)` + `CenterMaskScorer(CHIP_HISTONE)`
   - dsQTLs → `GeneMaskLFCScorer(DNASE)`
   - hQTLs → `CenterMaskScorer(CHIP_HISTONE)`

3. **Scorer Deduplication**
   - H3K27ac and H3K4me1 both map to `CHIP_HISTONE` scorer
   - We deduplicate to avoid redundant API calls

4. **Variant Aggregation**
   - Multiple genes per variant: take gene with max |quantile_score|
   - Ensures one prediction per variant for correlation analysis

See `IMPLEMENTATION_SUMMARY.md` for full technical details.

---

## Data Provenance

All QTL datasets are from published studies with public GEO accession numbers:

- **GSE86886**: CD4+ T cell caQTLs and eQTLs (Raj et al.)
- **GSE31388**: CD4+ T cell dsQTLs (Degner et al.)
- **GSE116193**: CD4+ T cell hQTLs (Pelikan et al.)

See `data/metadata/README.md` for detailed sample information and quality control metrics.

---

## Citation

If you use this benchmark in your research, please cite:

```bibtex
@software{stephenson2025alphagenome_qtl,
  author = {Stephenson, Grant},
  title = {AlphaGenome QTL Benchmark},
  year = {2025},
  url = {https://github.com/yourusername/alphagenome_phase3_qtl_benchmark}
}
```

And cite the original AlphaGenome paper:

```bibtex
@article{alphagenome2024,
  title={AlphaGenome: Predicting gene expression and regulatory effects},
  author={Google DeepMind},
  journal={Nature},
  year={2024}
}
```

---

## License

MIT License - see LICENSE file for details.

---

## Contact

**Grant Stephenson**  
Layer Laboratory  
University of Colorado Boulder  
grant.stephenson@colorado.edu

---

## Acknowledgments

- **Layer Laboratory** for computational resources
- **Google DeepMind** for AlphaGenome API access
- Original QTL study authors for public data release
