# AlphaGenome QTL Benchmark

**Systematic evaluation of AlphaGenome's ability to predict genetic variant effects on gene regulation using real human QTL data.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)

**Status:** Pipeline complete, ready for final predictions  
**Last Updated:** November 14, 2025

---

## 🎯 Project Overview

This repository benchmarks Google DeepMind's AlphaGenome model on **15,647 real genetic variants** affecting gene regulation in human cells. We test whether AlphaGenome can predict:

1. **Direction**: Does a variant increase or decrease gene activity?
2. **Magnitude**: How large is the regulatory effect?
3. **Cell-type specificity**: Do predictions match tissue-specific effects?

### Datasets Tested

| QTL Type | N Variants | What They Affect | Source |
|----------|-----------|------------------|--------|
| **eQTLs** | 424 | Gene expression (RNA-seq) | GSE86886 |
| **caQTLs** | 3,317 | Chromatin accessibility (ATAC-seq) | GSE86886 |
| **dsQTLs** | 6,070 | DNase sensitivity | GSE31388 |
| **hQTLs** | 6,261 | Histone modifications (H3K27ac/H3K4me1) | GSE116193 |

All from **CD4+ T cells** (primary human immune cells)

---

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/gsstephenson/alphagenome-qtl-benchmark.git
cd alphagenome-qtl-benchmark

# Create environment
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env
pip install -r requirements.txt

# Set API key
export ALPHA_GENOME_KEY="your_api_key_here"
```

### Run Pipeline

```bash
# Full pipeline on all datasets
python scripts/real/run_all.py --datasets all

# Individual dataset (faster)
python scripts/real/run_all.py --datasets eQTLs

# Test mode (no API calls)
python scripts/real/run_all.py --datasets caQTLs --mock --limit 100
```

---

## 📊 Methodology

### AlphaGenome API Usage

We use the **official high-level API** exactly as recommended:

```python
from alphagenome.models import dna_client, variant_scorers

# Initialize client
client = dna_client.create(api_key=API_KEY)

# Score variant with recommended scorers
variant_scores = client.score_variant(
    interval=variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
    variant=variant,
    variant_scorers=[
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'],      # eQTLs
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['ATAC'],         # caQTLs
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['DNASE'],        # dsQTLs
        variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_HISTONE']  # hQTLs
    ]
)

# Get gene-level predictions
tidy_df = variant_scorers.tidy_scores(variant_scores, match_gene_strand=True)
```

### Key Design Choices

1. **1MB genomic windows**: Standard AlphaGenome input (`SEQUENCE_LENGTH_1MB`)
2. **Gene-level aggregation**: Using `tidy_scores()` for proper gene matching
3. **GENCODE v46 annotations**: Accurate gene identification
4. **Target gene filtering (eQTLs)**: Match specific affected gene, not all genes in window
5. **Recommended scorers**: Using `RECOMMENDED_VARIANT_SCORERS` per AlphaGenome docs

---

## 📁 Repository Structure

```
alphagenome-qtl-benchmark/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
│
├── data/
│   ├── raw/                          # Original QTL datasets
│   │   ├── caQTLs_GSE86886/         # ATAC-seq QTLs
│   │   ├── dsQTLs_GSE31388/         # DNase-seq QTLs
│   │   ├── hQTLs_GSE116193/         # Histone ChIP-seq QTLs
│   │   └── reference_genome/        # GRCh38 + GENCODE
│   ├── annotations/                  # GENCODE v46 gene annotations
│   └── metadata/                     # Sample info & tissue mappings
│
├── scripts/real/                     # Production pipeline
│   ├── 00_normalize.py              # Standardize QTL data format
│   ├── 01_extract_sequences.py      # Extract 1MB genomic windows
│   ├── 02_predict.py                # Run AlphaGenome predictions
│   ├── 03_evaluate.py               # Compute correlation metrics
│   ├── 04_plots.py                  # Generate diagnostic plots
│   └── run_all.py                   # Master pipeline driver
│
└── results/                          # Pipeline outputs (generated)
    ├── sequences/                    # Extracted genomic sequences
    ├── predictions/                  # AlphaGenome predictions
    ├── tables/                       # Evaluation metrics
    └── plots/                        # Scatter plots, ROC curves, etc.
```

---

## 🔬 Pipeline Steps

### 1. Normalize QTL Data (`00_normalize.py`)

Convert raw QTL files to standardized format:
- Parse different QTL file formats (BED, TSV, CSV)
- Extract variant coordinates (chr, pos, ref, alt)
- Preserve effect sizes (beta) and metadata (tissue, gene)
- Output: `data/processed/variants/{dataset}_normalized.bed`

### 2. Extract Sequences (`01_extract_sequences.py`)

Generate 1MB genomic windows around variants:
- Center each variant in 1,048,576 bp window
- Extract reference sequence from GRCh38
- Add modality annotations (ATAC, DNase, RNA_SEQ, etc.)
- Output: `results/sequences/{dataset}/sequences.parquet`

### 3. Run Predictions (`02_predict.py`)

Call AlphaGenome API for variant effect predictions:
- Use `score_variant()` with recommended scorers
- Process ~16K variants (can take 3-4 hours)
- Get gene-level predictions via `tidy_scores()`
- Output: `results/predictions/{dataset}/ag_predictions.parquet`

**Fixed Nov 14, 2025:** Proper modality extraction from `output_type` field

### 4. Evaluate Predictions (`03_evaluate.py`)

Compute correlation metrics:
- **Spearman r**: Rank correlation (robust to outliers)
- **Pearson r**: Linear correlation
- **AUROC**: Classification accuracy (high vs low effects)
- **Sign concordance**: Directional accuracy (% correct)
- Output: `results/tables/{dataset}_metrics.tsv`

### 5. Generate Plots (`04_plots.py`)

Create diagnostic visualizations per modality:
- Scatter plots (predicted vs observed)
- ROC curves (classification performance)
- Precision-Recall curves
- Output: `results/plots/{dataset}/{modality}/*.png`

---

## 📈 Expected Results

Based on AlphaGenome paper (Nature, 2024):

| QTL Type | Expected r | Paper Reference |
|----------|-----------|-----------------|
| **eQTLs** | ~0.49 | Figure 4c,d (main benchmark) |
| **caQTLs** | ? | Not reported |
| **dsQTLs** | ? | Not reported |
| **hQTLs** | ? | Not reported |

**Three Possible Outcomes:**

1. ✅ **Success (r ≈ 0.49)**: Validates AlphaGenome for variant interpretation
2. ⚠️ **Partial (r = 0.1-0.3)**: Context-dependent performance
3. ❌ **Failure (r ≈ 0)**: High-value negative result (publishable)

---

## 🔧 Technical Implementation

### Scorer Mapping

| QTL Type | Modalities | AlphaGenome Scorer |
|----------|-----------|-------------------|
| eQTLs | RNA_SEQ | `GeneMaskLFCScorer(RNA_SEQ)` |
| caQTLs | ATAC, DNase, H3K27ac | `GeneMaskLFCScorer(ATAC/DNASE)`, `CenterMaskScorer(CHIP_HISTONE)` |
| dsQTLs | DNase, ATAC | `GeneMaskLFCScorer(DNASE/ATAC)` |
| hQTLs | H3K27ac, H3K4me1 | `CenterMaskScorer(CHIP_HISTONE)` |

### Critical Fixes

**November 14, 2025 - Modality Detection Bug:**
- **Problem:** All predictions labeled as "unknown" modality
- **Root Cause:** Not extracting modality from tidy_scores `output_type` column
- **Solution:** Parse `output_type` field (ATAC, DNASE, RNA_SEQ, CHIP_HISTONE)
- **Impact:** Previous results had correct predictions but wrong labels

**Previous - eQTL Gene Filtering:**
- **Problem:** AlphaGenome returns predictions for all genes in 1MB window
- **Solution:** Filter predictions to match target gene only
- **Code:** `df = df[df['gene_name'] == df['target_gene']]`

### Output Format

Each prediction contains:
```python
{
    'variant_id': 'eQTL_123',
    'chrom': 'chr1',
    'pos': 123456,
    'ref': 'A',
    'alt': 'G',
    'gene_id': 'ENSG00000123456',
    'gene_name': 'BRCA1',
    'modality': 'RNA_SEQ',
    'quantile_score': 0.73,      # Normalized (-1 to +1)
    'raw_score': 0.45,            # Log fold-change
    'beta': 0.82,                 # Observed QTL effect (ground truth)
    'tissue_ontology': 'CL:0000624',
    'biosample_name': 'CD4-positive T cell'
}
```

---

## 📊 Data Provenance

All datasets from published peer-reviewed studies:

### eQTLs & caQTLs - GSE86886
- **Study:** Gate et al. (2018)
- **Samples:** CD4+ T cells, 424 individuals
- **Assays:** RNA-seq + ATAC-seq
- **Variants:** 424 eQTLs, 3,317 caQTLs

### dsQTLs - GSE31388
- **Study:** Degner et al. (2012)
- **Samples:** CD4+ T cells, 70 individuals
- **Assay:** DNase-seq
- **Variants:** 6,070 dsQTLs

### hQTLs - GSE116193
- **Study:** Pelikan et al. (2018)
- **Samples:** CD4+ T cells, 195 individuals
- **Assays:** H3K27ac + H3K4me1 ChIP-seq
- **Variants:** 6,261 hQTLs

---

## 🎓 Scientific Value

### Why This Matters

1. **First systematic QTL benchmark** of genomic AI models
2. **Real-world validation** beyond training data
3. **Informs clinical use**: Should we trust these models for patient variants?
4. **Guides model development**: What's missing for variant prediction?

### Publication Potential

**If Positive (r ≈ 0.49):**
- "AlphaGenome Successfully Predicts Real Variant Effects"
- Target: *Nature Methods*, *Nature Biotechnology*

**If Negative (r ≈ 0):**
- "Genomic AI Models Fail to Predict Variant Effects from Sequence Alone"
- Target: *Genome Biology*, *PLoS Computational Biology*
- High-value negative result showing limitations

---

## 🐛 Known Issues & Limitations

### Fixed Issues
- ✅ Modality detection (Nov 14, 2025)
- ✅ eQTL gene filtering (Nov 13, 2025)
- ✅ Scorer deduplication (H3K27ac/H3K4me1)

### Current Limitations
- Single cell type (CD4+ T cells) - no tissue diversity
- Binary effect classification uses median threshold
- No allele frequency or linkage disequilibrium information
- API rate limits may affect runtime

---

## 📚 Usage Examples

### Quick Test on Subset

```bash
# Test on 10 variants (fast, for debugging)
python scripts/real/02_predict.py --datasets eQTLs --limit 10 --force
```

### Run Single Dataset

```bash
# Just eQTLs (fastest, ~30 min)
python scripts/real/run_all.py --datasets eQTLs
```

### Regenerate Plots Only

```bash
# If predictions already exist
python scripts/real/04_plots.py --datasets all
```

### Mock Mode (No API)

```bash
# Test pipeline without API calls
python scripts/real/run_all.py --mock --limit 100
```

---

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit pull request with clear description

---

## 📝 Citation

If you use this benchmark:

```bibtex
@software{stephenson2025alphagenome_qtl,
  author = {Stephenson, Grant S.},
  title = {AlphaGenome QTL Benchmark: Systematic Evaluation of Variant Effect Prediction},
  year = {2025},
  url = {https://github.com/gsstephenson/alphagenome-qtl-benchmark},
  note = {Layer Laboratory, University of Colorado Boulder}
}
```

Please also cite the original QTL studies and AlphaGenome paper.

---

## 📧 Contact

**Grant Stephenson**  
Layer Laboratory  
University of Colorado Boulder  
Email: grant.stephenson@colorado.edu  
GitHub: [@gsstephenson](https://github.com/gsstephenson)

---

## 📜 License

MIT License - See LICENSE file for details

---

## 🙏 Acknowledgments

- **Robin Layer** - Principal Investigator, mentorship
- **Layer Laboratory** - Computational resources
- **Google DeepMind** - AlphaGenome API access
- **Original study authors** - Public QTL data release

---

**Last Updated:** November 14, 2025  
**Current Status:** Ready for final prediction run (6 hours to completion)
