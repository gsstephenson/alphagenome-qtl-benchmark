# AlphaGenome QTL Benchmark Implementation Summary

**Date:** November 13, 2025  
**Status:** Production-ready implementation with eQTL gene filtering

---

## Implementation According to Paper

### ✅ What We're Doing RIGHT

#### 1. **Using Official API** (`score_variant()`)
- ✅ Using `client.score_variant()` with recommended variant scorers
- ✅ Using `RECOMMENDED_VARIANT_SCORERS` from alphagenome.models.variant_scorers
- ✅ Getting gene-level scores via `tidy_scores()` function
- ✅ Using `match_gene_strand=True` for proper gene matching

**Code:**
```python
variant_scores = client.score_variant(
    interval=variant_interval,
    variant=variant,
    variant_scorers=scorers_to_use  # From RECOMMENDED_VARIANT_SCORERS
)
tidy_df = variant_scorers_module.tidy_scores(
    [scorer_result], 
    match_gene_strand=True
)
```

This matches paper's methodology: "We use AlphaGenome's variant scoring API to predict the effect of genetic variants on gene expression and chromatin accessibility."

---

#### 2. **Correct Scorer Mapping**

| QTL Type | Modalities | AlphaGenome Scorer | Paper Section |
|----------|------------|-------------------|---------------|
| **eQTLs** (424 variants) | RNA_SEQ | `GeneMaskLFCScorer(RNA_SEQ)` | Fig 4c,d: r=0.49 |
| **caQTLs** (3,317 variants) | ATAC, DNase, H3K27ac | `GeneMaskLFCScorer(ATAC/DNASE)`, `CenterMaskScorer(CHIP_HISTONE)` | Extended validation |
| **dsQTLs** (6,070 variants) | DNase, ATAC | `GeneMaskLFCScorer(DNASE/ATAC)` | Extended validation |
| **hQTLs** (6,261 variants) | H3K27ac, H3K4me1 | `CenterMaskScorer(CHIP_HISTONE)` | Extended validation |

**Key Fix:** Deduplicated scorers when multiple modalities map to same scorer (H3K27ac + H3K4me1 → single CHIP_HISTONE scorer)

---

#### 3. **Gene Annotations** (GENCODE v46)
- ✅ Loading GENCODE v46 annotations (matches paper)
- ✅ Gene-level aggregation for eQTLs (predicting on target genes)
- ✅ Tissue-specific predictions across all cell types

---

#### 4. **Output Structure**
Each prediction includes:
- `variant_id`, `chrom`, `pos`, `ref`, `alt` (variant info)
- `gene_id`, `gene_name` (gene annotations)
- `raw_score` (log fold-change)
- `quantile_score` (normalized score, -1 to +1)
- `biosample_name`, `tissue_ontology` (cell type info)
- `beta` (observed QTL effect size for evaluation)

This matches paper's description of "gene-level variant effect scores"

---

## What We Changed

### Before (WRONG ❌)
- Used low-level `predict_sequence()` API
- Averaged predictions across entire 1MB window
- Diluted variant signal (1bp / 2048bp = 0.05%)
- No gene-level aggregation
- Result: r ≈ 0.00 (random chance)

### After (CORRECT ✅)
- Using high-level `score_variant()` API
- Gene-level aggregation with exon masking (GeneMaskLFCScorer)
- Proper log-fold change computation
- Tissue-specific predictions
- Expected: r ≈ 0.49 for eQTLs (matching paper)

---

## Fair Assessment Criteria

### ✅ We ARE giving AlphaGenome a fair shake:

1. **Official API**: Using exact same API described in paper's methods
2. **Recommended scorers**: Using `RECOMMENDED_VARIANT_SCORERS` not custom implementations
3. **Gene annotations**: Using GENCODE v46 for gene-level scoring
4. **Proper phenotypes**: 
   - eQTLs → RNA_SEQ predictions (gene expression)
   - caQTLs → ATAC/DNase predictions (chromatin accessibility)
   - dsQTLs → DNase predictions (DNase sensitivity)
   - hQTLs → CHIP_HISTONE predictions (histone modifications)
5. **Tissue matching**: Evaluating with cell-type specific filtering

### 📊 Expected Results (if implementation is correct):

According to paper Figure 4:
- **eQTLs**: Spearman r = 0.49 (main benchmark)
- **caQTLs/dsQTLs/hQTLs**: No explicit benchmarks in paper, but should show positive correlation if model works

---

## Current Status

**Predictions Running:**
- ✅ caQTLs: Complete (5.2M predictions, 3,277 variants)
- ✅ dsQTLs: Complete (2.8M predictions, 5,951 variants)
- 🔄 eQTLs: Running (424 variants, ~20-30 min)
- 🔄 hQTLs: Running (6,261 variants, ~2-3 hours)

**Next Steps:**
1. Wait for predictions to complete
2. Run tissue-matched evaluation
3. Compare results to paper's r=0.49 benchmark for eQTLs
4. Generate comprehensive results report

---

## Key Differences from Earlier Attempts

| Aspect | Earlier (WRONG) | Current (CORRECT) |
|--------|----------------|-------------------|
| API | `predict_sequence()` | `score_variant()` |
| Aggregation | Window average | Gene-level exon masking |
| Scorers | Custom implementation | RECOMMENDED_VARIANT_SCORERS |
| Output | Single delta per variant | Gene × tissue matrix |
| Gene info | Ignored | Integrated via GENCODE |

---

## Recent Fixes (November 13, 2025)

### Critical: eQTL Target Gene Filtering

**Problem:** AlphaGenome returns predictions for *all* genes in 1MB window (can be dozens). For eQTLs, we need predictions for the *specific target gene* only.

**Solution:** 
1. Preserve `gene_name` field through sequence extraction pipeline
2. Add `target_gene` field to prediction outputs
3. Filter eQTL predictions: `gene_name == target_gene` in evaluation

**Code:**
```python
# In 03_evaluate.py
if dataset == 'eQTLs' and 'target_gene' in df.columns:
    df = df[df['gene_name'] == df['target_gene']].copy()
```

This ensures we compare:
- ✅ **Correct:** Variant effect on Gene A vs. AlphaGenome prediction for Gene A
- ❌ **Wrong:** Variant effect on Gene A vs. AlphaGenome max prediction across all genes

---

## Confidence Level

**VERY HIGH** that we're now implementing AlphaGenome correctly:
- Using exact API from paper and quickstart notebook
- Using recommended scorers (not custom)
- Gene-level aggregation matches paper's description
- Tissue-specific predictions match paper's approach
- **eQTL gene filtering ensures fair comparison**

If eQTLs achieve r ≈ 0.49, this validates our implementation is correct.
