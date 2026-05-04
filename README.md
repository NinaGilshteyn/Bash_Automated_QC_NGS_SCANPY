# Bash-Automated SCANPY NGS Quality Control

**Author:** Nina Gilshteyn, M.S.
**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
**Status:** Manuscript in preparation — "The Biology is in The Tails: Skewness and Upregulation in scRNA-seq Data"
Based on: Gilshteyn, N. (2024). *The Biology is in The Tails: Skewness and Upregulation in scRNA-seq*. UCLA. ProQuest ID: Gilshteyn_ucla_0031N_23577. Merritt ID: ark:/13030/m5b96p3g. Retrieved from https://escholarship.org/uc/item/8p46s4ck
**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.com/NinaGilshteyn](https://github.com/NinaGilshteyn)

---

## Overview

This pipeline performs automated quality control (QC) and cell filtration for single-cell RNA-sequencing (scRNA-seq) data using [Scanpy](https://scanpy.readthedocs.io/). It operates in two stages: first generating QC visualizations to guide threshold selection, then applying sequential multi-step filtration to remove low-quality cells based on gene count, transcript count, and mitochondrial content.

The pipeline accepts raw gene expression matrices in CSV format and produces publication-ready QC figures alongside clean, filtered expression matrices ready for downstream analysis.

---

## Repository Structure

```
Bash_Automated_SCANPY_NGS_QC/
├── NEW_get_QC_plots.py              # STEP 1: Compute and visualize QC metrics
├── Data_filtering_run_this.py       # STEP 2: Apply sequential cell filtration
└── modules/
    ├── NEW_QC_mod.py                # Scanpy-based QC metric computation and 3D/2D plots
    ├── filter_mod.py                # Cell filtering logic with multiple filter strategies
    └── QCpltGenerator_mod_for_get_QC_plots.py  # Legacy QC scatter/histogram plotting
```

---

## Pipeline

### STEP 1 — QC Visualization (`NEW_get_QC_plots.py`)

Reads a raw gene expression CSV (genes × cells), transposes to cells × genes format, computes QC metrics via Scanpy, and generates histograms and scatter plots to guide threshold selection.

**Run:**
```bash
python NEW_get_QC_plots.py <path_to_expression_csv>
```

**What it does:**
1. Reads raw gene expression matrix (genes as rows, cells as columns)
2. Identifies mitochondrial genes by `mt-` prefix
3. Calculates per-cell QC metrics using Scanpy:
   - `n_genes_by_counts` — number of genes with non-zero counts
   - `total_counts` — total UMI/transcript count
   - `pct_counts_mt` — percentage of counts from mitochondrial genes
4. Generates log-scale and linear-scale histograms for each metric
5. Creates an interactive 3D Plotly scatter plot and a 2D scatter colored by mitochondrial percentage
6. Saves QC statistics to CSV for use in STEP 2

**Outputs:**
| File | Description |
|------|-------------|
| `QC_stats_output_[name].csv` | Per-cell QC metrics (n_genes, total_counts, pct_mt) |
| `[name]_n_genes_by_counts_QC_hist_log.png` | Log-scale histogram of gene counts |
| `[name]_total_counts_QC_hist_log.png` | Log-scale histogram of transcript counts |
| `[name]_pct_counts_mt_QC_hist_log.png` | Log-scale histogram of mitochondrial % |
| `QC_3D_plot_[name].html` | Interactive 3D QC scatter (Plotly) |
| `QC_3D_plot_gene_[name].png` | Static 3D QC scatter |
| `QC_2D_scatter_Colored_by_mito_percent_[name].png` | 2D scatter colored by mito % |

---

### STEP 2 — Data Filtering (`Data_filtering_run_this.py`)

Applies three sequential exclusive filters to remove low-quality cells. Filtering thresholds are set per dataset inside the script based on the QC distributions observed in STEP 1.

**Run:**
```bash
python Data_filtering_run_this.py <path_to_expression_csv>
```

**What it does:**
1. Reads gene expression matrix and QC statistics CSV (output of STEP 1)
2. Applies three sequential filtration steps using `filter_mod.filtercells()`:
   - **Filter 1 — Gene count:** Remove cells with gene counts outside specified bounds
   - **Filter 2 — Transcript count:** Remove cells with total counts outside specified bounds
   - **Filter 3 — Mitochondrial content:** Remove cells exceeding the mitochondrial percentage threshold
3. Generates post-QC histograms and 3D/2D QC plots for the filtered dataset
4. Saves three output files

**Outputs:**
| File | Description |
|------|-------------|
| `POST_QC_purified_gex_with_QC_stats_[name].csv` | Filtered expression matrix with QC columns |
| `POST_QC_QC_stats_from_Purified_gex[name].csv` | QC metrics for filtered cells only |
| `POST_QC_purified_gex_in_standard_format_[name].csv` | Filtered expression matrix without QC columns (ready for downstream analysis) |
| `[step]_removed_cells_*.csv` | Cells removed at each filtration step |
| Post-QC histograms and scatter plots | Same format as STEP 1 outputs |

---

## Input Format

- **CSV file** — gene expression matrix
- **Orientation for STEP 1:** genes as rows, cells as columns
- **QC columns required for STEP 2:** `n_genes_by_counts`, `total_counts`, `pct_counts_mt` (from STEP 1 output)
- **Mitochondrial gene naming:** genes prefixed with `mt-` are recognized as mitochondrial

---

## Setting Filtering Thresholds

After running STEP 1, inspect the histograms to choose appropriate cutoffs for your dataset. Edit the following variables in `Data_filtering_run_this.py`:

```python
gene_cutoff    = 4500    # upper bound on gene counts per cell
count_cutoff   = 20000   # upper bound on transcript counts per cell
mito_cutoff    = 30      # maximum mitochondrial percentage (%)
```

Thresholds are dataset-specific; values shown above are from the K6 dataset used in the associated manuscript.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `scanpy` | QC metric computation and AnnData objects |
| `pandas` | Data I/O and manipulation |
| `numpy` | Numerical operations |
| `matplotlib` | Static plotting |
| `seaborn` | Enhanced statistical plots |
| `plotly` | Interactive 3D/2D scatter plots |

**Install with conda:**
```bash
conda install scanpy pandas numpy matplotlib seaborn plotly
```

**Or with pip:**
```bash
pip install scanpy pandas numpy matplotlib seaborn plotly
```

---

## Module Descriptions

### `modules/NEW_QC_mod.py`
Core QC analytics module. Wraps Scanpy to calculate per-cell QC metrics, generates log- and linear-scale histograms for each metric, and produces interactive 3D and 2D Plotly visualizations.

### `modules/filter_mod.py`
Implements `filtercells()` — a flexible cell filtering function supporting multiple strategies (exclusive, inclusive, or exact-match bounds) for gene counts, transcript counts, and mitochondrial percentage. Produces filtered dataframes, removed-cell logs, and per-step histograms.

### `modules/QCpltGenerator_mod_for_get_QC_plots.py`
Legacy QC module retained for backward compatibility and cross-validation. Generates scatter plots (genes vs. transcripts) and histograms from raw expression data without Scanpy.

---

## Publication

If you use this pipeline, please cite:

> Gilshteyn, N. (2024). *The Biology is in The Tails: Skewness and Upregulation in scRNA-seq*. UCLA. ProQuest ID: Gilshteyn_ucla_0031N_23577. Merritt ID: ark:/13030/m5b96p3g. Retrieved from https://escholarship.org/uc/item/8p46s4ck

---

## License

Free for academic and non-commercial use. Commercial use requires a written agreement. See [CommercialLicense](CommercialLicense) for full terms or contact [ninagilshteyn@gmail.com](mailto:ninagilshteyn@gmail.com).
