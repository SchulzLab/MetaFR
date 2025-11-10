# MetaFR
We present the MetaFR approach to learn gene-specific models that link open-chromatin variation from scATAC-seq data to gene expression from scRNA-seq. Using efficient regression trees, we illustrate that accurate expression prediction models can be learned on the single-cell or meta-cell level.  Validation was done using fine-mapped eQTLs. Meta-cell models were found to outperform single-cell models for most genes. Comparison to the SOTA method SCARlink revealed advantages of MetaFR in terms of runtime and prediction performance. MetaFR thus allows time-efficient analysis and obtains reliable models of gene expression prediction, which can be used to study gene regulation in any organism for which scRNA-seq and scATAC-seq data is available.

# Installation & Requirements

conda env create -f metafr_env.yml
conda activate metafr

# Usage

Basic command: nextflow run main.nf -c metafr.config

# Parameters

All parameters can be set in yhe configuration file (e.g. `metafr.config`) or via `--param_name ` in the command line. 
Below is a detailed list of all available parameters, grouped by module.

### MetaCellaR

| Parameter | Type | Description |
|------------|------|-------------|
| `publishDir` | `str` | Output directory for MetaCellaR results |
| `input_file` | `str` | Path to Seurat `.rds` object containing scRNA-seq and scATAC-seq data |
| `gtf_file` | `str` | GTF annotation file |
| `output_dir` | `str` | Directory for output files |
| `rna_count_slot` | `str` | Slot in the Seurat object containing RNA counts (e.g., `"assays$RNA@counts"`) |
| `atac_count_slot` | `str` | Slot in the Seurat containing ATAC counts (e.g., `"assays$ATAC@counts"`) |
| `celltype_info` | `str` | Metadata field in Seurat object containing cell type annotations |
| `assay_slot` | `str` | Metadata field in Seurat object specifying assay origin (e.g., `"meta.data$orig.ident"`) |
| `expected_cells` | `int` | Number of cells aggregated into one meta-cell (k) |
| `umap_dim` | `int` | Number of UMAP dimensions after dimension reduction |
| `reduction` | `str` | Dimensionality reduction method (default: `"umap"`) |
| `metadata` | `list[str]` | Metadata fields to guide clustering of RNA cells, e.g. cell-type or condition |
| `normalization_flag` | `bool` | Whether to normalize RNA counts |
| `umap_flag` | `bool` | Whether UMAP dimension reduction should be performed |

### MetaCell filtering

| Parameter | Type | Description |
|------------|------|-------------|
| `publishDir` | `str` | Directory for filtered meta-cell outputs |
| `atac_threshold` | `int` | Minimum number of ATAC fragments required per meta-cell |
| `out_suffix` | `str` | Suffix for filtered output file names |

### Activity Binning

| Parameter | Type | Description |
|------------|------|-------------|
| `publishDir` | `str` | Output directory for binned LSI activity files |
| `fragment_file` | `str` | Path to fragment file (`.tsv.gz`) containing ATAC fragments |
| `gtf_file` | `str` | GTF annotation file for gene coordinates |
| `chrom_sizes_file` | `str` | Chromosome size file (tab-separated) |
| `tmp_dir` | `str` | Temporary directory for intermediate processing |
| `out_dir_activity` | `str` | Directory for LSI-normalized activity matrices (one per gene) |
| `out_dir_fragments` | `str` | Directory for intermediate fragment files split by chromosome |
| `out_extension` | `str` | File extension for activity output files (e.g., `_lsi_activity_1MB500bp.txt.gz`) |
| `test_partition_flag` | `bool` | Whether to automatically generate a train/test split |
| `test_cells_file` | `str` | Path to test cell/metacell file: used as input (if `test_partition_flag`=FALSE) or output (if `test_partition_flag` =TRUE) |
| `test_cells_percent` | `int` | Percentage of cells/metacells to assign to test set (used if `test_partition_flag`=TRUE) |
| `feature_select_flag` | `bool` | Whether to perform feature selection via Seurat's `FindTopFeatures` |
| `q_top_features` | `str` | Quantile cutoff for top feature selection (e.g., `'q85'` for top 15%): only used for `feature_select_flag`=TRUE|
| `normalize_flag` | `bool` | Whether to perform CPM normalization of RNA data |
| `scaling_factor` | `int` | Scaling factor for normalization (default: `10000`): only used for `normalize_flag` =TRUE |
| `mc_flag` | `bool` | Whether to run in meta-cell mode (`TRUE`) or single-cell mode (`FALSE`) |
| `extension_size` | `int` | Genomic window extension around TSS (e.g., `500000` for ±500 kb) |
| `bin_size` | `int` | Size of genomic bins in base pairs (e.g., `500`) |
| `n_outer_cores` | `int` | Number of outer-level parallel jobs (chromosomes) |
| `n_inner_cores` | `int` | Number of inner-level parallel jobs (genes per chromosome) |
| `n_threads_feature_counts` | `int` | Threads used by Signac `FeatureMatrix()` for counting |
| `seed` | `int` | Random seed for reproducibility |

### Random Forest

| Parameter | Type | Description |
|------------|------|-------------|
| `publishDir` | `str` | Output directory for trained per gene random forest regression results |
| `input_extension` | `str` | File extension for input LSI activity files (must match `out_extension` from Activity Binning Module)
| `seed_value` | `int` | Random seed for model initialization |
| `expr_threshold` | `float` or `'None'` | Minimum variance of log2-transformed gene expression required for training |
| `min_nonzero_expr` | `int` or `'None'` | Minimum percentage of non-zero expression entries required for a gene |
| `min_mc_total` | `int` | Minimum number of metacells required to train a gene model |
| `min_nonzero_feature_mat` | `int` or `'None'` | Minimum percentage of non-zero entries in feature matrix |
| `plot_dir` | `str` | Directory for scatter plots of predicted vs. observed expression |
| `log_dir` | `str` | Directory for training logs (one log file per gene) |
| `model_dir` | `str` | Directory for serialized Random Forest model objects (`.joblib`) |
| `performance_dir` | `str` | Directory for summary files containing performance metrics |
| `n_jobs` | `int` | Number of parallel Random Forest jobs to run |
| `timeout` | `str` | Max walltime per gene process (e.g., `"10m"`, `"1h"`) |



