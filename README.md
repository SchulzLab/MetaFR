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



