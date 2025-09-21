# MetaFR
We present the MetaFR approach to learn gene-specific models that link open-chromatin variation from scATAC-seq data to gene expression from scRNA-seq. Using efficient regression trees, we illustrate that accurate expression prediction models can be learned on the single-cell or meta-cell level.  Validation was done using fine-mapped eQTLs. Meta-cell models were found to outperform single-cell models for most genes. Comparison to the SOTA method SCARlink revealed advantages of MetaFR in terms of runtime and prediction performance. MetaFR thus allows time-efficient analysis and obtains reliable models of gene expression prediction, which can be used to study gene regulation in any organism for which scRNA-seq and scATAC-seq data is available.

# Installation & Requirements

conda env create -f metafr_env.yml
conda activate metafr

# Usage

nextflow run main.nf -c metafr.config
