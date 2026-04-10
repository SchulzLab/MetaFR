#!/usr/bin/env bash

#REGION_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/"
REGION_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/ct_bed_regions_affinity"
OUTPUT_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pastaa_results/trap_out_ct"

for bed_file in "${REGION_DIR}"/*_regions_10000.bed
do
    filename=$(basename "$bed_file")
    #ct=${filename#pbmc_mc_500bp_top_10000_abs_shap_zscores_}
    ct=${filename%_regions_10000.bed}
    echo ${bed_file}
    echo "run TRAP for ${ct}"
    /projects/triangulate/work/STARE/Legacy_Code/STARE_tillTRAP.sh -g /projects/abcpp/work/base_data/hg38.fa -s /projects/triangulate/work/STARE/PWMs/2.2/Jaspar_Hocomoco_Kellis_human_transfac.txt -b "$bed_file" -o "${OUTPUT_DIR}"/"${ct}" 

done

#/projects/triangulate/work/STARE/Legacy_Code/STARE_tillTRAP.sh -g /projects/abcpp/work/base_data/hg38.fa -s /projects/triangulate/work/STARE/PWMs/2.2/Jaspar_Hocomoco_Kellis_human_transfac.txt -b /projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pbmc_mc_500bp_top_1000_abs_shap_zscores_each_ct_train_ct_marker_genes_lfc_025_regions.bed  -o /projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pastaa_results/out



