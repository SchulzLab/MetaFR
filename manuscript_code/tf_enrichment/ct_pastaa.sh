#!/usr/bin/env bash

INPUT_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/"
OUTPUT_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pastaa_results/pastaa_out_ct"
#AFFINITY_FILE="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pastaa_results/out/out_Affinity.txt"
#AFFINITY_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/affinity_subsets"
#AFFINITY_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/ct_bed_regions_affinity/"
AFFINITY_DIR="/projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pastaa_results/trap_out_ct/"

PASTAA_BIN="/projects/single_cell_stitchit/work/tf_enrichment/PASTAA/PASTAA"  

for ct_dir in "${AFFINITY_DIR}"/*
do
    [ -d "$ct_dir" ] || continue
    #filename=$(basename "$affinity_file")
    #ct=${filename#affinity_subset_}
    #ct=${ct%.txt}
    ct=$(basename "$ct_dir")
    affinity_file="${ct_dir}"/"${ct}"_Affinity_1000.txt  
    echo "$affinity_file"
    #echo "${INPUT_DIR}"/pbmc_mc_500bp_top_10000_abs_shap_zscores_${ct}_regions.bed
    echo "${AFFINITY_DIR}"/"${ct}"_ranked_regions_1000.txt
    echo "run PASTAA for ${ct}"
    #"$PASTAA_BIN" "$affinity_file" "${INPUT_DIR}/pbmc_mc_500bp_top_10000_abs_shap_zscores_${ct}_regions.bed" #|  sort -k2,2 -g  > "${OUTPUT_DIR}/pbmc_mc_500bp_top_10000_abs_shap_zscores_${ct}_pastaa_out.txt"
    #"$PASTAA_BIN" "$affinity_file" "${AFFINITY_DIR}"/"${ct}"_affinity_regions_1000.bed |  sort -k2,2 -g  #> "${OUTPUT_DIR}/pbmc_mc_500bp_top_1000_abs_shap_zscores_${ct}_pastaa_out.txt"
    "$PASTAA_BIN" "$affinity_file" "${AFFINITY_DIR}"/"${ct}"_ranked_regions_1000.txt |  sort -k2,2 -g  > "${OUTPUT_DIR}/pbmc_mc_500bp_top_1000_abs_shap_zscores_${ct}_pastaa_out.txt"
    echo "finished ${ct}"
done

#./PASTAA /projects/single_cell_stitchit/work/tf_enrichment/affinity_subset_HSPC.txt /projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pbmc_mc_500bp_top_1000_abs_shap_zscores_HSPC_regions.bed

#./PASTAA  /projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pastaa_results/out/out_Affinity.txt   /projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pbmc_mc_500bp_top_1000_abs_shap_zscores_each_ct_train_ct_marker_genes_lfc_025_regions.bed  |  sort -k2,2 -g  >  /projects/single_cell_stitchit/work/metafr_revision_2025/ct_shap_analysis/pastaa_results/out/pastaa_results_pbmc_mc_500bp_top_1000_abs_shap_zscores_each_ct_train_ct_marker_genes_lfc_025_genes.tx
