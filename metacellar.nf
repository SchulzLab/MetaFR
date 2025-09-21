process metacellar {
    
    publishDir "${params.metacellar.publishDir}", mode: 'copy'
    
    input:
    path input_file
    path gtf_file
    path metacellar_script
    
    output:
    //path "${params.metacellar.output_dir}/results/cellSummarized.csv"
    path "${params.metacellar.output_dir}/results/cellSummarized_sum.csv", emit: rna_counts 
    //path "${params.metacellar.output_dir}/results/cellSummarized_ATAC.csv"
    path "${params.metacellar.output_dir}/results/cellSummarized_ATAC_sum.csv", emit: atac_counts
    //path "${params.metacellar.output_dir}/results/cellSummarized_normalized.csv", optional: true
    path "${params.metacellar.output_dir}/results/RNA_cellSummarized_normalized.csv", emit: rna_data, optional: true
    path "${params.metacellar.output_dir}/results/RNA_cell2metacell_info.csv", emit: mc_rna_info
    path "${params.metacellar.output_dir}/results/ATAC_cell2metacell_info.csv", emit: mc_atac_info
    path "${params.metacellar.output_dir}/results/RNA_metacell_umap.csv"
    path "${params.metacellar.output_dir}/plots/*.pdf"
    
    script:
    def metadata_str = params.metacellar.metadata.collect { "'$it'" }.join(', ')
    
    """
    #!/usr/bin/env Rscript
    source("${metacellar_script}")
    
    metacellar.run(
        file_name = "${input_file}",
        RNA_count_slot = "${params.metacellar.rna_count_slot}",
        celltype_info = "${params.metacellar.celltype_info}",
        output_file = "${params.metacellar.output_dir}",
        umap_flag = ${params.metacellar.umap_flag},
        assay_slot = "${params.metacellar.assay_slot}",
        ATAC_count_slot = "${params.metacellar.atac_count_slot}",
        expected_cells = ${params.metacellar.expected_cells},
        umap_dim = ${params.metacellar.umap_dim},
        reduction = "${params.metacellar.reduction}",
        gtf_file = "${gtf_file}",
        metadata = c(${metadata_str}),
        normalization_flag = ${params.metacellar.normalization_flag}
    )
    """
}
