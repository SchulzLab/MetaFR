process metacell_filtering {

    publishDir "${params.metacell_filtering.publishDir}", mode: 'copy'

    input:
    path rna_counts_file // can be normalized or raw RNA counts
    path atac_counts_file
    path atac_cell2metacell_file
    //val threshold
    //val out_suffix

    output:
    path "cellSummarized_sum_RNA_valid_metacells.csv", emit: rna_counts, optional: true
    path "RNA_cellSummarized_normalized_RNA_valid_metacells.csv", emit: rna_data, optional: true
    path "ATAC_cell2metacell_info_valid_metacells.csv", emit: atac_cell2metacell

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.metacell_filtering.filtering_script}")

    filter_metacells(
        rna.counts.file = "${rna_counts_file}",
        atac.counts.file = "${atac_counts_file}",
        atac.cell2metacell.file = "${atac_cell2metacell_file}",
        atac.threshold = ${params.metacell_filtering.atac_threshold},
        output.suffix = "${params.metacell_filtering.out_suffix}"
    )
    """
}