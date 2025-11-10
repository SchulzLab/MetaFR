process activity_binning {

    //errorStrategy 'ignore' //TODO: remove this line when debugging is done
    errorStrategy 'finish'
    publishDir "${params.activity_binning.publishDir}", mode: 'copy'

    input:
    path cb_mc_file
    path rna_file
   
    output:
    path "${params.activity_binning.out_dir_activity}/*${params.activity_binning.out_extension}", emit: lsi_activity_files //wait until all files are generated?
    path "${params.activity_binning.out_dir_activity}", emit: lsi_activity_path //needed as input for random forest
    //path "${params.activity_binning.test_cells_file}", emit: test_cells_file ///needed as input for random forest
    path "test_*.txt", optional: true, emit: test_cells_file ///needed as input for random forest

    script:
    """
    #!/usr/bin/env Rscript
    source("${params.activity_binning.activity_script}")
 
    activity.bin(
        cb_mc_file = "${cb_mc_file}",
        rna_file = "${rna_file}",
        fragment_file = "${params.activity_binning.fragment_file}", 
        gtf_file = "${params.gtf_file}", 
        chrom_sizes_file = "${params.activity_binning.chrom_sizes_file}",
        tmp_dir = "${params.activity_binning.tmp_dir}",
        out_dir_activity = "${params.activity_binning.out_dir_activity}",
        out_dir_fragments = "${params.activity_binning.out_dir_fragments}",
        out_extension = "${params.activity_binning.out_extension}",
        test_partition_flag = ${params.activity_binning.test_partition_flag},
        test_cells_file = "${params.activity_binning.test_cells_file}",
        test_cells_percent = ${params.activity_binning.test_cells_percent},
        feature_select_flag = ${params.activity_binning.feature_select_flag},
        q_top_features = "${params.activity_binning.q_top_features}",
        normalize_flag = ${params.activity_binning.normalize_flag},
        scaling_factor = ${params.activity_binning.scaling_factor},
        mc_flag = ${params.activity_binning.mc_flag},
        extension_size = ${params.activity_binning.extension_size},
        bin_size = ${params.activity_binning.bin_size},
        n_outer_cores = ${params.activity_binning.n_outer_cores},
        n_inner_cores = ${params.activity_binning.n_inner_cores},
        n_threads_feature_counts = ${params.activity_binning.n_threads_feature_counts},
        seed = ${params.activity_binning.seed})
    """
}

 
