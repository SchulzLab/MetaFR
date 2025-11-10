nextflow.enable.dsl = 2

include { metacellar } from './metacellar.nf'
include { metacell_filtering } from './metacell_filtering.nf'
include { activity_binning } from './activity_binning.nf'
include { random_forest } from './random_forest.nf'

workflow {

    seurat_ch = channel.fromPath(params.metacellar.input_file, checkIfExists: true)
    gtf_ch = channel.fromPath(params.metacellar.gtf_file, checkIfExists: true)
    rscript_ch = channel.fromPath(params.metacellar.metacellar_script, checkIfExists: true)

    metacellar_out = metacellar(
        seurat_ch,
        gtf_ch,
        rscript_ch
    )

    //metacellar_out.mc_rna_info.view()
    //metacellar_out.mc_atac_info.view()
    //metacellar_out.rna_counts.view()
    
    def normalize = params.metacellar.normalization_flag == 'T' || params.metacellar.normalization_flag == true

    if (normalize) { // if normalization is enabled, use normalized RNA counts
        metacell_filtering_out = metacell_filtering(
        metacellar_out.rna_data,
        metacellar_out.atac_counts, 
        metacellar_out.mc_atac_info 
       //params.metacell_filtering.atac_threshold,
       // params.metacell_filtering.out_suffix      
    )
        metacell_filtering_out.rna_data.view()
        metacell_filtering_out.atac_cell2metacell.view()

        activity_binning_out = activity_binning(
            metacell_filtering_out.atac_cell2metacell,
            metacell_filtering_out.rna_data
            
        )

        //random_forest(
        //    activity_binning_out.lsi_activity_files,
        //    activity_binning_out.test_cells_file
        //)

    } else {
        metacell_filtering_out = metacell_filtering(
            metacellar_out.rna_counts,
            metacellar_out.atac_counts, 
            metacellar_out.mc_atac_info 
            //params.metacell_filtering.atac_threshold,
            //params.metacell_filtering.out_suffix
        )
        metacell_filtering_out.rna_counts.view()
        metacell_filtering_out.atac_cell2metacell.view()

        activity_binning_out = activity_binning(
            metacell_filtering_out.atac_cell2metacell,
            metacell_filtering_out.rna_counts
            
        )

        //activity_binning_out.lsi_activity_files.view()
        //activity_binning_out.lsi_activity_path.view()
        //activity_binning_out.test_cells_file.view()

        //def test_cells_val = activity_binning_out.test_cells_file.first() //get the single value from the channel

        //random_forest(
        //   activity_binning_out.lsi_activity_files,
        //    test_cells_val
        //)
       
    }

    test_cells_ch = Channel.empty()
    if (params.activity_binning.test_partition_flag.toString().toLowerCase() == 'true') {
        test_cells_ch = activity_binning_out.test_cells_file
    } else {
        test_cells_ch = Channel.fromPath(params.activity_binning.test_cells_file, checkIfExists: true)
    }

        activity_binning_out.lsi_activity_files.view()
        activity_binning_out.lsi_activity_path.view()
        test_cells_ch.view()


    activity_binning_out.lsi_activity_files
        .combine(test_cells_ch)
        | random_forest
    
}

