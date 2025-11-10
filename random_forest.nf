process random_forest {
    //label 'high_mem'
    // errorStrategy 'ignore' 
    errorStrategy 'finish'
    time params.rf_learning.timeout ?: '1h'

    publishDir "${params.rf_learning.model_dir}", mode: 'copy', pattern: "models_*", overwrite:true
    publishDir "${params.rf_learning.plot_dir}", mode: 'copy', pattern: "results_*.pdf", overwrite:true
    publishDir "${params.rf_learning.performance_dir}", mode: 'copy', pattern: "results_*", overwrite:true
    publishDir "${params.rf_learning.log_dir}", mode: 'copy', pattern: "*.log", overwrite:true
    
    input:
    //path activity_file  
    //val test_cells_file
    tuple path(activity_file), path(test_cells_file)
    
    output:
    //path "${params.rf_learning.model_dir}/models_*", emit: models
    //path "${params.rf_learning.plot_dir}/results_*.pdf", emit: plots
    //path "${params.rf_learning.performance_dir}/results_*", emit: performance
    //path "${params.rf_learning.log_dir}/*.log", emit: logs
    path "models_*", emit: models
    path "results_*.pdf", emit: plots
    path "results_*", emit: performance
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    
    set -euo pipefail
    gene_name=\$(basename "${activity_file}" | cut -d'_' -f1)
    TAG="\${gene_name}_\${RANDOM}"

    python3 ${params.rf_learning.rf_python_script} \\
        "${activity_file}" \\
        "${params.rf_learning.seed_value}" \\
        "${params.rf_learning.expr_threshold}" \\
        "${params.rf_learning.min_nonzero_expr}" \\
        "${params.rf_learning.min_mc_total}" \\
        "${params.rf_learning.min_nonzero_feature_mat}" \\
        "${test_cells_file}" \\
        "${params.rf_learning.plot_dir}" \\
        "${params.rf_learning.log_dir}/\${gene_name}.log" \\
        "${params.rf_learning.model_dir}" \\
        "${params.rf_learning.performance_dir}"

    ls -R    
    """
}

// process rf_learning{
//     input:
//         path lsi_activity_path   
//         path test_cells_file

//     output:
//         path "${params.rf_learning.model_dir}", emit: rf_results
//         path "${params.rf_learning.plot_dir}", emit: rf_plots
//         path "${params.rf_learning.log_dir}", emit: rf_logs

//     //when:
//     //   file(test_cells_file) && file(lsi_activity_path)

//     script:
//     """ 
//     ${params.rf_learning.rf_bash_script} \
//         --input_dir "${lsi_activity_path}" \
//         --test_cells "${test_cells_file}" \
//         --input_extension "${params.rf_learning.input_extension}" \
//         --model_type "${params.rf_learning.model_type}" \
//         --seed_value ${params.rf_learning.seed_value} \
//         --expr_threshold ${params.rf_learning.expr_threshold} \
//         --min_nonzero_expr ${params.rf_learning.min_nonzero_expr} \
//         --min_mc_total ${params.rf_learning.min_mc_total} \
//         --min_nonzero_feature_mat ${params.rf_learning.min_nonzero_feature_mat} \
//         --plot_dir "${params.rf_learning.plot_dir}" \
//         --log_dir "${params.rf_learning.log_dir}" \
//         --model_dir "${params.rf_learning.model_dir}" \
//         --performance_dir "${params.rf_learning.performance_dir}" \
//         --n_jobs ${params.rf_learning.n_jobs} \
//         --timeout "${params.rf_learning.timeout}"
//     """
// }

