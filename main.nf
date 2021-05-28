#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { nxfvars; file_tuple } from "./nxfvars.nf"

process P11_easier {
    publishDir "${params.result_dir}/10_mofa/${id}", mode: params.publish_dir_mode 
    // conda "/data/scratch/sturm/conda/envs/2021-nsclc_heterogeneity-easier"
    container "containers/2021-nsclc_heterogeneity-easier.sif"
    input: 
        tuple val(id), path(script)
        path("NSCLC_expr_data_sel.RData") 
        path("GTEx_expr_data.RData")
    
    output:
        path("*.rds")

    """
    Rscript ${script}
    """
}

process P12_prepare_mofa_data {
    publishDir "${params.result_dir}/10_mofa/${id}", mode: params.publish_dir_mode 
    conda "/data/scratch/sturm/conda/envs/2021-nsclc_heterogeneity-mofa"
    cpus 1
    input:
        tuple val(id), path(notebook)
        path("helper_functions.R")
        path(datasets)

    output:
        path("*.html"), emit: report
        path("datasets/*.rds"), emit: datasets 

    script:
    data_dir = "./"
    out_dir = "datasets"
    """
    ${nxfvars(task)}
    execute_rmd.r ${notebook}
    """
}


workflow W10_mofa {
    def dir = "analyses/10_mofa"
    P11_easier(
        file_tuple("$dir/11_easier.R"),
        file("data/01_processed/bulk_rna_seq/NSCLC_expr_data_sel.RData"),
        file("data/01_processed/bulk_rna_seq/GTEx_expr_data.RData")
    )
    P12_prepare_mofa_data(
        file_tuple("$dir/12_prepare_mofa_data.Rmd"),
        file("analyses/10_mofa/helper_functions.R"),
        P11_easier.out
    )
}

workflow {
    W10_mofa()
}
