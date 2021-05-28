#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { nxfvars } from "./nxfvars.nf"

process P11_easier {
    def id = "11_easier" 
    publishDir "${params.result_dir}/10_mofa/${id}", mode: params.publish_dir_mode 
    // conda "/data/scratch/sturm/conda/envs/2021-nsclc_heterogeneity-easier"
    container "containers/2021-nsclc_heterogeneity-easier.sif"
    input: 
        path("${id}.R")
        path("NSCLC_expr_data_sel.RData") 
        path("GTEx_expr_data.RData")
    
    output:
        path("*.rds")

    """
    Rscript ${id}.R
    """
}

process P12_mofa_analysis {
    def id = "12_mofa_analysis"
    publishDir "${params.result_dir}/10_mofa/${id}", mode: params.publish_dir_mode 
    conda "/data/scratch/sturm/conda/envs/2021-nsclc_heterogeneity-mofa"
    cpus 10
    input:
        path("${id}.Rmd")
        path("helper_functions.R")
        path("*.rds")

    output:
        path("${id}.html"), emit: report
        path("models/*.h5ad"), emit: models

    script:
    def data_dir = "."
    def model_dir = "models"
    """
    ${nxfvars(task)}
    execute_rmd.r ${id}.Rmd ${id}.html
    """
}

workflow W10_mofa {
  P11_easier(
      file("analyses/10_mofa/11_easier.R"),
      file("data/01_processed/bulk_rna_seq/NSCLC_expr_data_sel.RData"),
      file("data/01_processed/bulk_rna_seq/GTEx_expr_data.RData")
  )
  P12_mofa_analysis(
      file("analyses/10_mofa/12_mofa_analysis.Rmd"),
      file("analyses/10_mofa/helper_functions.R"),
      P11_easier.out
  )
}

workflow {
    W10_mofa()
}
