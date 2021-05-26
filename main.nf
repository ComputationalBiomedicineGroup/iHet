#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
        "*.rds"

    """
    Rscript ${id}.R
    """
}

// process P12_mofa_analysis {
//     def id = "12_mofa_analysis"
//     conda "/data/scratch/sturm/conda/envs/2021-nsclc_heterogeneity-mofa"
//     input:
//         path("${id}.Rmd")

//     output:
//         "${id}.html", emit: report

//     script:
//     """
//     """

// }

workflow W10_mofa {
  P11_easier(
      file("analyses/10_mofa/11_easier.R"),
      file("data/01_processed/bulk_rna_seq/NSCLC_expr_data_sel.RData"),
      file("data/01_processed/bulk_rna_seq/GTEx_expr_data.RData")
  )
}

workflow {
    W10_mofa()
}
