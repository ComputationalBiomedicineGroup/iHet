#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process run_cellphonedb {
    tag { id }
    conda "/data/scratch/sturm/conda/envs/2021-nsclc_heterogeneity-cellphonedb"

    publishDir "../../data/60_cellphonedb/20_run_cellphonedb/${id}", mode: 'copy'
    cpus 22

    input:
    tuple val(id), path(adata), path(obs)

    output:
    path("out/*")

    script:
    """
    cellphonedb method statistical_analysis ${obs} ${adata} \\
        --iterations 1000  \\
        --threads ${task.cpus} \\
        --counts-data hgnc_symbol && \\
    cellphonedb plot dot_plot && \\
    cellphonedb plot heatmap_plot ${obs}
    """
}


workflow {
    run_cellphonedb(
        Channel.from([
            [
                "all_cells",
                file("../../data/60_cellphonedb/01_prepare_input_data/adata.h5ad"),
                file("../../data/60_cellphonedb/01_prepare_input_data/adata.obs.csv")
            ]
        ])
    )
}
