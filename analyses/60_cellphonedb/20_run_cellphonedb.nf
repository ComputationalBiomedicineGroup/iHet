#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process run_cellphonedb {
    tag { id }
    conda "/data/scratch/sturm/conda/envs/2021-nsclc_heterogeneity-cellphonedb"

    publishDir "../../data/60_cellphonedb/20_run_cellphonedb/${id}", mode: 'copy'
    cpus 8

    input:
    tuple val(id), path(adata), path(obs)

    output:
    path("${id}")

    script:
    """
    cellphonedb method statistical_analysis ${obs} ${adata} \\
        --iterations 100  \\
        --threads ${task.cpus} \\
        --counts-data hgnc_symbol && \\
    cellphonedb plot dot_plot && \\
    cellphonedb plot heatmap_plot ${obs}

    mv out ${id}
    """
}


workflow {
    run_cellphonedb(
        Channel.from([
            [
                "myeloid_cells",
                file("../../data/60_cellphonedb/01_prepare_input_data/adata_myeloid.h5ad"),
                file("../../data/60_cellphonedb/01_prepare_input_data/adata_myeloid.obs.csv")
            ],
            [
                "all_cells",
                file("../../data/60_cellphonedb/01_prepare_input_data/adata.h5ad"),
                file("../../data/60_cellphonedb/01_prepare_input_data/adata.obs.csv")
            ]
        ])
    )
}
