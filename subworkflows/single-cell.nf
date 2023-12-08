
include {
    JUPYTERNOTEBOOK as JUPYTER_SUBSET_ATLAS;
    JUPYTERNOTEBOOK as JUPYTER_ANNOTATE_MYELOID;
    JUPYTERNOTEBOOK as JUPYTER_TF_PW;
    JUPYTERNOTEBOOK as JUPYTER_OVERVIEW_PLOTS
} from "../modules/local/jupyternotebook/main"

def get_notebook_channel(id) {
    return Channel.value([["id": id], file("${projectDir}/analyses/20_single-cell/${id}.py")])
}

workflow W20_single_cell {

    ch_lung_atlas = Channel.fromPath("$projectDir/data/01_processed/sc_rna_seq/full_atlas_merged.h5ad")

    if (params.run_clustering) {
        JUPYTER_SUBSET_ATLAS(
            get_notebook_channel("21_subset_atlas"),
            ch_lung_atlas.map { adata -> ["atlas_adata": adata.name]},
            ch_lung_atlas
        )

        JUPYTER_ANNOTATE_MYELOID(
            get_notebook_channel("22_annotate_myeloid"),
            Channel.value([
                "adata_m": "adata_m.h5ad",
                "adata_nsclc": "adata_nsclc.h5ad",
                "hlca_markers": "hlca_cell_type_signatures.csv"
            ]),
            JUPYTER_SUBSET_ATLAS.out.artifacts.concat(
                Channel.fromPath("$projectDir/tables/gene_annotations/hlca_cell_type_signatures.csv", checkIfExists: true)
            ).collect()
        )
        ch_adata_annotated =JUPYTER_ANNOTATE_MYELOID.out.artifacts
    } else {
        ch_adata_annotated = Channel.fromPath("${ params.precomputed_adata_annotated }/*", checkIfExists: true).collect()
    }


    JUPYTER_TF_PW(
        get_notebook_channel("23_tf_pw"),
        Channel.value([
            "adata_m": "adata_myeloid_reannotated.h5ad",
            "adata_nsclc": "adata_nsclc_reannotated.h5ad"
        ]),
        ch_adata_annotated
    )
    JUPYTER_OVERVIEW_PLOTS(
        get_notebook_channel("29_overview_plots"),
        Channel.value([
            "adata_m": "adata_myeloid_reannotated.h5ad",
            "adata_nsclc": "adata_nsclc_reannotated.h5ad"
        ]),
        ch_adata_annotated
    )
}
