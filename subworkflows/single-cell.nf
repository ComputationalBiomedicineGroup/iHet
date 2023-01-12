
include { JUPYTERNOTEBOOK as JUPYTER_SUBSET_ATLAS } from "../modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_ANNOTATE_MYELOID } from "../modules/local/jupyternotebook/main"

def get_notebook_channel(id) {
    return Channel.value([["id": id], file("${projectDir}/analyses/20_single-cell/${id}.py")])
}

workflow W20_single_cell {

    ch_lung_atlas = Channel.fromPath("$projectDir/data/01_processed/sc_rna_seq/full_atlas_merged.h5ad")

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
}
