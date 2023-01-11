
include { JUPYTERNOTEBOOK as JUPYTER_SUBSET_ATLAS } from "../modules/local/jupyternotebook/main"

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
}
