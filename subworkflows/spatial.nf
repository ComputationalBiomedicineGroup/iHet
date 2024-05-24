include {
    JUPYTERNOTEBOOK as JUPYTER_SPATIAL
} from "../modules/local/jupyternotebook/main"

def get_notebook_channel(id) {
    return Channel.value([["id": id], file("${projectDir}/analyses/30_spatial/${id}.py")])
}

workflow W30_spatial {
    visium_slide = file("${projectDir}/data/01_processed/spatial/lungCancerVisium/outs", checkIfExists: true)
    cell2location_estimates = file("${projectDir}/data/01_processed/spatial/cell2location_estimates.csv", checkIfExists: true)
    JUPYTER_SPATIAL(
        get_notebook_channel("31_make_plots"),
        Channel.value([
            "visium_slide": visium_slide.name,
            "cell2location": cell2location_estimates.name
        ]),
        Channel.from(visium_slide, cell2location_estimates).collect()
    )
}
