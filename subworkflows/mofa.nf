
nextflow.enable.dsl = 2
def modules = params.modules.clone()

def file_tuple(filename) {
    f = file(filename, checkIfExists: true)
    return tuple([id: f.simpleName], f)
}

include { RMARKDOWNNOTEBOOK as P12_prepare_mofa_data } from "../modules/local/rmarkdownnotebook/main.nf" addParams(
    options: modules["P12_prepare_mofa_data"]
)
include { RMARKDOWNNOTEBOOK as P14_mofa_analysis } from "../modules/local/rmarkdownnotebook/main.nf" addParams(
    options: modules["P14_mofa_analysis"]
)

process P11_easier {
    input:
        tuple val(meta), path(script)
        path(expr_data)

    output:
        path("*.features.rds")

    script:
    """
    Rscript ${script} \\
        ${expr_data} \\
        ${baseDir}/tables/easier_Zscores/immscoreZ.rds \\
        ${expr_data.baseName}.features.rds
    """
}

process P13_run_mofa {
    input:
        tuple val(meta), path(script)
        each path(dataset)

    output:
        path("*.hdf5"), emit: models

    script:
    """
    Rscript ${script} ${dataset}
    """
}

workflow W10_mofa {
    def dir = "analyses/10_mofa"
    P11_easier(
        file_tuple("$dir/11_easier.R"),
        Channel.fromPath("${ params.bulk_input_dir }/*_expr_data*.rds")
    )
    P12_prepare_mofa_data(
        file_tuple("$dir/12_prepare_mofa_data.Rmd"),
        [:],
        P11_easier.out.collect()
    )
    P13_run_mofa(
        file_tuple("$dir/13_run_mofa.R"),
        P12_prepare_mofa_data.out.artifacts.map{ meta, it -> it }.flatten().filter{
            it -> it.simpleName.startsWith("mofa_")
        }
    )
    P14_mofa_analysis(
        file_tuple("$dir/14_mofa_analysis.Rmd"),
        [features_dir: "./", easier_dir: "./", tmb_dir: "./"],
        P12_prepare_mofa_data.out.artifacts.map{ meta, it -> it}.mix(
            Channel.of(file("$dir/helper_functions.R")),
            P11_easier.out.collect(),
            Channel.fromPath("${ params.bulk_input_dir }/*.{txt,tsv}").collect(),
            P13_run_mofa.out.models.collect()
        ).collect()
    )
}
