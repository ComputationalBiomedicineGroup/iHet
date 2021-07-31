
nextflow.enable.dsl = 2
def modules = params.modules.clone()

def file_tuple(filename) {
    f = file(filename)
    return tuple([id: f.simpleName], f)
}

include { RMARKDOWN as P12_prepare_mofa_data } from "../modules/local/rmarkdown/main.nf" addParams(
    options: modules["P12_prepare_mofa_data"]
)
include { RMARKDOWN as P14_mofa_analysis } from "../modules/local/rmarkdown/main.nf" addParams(
    options: modules["P14_mofa_analysis"]
)

process P11_easier {
    publishDir "${params.outdir}/10_mofa/${meta.id}", mode: params.publish_dir_mode

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
    publishDir "${params.outdir}/10_mofa/${meta.id}", mode: params.publish_dir_mode

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
        Channel.fromPath("data/01_processed/bulk_rna_seq/*_expr_data*.rds")
    )
    P12_prepare_mofa_data(
        file_tuple("$dir/12_prepare_mofa_data.Rmd"),
        [:],
        P11_easier.out.collect()
    )
    P13_run_mofa(
        file_tuple("$dir/13_run_mofa.R"),
        P12_prepare_mofa_data.out.artifacts.flatten().filter{
            it -> it.simpleName.startsWith("mofa_")
        }
    )
    P14_mofa_analysis(
        file_tuple("$dir/14_mofa_analysis.Rmd"),
        [features_dir: "./", easier_dir: "./", tmb_dir: "./"],
        P12_prepare_mofa_data.out.artifacts.mix(
            Channel.of(file("$dir/helper_functions.R")),
            P11_easier.out.collect(),
            Channel.fromPath("data/01_processed/bulk_rna_seq/*TMB.rds").collect(),
            P13_run_mofa.out.models.collect()
        ).collect()
    )
}
