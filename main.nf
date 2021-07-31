#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { W10_mofa } from "./subworkflows/mofa.nf"

workflow {
    W10_mofa()
}
