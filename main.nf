#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { W10_mofa } from "./subworkflows/mofa.nf"
include { W20_single_cell } from "./subworkflows/single-cell.nf"

workflow {
    // W10_mofa()
   W20_single_cell()
}
