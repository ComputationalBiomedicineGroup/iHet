#!/bin/bash
export NXF_VER=22.04.5
nextflow run main.nf -profile icbi_long -w $(readlink -f work)
