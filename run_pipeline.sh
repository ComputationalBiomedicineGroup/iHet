#!/bin/bash
# workaround conda issues on my system
source $HOME/.bashrc

export NXF_VER=22.04.5

nextflow run main.nf -w /scratch/nsclc_heterogeneity/work
