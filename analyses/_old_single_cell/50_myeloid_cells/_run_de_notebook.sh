#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_cpus=1
export MKL_NUM_cpus=1
export OPENBLAS_NUM_cpus=1

jupyter nbconvert 54_de_analysis.ipynb --execute --to html
