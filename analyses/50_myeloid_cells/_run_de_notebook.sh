#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_cpus=1
export MKL_NUM_cpus=1
export OPENBLAS_NUM_cpus=1

jupyter nbconvert myeloid_cluster.ipynb --execute --to html
