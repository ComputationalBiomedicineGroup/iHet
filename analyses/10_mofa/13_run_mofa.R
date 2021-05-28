#!/usr/bin/env Rscript

source('helper_functions.R')

args = commandArgs(trailingOnly=TRUE)

data_file = args[1]

mofa_data = readRDS(data_file)

runMOFA(
  mofa_data %>% select(sample, feature, value, group, view),
  str_replace(data_file, "\\.rds", ".hdf5")
)
