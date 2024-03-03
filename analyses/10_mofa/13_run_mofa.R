#!/usr/bin/env Rscript

library(conflicted)
library(dplyr)
library(stringr)
library(MOFA2)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]


#' train the MOFA model
runMOFA <- function(MOFAdata, out_file) {
  MOFAobject <- create_mofa(MOFAdata)

  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  train_opts <- get_default_training_options(MOFAobject)
  model_opts$num_factors <- 7 # number of samples is very small for learning 10 factors, it should not exceed ~7
  train_opts$convergence_mode <- "medium"
  train_opts$maxiter <- 7000
  # train_opts$drop_factor_threshold <- 0

  set.seed(1234)

  # Run MOFA
  MOFAobject <- prepare_mofa(object = MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts)

  model <- run_mofa(MOFAobject, outfile = out_file, use_basilisk = FALSE)

  set.seed(NULL)

  return(model)
}


mofa_data <- readRDS(data_file)

runMOFA(
  mofa_data %>% select(sample, feature, value, group, view),
  str_replace(data_file, "\\.rds", ".hdf5")
)
