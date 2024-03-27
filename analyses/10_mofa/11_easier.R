#!/usr/bin/env Rscript
"
Prepare Feature data with EASIER.

Usage:
  11_easier.R <INPUT_FILE.rds> <OUTPUT_FILE.rds> <regulon_net> <TCGA_file.rds>

" -> doc

library(easier)
library(immunedeconv)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
regulon_net <- args[3]
tcga_file <- args[4]

getFeat <- function(dataobj,
                    cancertype = "LUAD",
                    remove.genes.ICB_proxies = FALSE,
                    onlyEasier = TRUE,
                    epic = FALSE,
                    proxies = TRUE,
                    tcga_obj) {
  # List of datasets in the input object
  datasets <- names(dataobj[[1]])

  # Initialization of the output object
  dataobj$tf <- dataobj$pathway <- vector(mode = "list", length = length(datasets))
  names(dataobj$tf) <- names(dataobj$pathway) <- datasets

  for (dataset in datasets) {
    cat("Computing cell fractions, TF scores, pathway scores, and proxies of response from ", dataset, " data...\n", sep = "")

    tcga_mean <- apply(dplyr::bind_cols(tcga_obj$tpm), 1, mean)
    tcga_sd <- apply(dplyr::bind_cols(tcga_obj$tpm), 1, sd)
    genes_intersection <- intersect(rownames(tcga_obj$tpm[[1]]), rownames(dataobj$tpm[[dataset]]))

    log_tpm <- log1p(dataobj$tpm[[dataset]])
    log_tpm_scaled <- log_tpm[genes_intersection, ] - tcga_mean[genes_intersection] / tcga_sd[genes_intersection]

    # Quantification of immune cell fractions with quanTIseq
    cellfrac <- compute_cell_fractions(
      RNA_tpm = dataobj$tpm[[dataset]]
    )

    # Optional quantification of CAFs and endothelial cells qith EPIC
    if (epic) {
      epic.cellfrac <- deconvolute_epic(dataobj$tpm[[dataset]],
        tumor = TRUE,
        scale_mrna = TRUE
      )

      cellfrac <- cbind(
        cellfrac,
        t(epic.cellfrac[c("CAFs", "Endothelial"), rownames(cellfrac)])
      )
    }

    # Quantification of pathway activity scores with RPOGENy
    pathway <- compute_pathway_activity(
      RNA_tpm = log_tpm_scaled,
      remove_sig_genes_immune_response = remove.genes.ICB_proxies
    )

    # Quantification of TF activity scores with DoRothEA
    TF <- compute_TF_activity(
      RNA_tpm = log_tpm_scaled,
      regulon_net = regulon_net
    )

    immscore <- compute_scores_immune_response(
      RNA_tpm = dataobj$tpm[[dataset]],
    )

    # Storing of the computed features in the output object
    dataobj$cellfrac[[dataset]] <- cellfrac
    dataobj$pathway[[dataset]] <- pathway
    dataobj$tf[[dataset]] <- TF
    dataobj$immresp[[dataset]] <- immscore
  }

  # Removal of gene expression data
  if (onlyEasier == TRUE) {
    dataobj <- dataobj[c("cellfrac", "pathway", "tf", "immresp")]
  }

  return(dataobj)
}

cancertype <- "NSCLC"
expr_obj <- readRDS(input_file)
tcga_obj <- readRDS(tcga_file)

easier_obj <- getFeat(expr_obj, cancertype = cancertype, epic = TRUE, onlyEasier = FALSE, tcga_obj = tcga_obj)

saveRDS(easier_obj, file = output_file)
