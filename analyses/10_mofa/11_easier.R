#!/usr/bin/env Rscript
"
Prepare Feature data with EASIER.

Usage:
  11_easier.R <INPUT_FILE.rds> <OUTPUT_FILE.rds>

" -> doc

library(easier)
library(immunedeconv)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

range01 <- function(x) {
  x <- (x - min(x)) / (max(x) - min(x))

  return(x)
}

getFeat <- function(dataobj, cancertype = "LUAD", remove.genes.ICB_proxies = FALSE, onlyEasier = TRUE, epic = FALSE) {
  goldstd <- c(
    "IMPRES",
    "MSI",
    "CTLA4",
    "PD1",
    "PDL1",
    "CYT",
    "Roh_IS",
    "chemokines",
    "Davoli_IS",
    "IFNy",
    "Ayers_expIS",
    "Tcell_inflamed",
    "RIR",
    "TLS"
  )

  datasets <- names(dataobj[[1]])

  dataobj$tf <- dataobj$pathway <- vector(mode = "list", length = length(datasets))
  names(dataobj$tf) <- names(dataobj$pathway) <- datasets

  for (dataset in datasets) {
    cat("Computing cell fractions, TF scores, pathway scores, and proxies of response from ",
      dataset, " data...\n",
      sep = ""
    )

    cellfrac <- compute_cell_fractions(
      RNA.tpm = dataobj$tpm[[dataset]]
    )

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

    pathway <- compute_pathways_scores(
      RNA.counts = dataobj$count[[dataset]],
      remove.genes.ICB_proxies = remove.genes.ICB_proxies
    )

    TF <- compute_TF_activity(
      RNA.tpm = dataobj$tpm[[dataset]],
      remove.genes.ICB_proxies = remove.genes.ICB_proxies
    )

    # Compute and aggregate proxies
    proxies <- compute_gold_standards(
      RNA.tpm =
        as.data.frame(dataobj$tpm[[dataset]]),
      list_gold_standards = goldstd,
      cancertype = cancertype
    )
    proxies.mat <- as.data.frame(matrix(unlist(proxies),
      nrow = length(proxies),
      byrow = T
    ))
    rownames(proxies.mat) <- names(proxies)
    colnames(proxies.mat) <- colnames(proxies[[1]])
    proxies.mat <- t(proxies.mat)

    # Rotate Chemokines proxy if not correlated with CYT
    chemosign <- sign(cor(proxies.mat[, "CYT"], proxies.mat[, "chemokines"]))
    if (chemosign < 0) proxies.mat[, "chemokines"] <- -proxies.mat[, "chemokines"]

    # Compute median scaled response
    selscores <- c(
      "CYT",
      "Roh_IS",
      "chemokines",
      "Davoli_IS",
      "IFNy",
      "Ayers_expIS",
      "Tcell_inflamed",
      "RIR",
      "TLS"
    )

    response <- proxies.mat[, match(selscores, colnames(proxies.mat))]
    response <- apply(response, 2, range01)
    response <- apply(response, 1, median)
    proxies.mat <- as.data.frame(proxies.mat)
    proxies.mat$response <- response

    dataobj$cellfrac[[dataset]] <- cellfrac
    dataobj$pathway[[dataset]] <- pathway$scores
    dataobj$tf[[dataset]] <- TF$scores
    dataobj$immresp[[dataset]] <- proxies.mat
  }

  if (onlyEasier == TRUE) {
    dataobj <- dataobj[c("cellfrac", "pathway", "tf", "immresp")]
  }

  return(dataobj)
}

cancertype <- "NSCLC"
expr_obj <- readRDS(input_file)

easier_obj <- getFeat(expr_obj, cancertype = cancertype, epic = TRUE)

saveRDS(easier_obj, file = output_file)
