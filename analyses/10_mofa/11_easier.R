#!/usr/bin/env Rscript
"
Prepare Feature data with EASIER.

Usage:
  11_easier.R <INPUT_FILE.rds> <immscore_file.rds> <OUTPUT_FILE.rds>

" -> doc

library(easier)
library(immunedeconv)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
immscore_file <- args[2]
output_file <- args[3]

getFeat <- function(dataobj,
                    cancertype = "LUAD",
                    remove.genes.ICB_proxies = FALSE,
                    onlyEasier = TRUE,
                    epic = FALSE,
                    proxies = TRUE,
                    zscoredata = "../../tables/easier_Zscores/immscoreZ.rds") {
  immscoreZ <- readRDS(zscoredata)

  # List of immune response scores to be computed
  selscores <- c(
    "CYT",
    "Roh_IS",
    "chemokines",
    "Davoli_IS",
    "IFNy",
    "Ayers_expIS",
    "Tcell_inflamed",
    "resF_down",
    "TLS"
  )

  # List of datasets in the input object
  datasets <- names(dataobj[[1]])

  # Initialization of the output object
  dataobj$tf <- dataobj$pathway <- vector(mode = "list", length = length(datasets))
  names(dataobj$tf) <- names(dataobj$pathway) <- datasets

  for (dataset in datasets) {
    cat("Computing cell fractions, TF scores, pathway scores, and proxies of response from ", dataset, " data...\n", sep = "")

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
      RNA_counts = dataobj$count[[dataset]],
      remove_sig_genes_immune_response = remove.genes.ICB_proxies
    )

    # Quantification of TF activity scores with DoRothEA
    TF <- compute_TF_activity(
      RNA_tpm = dataobj$tpm[[dataset]],
    )

    # Storing of the computed features in the output object
    dataobj$cellfrac[[dataset]] <- cellfrac
    dataobj$pathway[[dataset]] <- pathway$scores
    dataobj$tf[[dataset]] <- TF$scores

    proxies.mat <- compute_scores_immune_response(
      RNA_tpm = dataobj$tpm[[dataset]],
      selected_scores = selscores
    )

    # Rotation of the "Chemokines" score to agree with "CYT" direction
    chemosign <- sign(cor(proxies.mat[, "CYT"], proxies.mat[, "chemokines"]))
    if (chemosign < 0) proxies.mat[, "chemokines"] <- -proxies.mat[, "chemokines"]

    # Computation of the median scaled response
    immscoreZ <- immscoreZ[match(selscores, immscoreZ$score), ]
    response <- proxies.mat[, match(selscores, colnames(proxies.mat))]
    response <- ((t(response) - immscoreZ$mean) / immscoreZ$sd)
    response <- apply(response, 2, median)
    proxies.mat <- as.data.frame(proxies.mat)
    proxies.mat$response <- response

    dataobj$immresp[[dataset]] <- proxies.mat
  }

  # Removal of gene expression data
  if (onlyEasier == TRUE) {
    dataobj <- dataobj[c("cellfrac", "pathway", "tf", "immresp")]
  }

  return(dataobj)
}

cancertype <- "NSCLC"
expr_obj <- readRDS(input_file)

easier_obj <- getFeat(expr_obj, cancertype = cancertype, epic = TRUE, onlyEasier = FALSE, zscoredata = immscore_file)

saveRDS(easier_obj, file = output_file)
