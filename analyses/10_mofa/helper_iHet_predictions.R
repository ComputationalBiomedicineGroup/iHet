
# packages #
suppressPackageStartupMessages({
  library("remotes")
  library("conflicted")
  library("corrplot")
  library("gplots")
  library("RColorBrewer")
  library("Matrix")
  library("reshape2")
  library("ggrepel")
  library("GGally")
  library('easier')
  require('MOFA2') # read mofa models
  require('reshape2')
  require('immunedeconv')
  #library('dorothea')
  library("ROCR")
  library('dplyr')
  conflict_prefer("summarise", "dplyr")
  conflict_prefer("mutate", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("arrange", "dplyr")
  library('ggpubr')
  library('tibble')
  library('magrittr')
  library('ggtern')
  library("ggplot2")
  conflict_prefer("aes", "ggplot2")
  conflict_prefer("ggsave", "ggplot2")
  conflict_prefer("theme_bw", "ggplot2")
  conflict_prefer("theme_minimal", "ggplot2")
  library('ggforce')
  library('svglite')
  library('tidyr')
  library('plyr')
  library('purrr')
  library('ggfortify')
  library('stringr')
  library('rstatix')
  library('ComplexHeatmap')
})

rdbu10_palette <- c(
  "#67001F", "#B2182B", "#D6604D", "#F4A582",
  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
  "#4393C3", "#2166AC", "#053061"
)

colors_mfp <- c("#64a860", "#9970c1", "#cc545e", "#b98d3e")
names(colors_mfp) <- c("D", "F", "IE/F", "IE")

#' Rotate factors such that they are positively associated with CD8+ T-cell infiltration.
get_rotation <- function(weights) {
  diag(sign(weights[["Immune cells quantification"]]["CD8 T", ]))
}

#' Pre-process gene names
solve_annotation_issues <- function(gene_expr){

  # Some genes are causing issues due to approved symbols matching more than one gene
  genes_info <- easier:::reannotate_genes(cur_genes = rownames(gene_expr))
  ## Remove non-approved symbols
  non_na <- !is.na(genes_info$new_names)
  gene_expr <- gene_expr[non_na, ]
  genes_info <- genes_info[non_na, ]
  ## Remove entries that are withdrawn
  if (any(genes_info$new_names == "entry withdrawn")){
    gene_expr <- gene_expr[-which(genes_info$new_names == "entry withdrawn"), ]
    genes_info <- genes_info[-which(genes_info$new_names == "entry withdrawn"), ]
  }
  ## Identify duplicated new genes
  newnames_dup <- unique(genes_info$new_names[duplicated(genes_info$new_names)])
  newnames_dup_ind <- do.call(c, lapply(newnames_dup, function(X) which(genes_info$new_names == X)))
  newnames_dup <- genes_info$new_names[newnames_dup_ind]
  if(!is.null(newnames_dup_ind)) {
    ## Retrieve data for duplicated genes
    tmp <- gene_expr[genes_info$old_names[genes_info$new_names %in% newnames_dup],]
    ## Remove data for duplicated genes
    gene_expr <- gene_expr[-which(rownames(gene_expr) %in% rownames(tmp)),]
    ## Aggregate data of duplicated genes
    dup_genes <- genes_info$new_names[which(genes_info$new_names %in% newnames_dup)]
    names(dup_genes) <- rownames(tmp)
    if (anyDuplicated(newnames_dup)){
      tmp2 <- stats::aggregate(tmp, by = list(dup_genes), FUN = "mean")
      rownames(tmp2) <- tmp2$Group.1
      tmp2$Group.1 <- NULL
    }
    # Put data together
    gene_expr <- rbind(gene_expr, tmp2)
  }else{
    ## Retrieve data for duplicated genes
    rownames(gene_expr[genes_info$old_names,]) <- genes_info$new_names
  }
  return(gene_expr)
}

#' Compute correlation between features and HIF Density Fibroblasts
compute_corr <- function(my_data, hif, cor_method="spearman", dataset_names, pancancer.hif.mfp.cols) {

  # combine the datasets of interest beforehand (NSCLC= LUAD + LUSC)
  if(dataset_names == "pan-cancer"){
    my_data.df <- do.call(rbind, lapply(my_data, function(X) return(X)))
  }else{
    my_data.df <- do.call(rbind, my_data[dataset_names])
  }

  N <- ncol(my_data.df)

  zz <- strsplit(rownames(my_data.df)[1], '.', fixed = TRUE)
  if (length(zz[[1]])==5){
    my_data.df$patient <- sapply(rownames(my_data.df), function(x){
      paste(strsplit(x, '.', fixed = TRUE)[[1]][2:4],collapse ="-")
    })
  }else{
    my_data.df <- as.data.frame(my_data.df)
    my_data.df$patient <- sapply(rownames(my_data.df), function(x){
      paste(strsplit(x, '.', fixed = TRUE)[[1]][1:3],collapse ="-")
    })
  }

  my_data_hif.df <- join(my_data.df, pancancer.hif.mfp.cols[,c("patient", hif)])

  res <- do.call(rbind, lapply(colnames(my_data.df)[1:N], function(my_feature){
    # cat(my_feature, "\n")
    res_cor_test <- cor.test(my_data_hif.df[,my_feature], my_data_hif.df[,hif], method = cor_method)
    res_summary <- data.frame(cor=res_cor_test$estimate, p_val=res_cor_test$p.value, feature=my_feature, hif=hif, method=cor_method)
  }))

  res$FDR <- p.adjust(res$p_val, method = "BH")

  return(res)
}

#' Import validation datasets
import_dataset <- function(dataset) {
  myfilename <- file.path(params$validation_data, paste0("all_", dataset, ".RData"))

  if (length(grep("ontreatment", myfilename))>0){
    myfilename <- myfilename[-grep("ontreatment", myfilename)]
  }

  all_dataset <- loadRData(myfilename)

  prediction_easier <- apply(all_dataset$pred_easier,1,median)
  prediction_tasks <- apply(scale(all_dataset$pred_tasks),1,median)
  #prediction_tasks <- all_dataset$pred_tasks$response
  names(prediction_tasks) <- rownames(all_dataset$pred_tasks)

  TGFb <- all_dataset$computed_features[["pathway_activity"]][,"TGFb"]
  names(TGFb) <- rownames(all_dataset$computed_features[["pathway_activity"]])

  EGFR <- all_dataset$computed_features[["pathway_activity"]][,"EGFR"]
  names(EGFR) <- rownames(all_dataset$computed_features[["pathway_activity"]])

  if (dataset=="kim"){
    response <- all_dataset$clinical_response$Best_Response
    names(response) <- rownames(all_dataset$clinical_response)
    # TMB <- all_dataset$TMB$nonsyn
    # names(TMB) <- rownames(all_dataset$TMB)
    # TMB <- TMB[!is.na(TMB)]

    TMB <- all_dataset$TMB$Num_SNVs
    names(TMB) <- rownames(all_dataset$TMB)
    TMB[TMB=="Low"]=1
    TMB[TMB=="Mod"]=2
    TMB[TMB=="High"]=3
  }else if(dataset %in% c("gideauslanderpd1", "gideauslanderpd1on")){
    response <- all_dataset$clinical_response$BOR
    names(response) <- all_dataset$clinical_response$Sample
    response <- gsub(" ","", response)
    response <- gsub("CR|PR|MR|PRCR","R", response)
    response <- gsub("SD|PD|PD ","NR", response)
    TMB <-rep(NA, length(response))
    names(TMB) <- names(response)
  }else if(dataset=="zhao"){
    response <- all_dataset$clinical_response$Response.to.PD1.
    names(response) <- rownames(all_dataset$clinical_response)
    response <- gsub("Yes","R", response)
    response <- gsub("No","NR", response)
    TMB <-as.numeric(all_dataset$TMB[[1]])
    names(TMB) <- rownames(all_dataset$TMB)
    TMB <- TMB[!is.na(TMB)]
  }else if(dataset=="riaz"){
    response <- all_dataset$clinical_response$Response
    names(response) <- rownames(all_dataset$clinical_response)
    # response <- response[all_dataset$clinical_response$Cohort == "NIV3-NAIVE"]
    response <- gsub(" ","", response)
    response <- gsub("CR|PR|MR|PRCR","R", response)
    response <- gsub("SD|PD|PD ","NR", response)
    TMB <- all_dataset$TMB[[1]]
    names(TMB) <- rownames(all_dataset$TMB)
    TMB <- TMB[!is.na(TMB)]
  }else if(dataset=="hugo"){
    response <- all_dataset$clinical_response[[1]]
    names(response) <- rownames(all_dataset$clinical_response)
    levels(response)[levels(response)=="Progressive Disease"] <- "NR"
    levels(response)[levels(response)=="Partial Response"] <- "R"
    levels(response)[levels(response)=="Complete Response"] <- "R"
    response <- as.character(response)
    names(response) <- rownames(all_dataset$clinical_response)

    TMB <- all_dataset$TMB$TotalNonSyn
    names(TMB) <- rownames(all_dataset$TMB)
    # TMB <- all_dataset$TMB$TotalNonSyn_Exp
  }else if(dataset=="mariathasan"){
    response <- all_dataset$clinical_response$`Best Confirmed Overall Response`
    names(response) <- rownames(all_dataset$clinical_response)
    response <- response[response!="SD"]
    response <- response[response!="PR"]
    response <- gsub("CR|PR|MR|PRCR","R", response)
    response <- gsub("SD|PD|PD ","NR", response)

    TMB <- all_dataset$TMB[[1]]
    names(TMB) <- rownames(all_dataset$TMB)
    TMB <- TMB[!is.na(TMB)]
    # ImmPheno <- all_dataset$clinical_response$`Immune phenotype`
    # names(ImmPheno) <- rownames(all_dataset$clinical_response)
    # ImmPheno <- ImmPheno[!is.na(ImmPheno)]
  }else if(dataset=="liu"){
    response <- all_dataset$clinical_response$BR
    response <- gsub("CR|PR|MR|PRCR","R", response)
    response <- gsub("SD|PD|PD ","NR", response)
    names(response) <- rownames(all_dataset$clinical_response)

    TMB <- all_dataset$TMB[[1]]
    names(TMB) <- rownames(all_dataset$TMB)
  }else if(dataset=="snyder"){
    # response <- all_dataset$clinical_response[, "Best Response RECIST 1.1"]
    response <- all_dataset$clinical_response[, "Best response mRECIST"]
    names(response) <- all_dataset$clinical_response[,"ID"]
    response <- response[response!="SD"]
    response <- response[response!="PR"]
    response <- gsub("CR|PR|MR|PRCR","R", response)
    response <- gsub("SD|PD|PD ","NR", response)

    response <- response[response=="R" | response=="NR"]
    TMB <- all_dataset$TMB[,"missense_snv_count"]
    names(TMB) <- rownames(all_dataset$TMB)
  }else if(dataset=="jung"){
    # response <- all_dataset$clinical_response[, "Best Response RECIST 1.1"]
    response <- all_dataset$clinical_response[, "Clinical.benefit"]
    names(response) <- rownames(all_dataset$clinical_response)
    response <- gsub(" ","", response)
    response <- gsub("DCB","R", response)
    response <- gsub("NDB","NR", response)

    TMB <- as.vector(all_dataset$TMB[["Mutation.burden"]])
    names(TMB) <- rownames(all_dataset$TMB)
  }else if(dataset=="cho"){
    # response <- all_dataset$clinical_response[, "Best Response RECIST 1.1"]
    response <- all_dataset$clinical_response[, "Best.response"]
    names(response) <- rownames(all_dataset$clinical_response)
    response <- gsub("CR|PR|MR|PRCR","R", response)
    response <- gsub("SD|PD|PD ","NR", response)

    TMB <-rep(NA, length(response))
    names(TMB) <- names(response)

    HIST <- all_dataset$clinical_response$Histology
    HIST <- gsub("Squamous cell carcinoma ", "SQUAMOUS", HIST)
    HIST <- gsub("Adenocarcinoma ", "NON-SQUAMOUS", HIST)
    names(HIST) <- rownames(all_dataset$clinical_response)
    HIST <- HIST[names(response)]

  }else if(dataset=="mcdermott"){
    response <- all_dataset$clinical_response[, "BestResponse"]
    names(response) <- rownames(all_dataset$clinical_response)
    response <- gsub("CR|PR|MR|PRCR","R", response)
    response <- gsub("SD|PD|PD ","NR", response)

    TMB <- as.numeric(all_dataset$TMB$TMB)
    names(TMB) <- names(response)

  }else if(dataset %in% c("poplar", "oak")){
    response <- all_dataset$clinical_response[, "BCOR"]
    names(response) <- rownames(all_dataset$clinical_response)

    # remove NE, SD, ""
    response <- response[!response %in% c("NE", "")]
    response <- gsub("CR|PR","R", response)
    response <- gsub("PD","NR", response)

    TMB <- all_dataset$TMB[, "tTMB"]
    names(TMB) <- rownames(all_dataset$TMB)
    TMB <- TMB[names(response)]

    HIST <-all_dataset$clinical_response$HIST
    names(HIST) <- rownames(all_dataset$clinical_response)
    HIST <- HIST[names(response)]
  }

  if (dataset %in% c("cho", "poplar", "oak")){
    return(list(all_dataset=all_dataset,
                response=response,
                prediction_tasks=prediction_tasks,
                prediction_easier=prediction_easier,
                TMB=TMB,
                TGFb=TGFb,
                HIST=HIST,
                EGFR=EGFR))
  }else{
    return(list(all_dataset=all_dataset,
                response=response,
                prediction_tasks=prediction_tasks,
                prediction_easier=prediction_easier,
                TMB=TMB,
                TGFb=TGFb,
                EGFR=EGFR))
  }
}

#' Arrange a numeric vector in three tertiles
categorize <- function(x, thresholds=NULL, labels = c(1,2,3)) {
  vTert = quantile(x, c(0:3/3))
  # print(vTert)
  if (is.null(thresholds)){
    x.out = cut(x, vTert, include.lowest = TRUE, labels = labels)
  }else if ((is.numeric(thresholds)) & (length(thresholds)==2)){
    # x.out = cut(x, c(min(x, na.rm = T), 100, 400, max(x, na.rm = T)), include.lowest = T, labels = c(1, 2, 3))
    x.out = cut(x, c(min(x, na.rm = T), thresholds[1], thresholds[2], max(x, na.rm = T)),
                include.lowest = TRUE, labels = labels)
  }
  x.out <- as.numeric(as.character(x.out))
  names(x.out) <- names(x)
  return(x.out)
}

get_bootstrap_samples <- function(dir, model){
  bootstrap_samples <- lapply(as.character(1:100), function(run){

    # _id refers to the row number of that sample
    run_id <- paste0(paste(rep("0", 4-nchar(run)), collapse = ""), run)
    file_to_load <- paste0(paste("mofa_boot", model, run_id, sep = "_"), ".rds")
    model_bootstrap_run = readRDS(file.path(dir, file_to_load))

    if (model == "NSCLC"){
      boot_sample1 <- unique(subset(model_bootstrap_run, group == "LUAD")$sample)
      boot_sample2 <- unique(subset(model_bootstrap_run, group == "LUSC")$sample)
      boot_sample <- list(LUAD = boot_sample1, LUSC= boot_sample2)
    }else{
      boot_sample <- unique(model_bootstrap_run$sample)
    }
    return(boot_sample)
  })

  return(bootstrap_samples)
}

get_iHet_associated_weights <- function(use_bootstrap_weights=FALSE,
                                        model="NSCLC",
                                        scale_mofa_weights = FALSE,
                                        feat_cor_hifs){

  # Load knowledge of features correlated with information about the density of fibroblasts in tumor region
  feat_cor_tmp <- feat_cor_hifs
  rownames(feat_cor_tmp) <- feat_cor_tmp$feature
  feat_cor_tmp <- feat_cor_tmp[!feat_cor_tmp$feature %in% c("Other", "DC", "Monocyte"),]
  feat_sign <- feat_cor_tmp$feature[which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<=0.05)]

  # 1. get MOFA weights
  # 2. Adjust sign of the factor 1 to have consistent positive
  # association with immune response across bootstrap runs (get_rotation())

  mofa_boot_weights <- lapply(as.character(1:100), function(run){

    run_id <- paste0(paste(rep("0", 4-nchar(run)), collapse = ""), run)
    file <- paste0(paste("mofa_boot", model, run_id, sep = "_"), ".hdf5")

    boot_model <- MOFA2::load_model(file = file.path(params$mofa_dir, file), load_data = FALSE)
    boot_weights <- MOFA2::get_weights(boot_model, scale = scale_mofa_weights)

    return(boot_weights)
  })

  ## Bootstrap
  if (use_bootstrap_weights){

    F1_weights_df <- lapply(1:length(mofa_boot_weights), function(dataset_id) {
      weights <- mofa_boot_weights[[dataset_id]]
      tmp <- lapply(names(weights), function(modality) {
        weights[[modality]] %*% get_rotation(weights) %>%
          as_tibble(rownames = "feature") %>%
          mutate(view = modality, run = dataset_id) %>%
          rename_with(function(x) {
            str_replace(x, "V", "factor ")
          }, starts_with("V"))
      }) %>%
        bind_rows() %>%
        arrange(view, feature)
    }) %>%
      bind_rows() %>%
      select("feature", "factor 1", "view", "run") %>%
      dplyr::rename("weight" = "factor 1")

  }else{ ## Median

    median_weights <- lapply(1:length(mofa_boot_weights), function(dataset_id) {
      weights <- mofa_boot_weights[[dataset_id]]
      tmp <- lapply(names(weights), function(modality) {
        weights[[modality]] %*% get_rotation(weights) %>%
          as_tibble(rownames = "feature") %>%
          mutate(view = modality) %>%
          rename_with(function(x) {
            str_replace(x, "V", "factor ")
          }, starts_with("V"))
      }) %>%
        bind_rows() %>%
        arrange(view, feature)
    }) %>%
      bind_rows() %>%
      group_by(view, feature) %>%
      summarise_all(list(median = median, mad = mad)) %>%
      ungroup() %>%
      rename_with(~ gsub("_median", "", .x), ends_with("_median"))

    F1_weights_df <- median_weights %>%
      dplyr::select("feature", "factor 1", "view") %>%
      dplyr::rename("weight" = "factor 1")
  }

  # iHet_excl weights (set cor features to 0)
  weights_excl <- F1_weights_df %>%
    dplyr::filter(feature %in% feat_sign) %>%
    dplyr::mutate(weight=0, score = "iHet_excl")

  # iHet_rev (negate only cor features)
  weights_rev <- F1_weights_df %>%
    dplyr::filter(feature %in% feat_sign) %>%
    dplyr::mutate(weight=-weight, score = "iHet_rev")

  # iHet
  weights_iHet <- F1_weights_df %>%
    dplyr::mutate(score = "iHet")

  all_weights <- rbind(weights_iHet, weights_excl, weights_rev)
  all_weights_df <- all_weights %>%
    tidyr::pivot_wider(id_cols = c(feature, view), names_from = score, values_from = weight)

  all_weights_df <- all_weights_df %>%
    dplyr::rename("feature_type" = "view",
                  "iHet_weights" = "iHet",
                  "iHet_excl_weights" = "iHet_excl",
                  "iHet_rev_weights" = "iHet_rev")

  all_weights_df <- all_weights_df %>%
    dplyr::mutate(iHet_excl_weights=ifelse(is.na(iHet_excl_weights), iHet_weights, iHet_excl_weights),
                  iHet_rev_weights=ifelse(is.na(iHet_rev_weights), iHet_weights, iHet_rev_weights))

  return(all_weights_df)
}


#' Compute iHet scores
compute_iHet <- function(dataset,
                         model = c("BLCA", "BRCA", "CESC","CRC", "GBM", "HNSC", "KIRC", "KIRP",
                                   "LIHC", "LUAD", "LUSC", "NSCLC", "OV", "PAAD", "PRAD",
                                   "SKCM", "STAD", "THCA", "UCEC"),
                         use_bootstrap_weights=FALSE,
                         method=c("iHet", "iHet_rev", "iHet_excl"),
                         scale_mofa_weights = FALSE,
                         regulon_net = c("collectri", "dorothea")){

  # match arguments
  model <- match.arg(model, c("BLCA", "BRCA", "CESC","CRC", "GBM", "HNSC", "KIRC", "KIRP",
                              "LIHC", "LUAD", "LUSC", "NSCLC", "OV", "PAAD", "PRAD",
                              "SKCM", "STAD", "THCA", "UCEC", "JiaSharma"))

  # 1. get MOFA weights
  # 2. Adjust sign of the factor 1 to have consistent positive
  # association with immune response across bootstrap runs (get_rotation())

  ## Bootstrap
  if (use_bootstrap_weights){

    mofa_boot_weights <- lapply(as.character(1:100), function(run){

      run_id <- paste0(paste(rep("0", 4-nchar(run)), collapse = ""), run)
      file <- paste("mofa_boot", model, run_id, sep = "_")
      file_to_load <- paste0(file, ".hdf5")

      boot_model <- MOFA2::load_model(file = file.path(params$mofa_dir, file_to_load), load_data = FALSE)
      boot_weights <- MOFA2::get_weights(boot_model, scale = scale_mofa_weights)

      return(boot_weights)
    })

    F1_weights_df <- lapply(1:length(mofa_boot_weights), function(dataset_id) {
      weights <- mofa_boot_weights[[dataset_id]]
      tmp <- lapply(names(weights), function(modality) {
        weights[[modality]] %*% get_rotation(weights) %>%
          as_tibble(rownames = "feature") %>%
          mutate(view = modality, run = dataset_id) %>%
          rename_with(function(x) {
            str_replace(x, "V", "factor ")
          }, starts_with("V"))
      }) %>%
        bind_rows() %>%
        arrange(view, feature)
    }) %>%
      bind_rows() %>%
      select("feature", "factor 1", "view", "run") %>%
      dplyr::rename("weight" = "factor 1")

  ## Median
  }else{
    ## FF provided file: (paste0(repo_dir, "/data/iHet_data/median_weights.rds"))
    ## r=0.92 between these computed weights and FF's
    median_weights <- lapply(1:length(mofa_boot_weights), function(dataset_id) {
      weights <- mofa_boot_weights[[dataset_id]]
      tmp <- lapply(names(weights), function(modality) {
        weights[[modality]] %*% get_rotation(weights) %>%
          as_tibble(rownames = "feature") %>%
          mutate(view = modality) %>%
          rename_with(function(x) {
            str_replace(x, "V", "factor ")
          }, starts_with("V"))
      }) %>%
        bind_rows() %>%
        arrange(view, feature)
    }) %>%
      bind_rows() %>%
      group_by(view, feature) %>%
      summarise_all(list(median = median, mad = mad)) %>%
      ungroup() %>%
      rename_with(~ gsub("_median", "", .x), ends_with("_median"))

    F1_weights_df <- median_weights$`factor 1`
    names(F1_weights_df) <- median_weights$feature
  }

  #  --- compute features --- #

  # immune cell fractions ------------------------------------------------------

  ## quantiseq (easier returns already CD4 T = CD4 T + Treg)
  cellfrac <- easier::compute_cell_fractions(RNA_tpm = dataset$tpm)

  ### Change cell type 'CD8+ T' for 'CD8 T'
  cellfrac <- cellfrac %>% dplyr::rename("CD8 T" = "CD8+ T")

  ## Remove some cell types
  if(all(c("DC", "Monocyte") %in% F1_weights_df$feature)){
    cellfrac <- cellfrac[,colnames(cellfrac) != "Other"]
  }else{
    cellfrac <- cellfrac[,!(colnames(cellfrac) %in% c("Other", "DC", "Monocyte"))]
  }

  ## epic (quantification of cafs and endothelial cells)
  epic_cellfrac <- immunedeconv::deconvolute_epic(gene_expression_matrix = dataset$tpm,
                                                  tumor = TRUE,
                                                  scale_mrna = TRUE)
  cellfrac <- cbind(cellfrac,
                    t(epic_cellfrac[c("CAFs", "Endothelial"), rownames(cellfrac)]))

  ## log transformation
  cellfrac <- log10(cellfrac * 100 + 0.001)

  # TF activity ----------------------------------------------------------------

  tf <- easier::compute_TF_activity(RNA_tpm = dataset$tpm, regulon_net = regulon_net)

  # pathway activity -----------------------------------------------------------

  path <- easier::compute_pathway_activity(RNA_counts = dataset$counts,
                                           remove_sig_genes_immune_response = FALSE)

  #  --- feature scaling based on TCGA MOFA training values --- #

  mofa_boot_feature_stats <- readRDS(file =
                                       file.path(params$artifact_dir,
                                                 paste0(regulon_net, "_TCGA_boot_features_stats.rds")))

  boot_features_stats <- mofa_boot_feature_stats[[model]]

  if (model == "NSCLC"){

    norm_boot_features <- lapply(c("LUAD", "LUSC"), function(XX){

      norm_boot_features <- lapply(1:length(boot_features_stats), function(run){

        boot_features_stats_run <- subset(boot_features_stats[[run]], group == XX)

        # cell fractions
        stats_cellfrac <- boot_features_stats_run %>%
          dplyr::filter(feature %in% colnames(cellfrac)) %>%
          dplyr::arrange(factor(feature, levels = colnames(cellfrac)))

        cellfrac_centered <- sweep(cellfrac, 2, stats_cellfrac$feature_mean, FUN = "-")

        # tf activity
        stats_tf <- boot_features_stats_run %>%
          dplyr::filter(feature %in% colnames(tf)) %>%
          dplyr::arrange(factor(feature, levels = colnames(tf)))

        tf_centered <- sweep(tf, 2, stats_tf$feature_mean, FUN = "-")

        # pathway activity
        stats_path <- boot_features_stats_run %>%
          dplyr::filter(feature %in% colnames(path)) %>%
          dplyr::arrange(factor(feature, levels = colnames(path)))

        path_centered <- sweep(path, 2, stats_path$feature_mean, FUN = "-")
        path_scaled <- sweep(path_centered, 2, stats_path$feature_std, FUN = "/")

        return(list(path = as.matrix(path_scaled),
                    cellfrac = as.matrix(cellfrac_centered),
                    tf = as.matrix(tf_centered)))
      })
      return(norm_boot_features)
    })
    names(norm_boot_features) <- c("LUAD", "LUSC")

  }else{

    norm_boot_features <- lapply(1:length(boot_features_stats), function(run){

      # cell fractions
      stats_cellfrac <- boot_features_stats[[run]] %>%
        dplyr::filter(feature %in% colnames(cellfrac)) %>%
        dplyr::arrange(factor(feature, levels = colnames(cellfrac)))

      cellfrac_centered <- sweep(cellfrac, 2, stats_cellfrac$feature_mean, FUN = "-")

      # tf activity
      stats_tf <- boot_features_stats[[run]] %>%
        dplyr::filter(feature %in% colnames(tf)) %>%
        dplyr::arrange(factor(feature, levels = colnames(tf)))

      tf_centered <- sweep(tf, 2, stats_tf$feature_mean, FUN = "-")

      # pathway activity
      stats_path <- boot_features_stats[[run]] %>%
        dplyr::filter(feature %in% colnames(path)) %>%
        dplyr::arrange(factor(feature, levels = colnames(path)))

      path_centered <- sweep(path, 2, stats_path$feature_mean, FUN = "-")
      path_scaled <- sweep(path_centered, 2, stats_path$feature_std, FUN = "/")

      return(list(path = as.matrix(path_scaled),
                  cellfrac = as.matrix(cellfrac_centered),
                  tf = as.matrix(tf_centered)))
    })

  }

  feat_cor_hifs <- readRDS(file = file.path(params$artifact_dir,
                                            paste0(regulon_net, "_TCGA_NSCLC_cor_features_HIF_Fib_tumor.rds")))

  if (use_bootstrap_weights == FALSE){

    if (model == "NSCLC"){

      norm_boot_features <- reshape2::melt(norm_boot_features)
      colnames(norm_boot_features) <- c("patient", "feature", "value", "view", "run", "model")
      norm_boot_features <- norm_boot_features %>%
        dplyr::group_by(patient, feature, view, model) %>%
        dplyr::summarise(median_value = median(value))

    }else{

      norm_boot_features <- reshape2::melt(norm_boot_features)
      colnames(norm_boot_features) <- c("patient", "feature", "value", "view", "run")
      norm_boot_features <- norm_boot_features %>%
        dplyr::group_by(patient, feature, view) %>%
        dplyr::summarise(median_value = median(value))
    }
  }

  cat("Computing iHet score... \n")
  if (use_bootstrap_weights){

    # Just for NSCLC mofa model
    if (model == "NSCLC"){

      iHet_bootstrap <- lapply(c("LUAD", "LUSC"), function(XX){

        # compute iHet score
        iHet_bootstrap <- lapply(1:100, function(bot_sample){

          F1 <- subset(F1_weights_df, run %in% bot_sample)$weight
          names(F1) <- subset(F1_weights_df, run %in% bot_sample)$feature

          # combine features
          features <- cbind(norm_boot_features[[XX]][[bot_sample]]$cellfrac,
                            norm_boot_features[[XX]][[bot_sample]]$tf,
                            norm_boot_features[[XX]][[bot_sample]]$path)

          # match features
          cfeatures <- intersect(colnames(features), names(F1))
          features <- features[, cfeatures]
          F1_vec <- F1[cfeatures]
          features <- t(features)

          iHet_method <- lapply(method, function(method_x){
            # modify the F1 based on the defined method (the features association with hifs density CAFs)
            if (method_x != "iHet") {

              feat_cor_tmp <- feat_cor_hifs
              rownames(feat_cor_tmp) <- feat_cor_tmp$feature
              feat_cor_tmp <- feat_cor_tmp[cfeatures,]
              ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<=0.05)

              if (method_x == "iHet_excl") {
                iHet <- F1_vec[ix_sign] %*% features[ix_sign,]
              } else if (method_x == "iHet_rev") {
                F1_vec[ix_sign] <- - F1_vec[ix_sign]
                iHet <- F1_vec %*% features
              }

              iHet <- t(iHet[1,,drop=FALSE])

              # just iHet method

            } else{
              iHet <- F1_vec %*% features
              #iHet <- qr.solve(F1_vec, features)
              iHet <- t(iHet[1,,drop=FALSE])
            }
            return(iHet)
          })
          names(iHet_method) <- method
          return(iHet_method)
        })
        names(iHet_bootstrap) <- 1:100
        return(iHet_bootstrap)
      })
      names(iHet_bootstrap) <- c("LUAD", "LUSC")

      # For the rest of mofa models
    } else {

      # compute iHet score
      iHet_bootstrap <- lapply(1:100, function(bot_sample){

        F1 <- subset(F1_weights_df, run %in% bot_sample)$weight
        names(F1) <- subset(F1_weights_df, run %in% bot_sample)$feature

        # combine features
        features <- cbind(norm_boot_features[[bot_sample]]$cellfrac,
                          norm_boot_features[[bot_sample]]$tf,
                          norm_boot_features[[bot_sample]]$path)

        # match features
        cfeatures <- intersect(colnames(features), names(F1))
        features <- features[, cfeatures]
        F1_vec <- F1[cfeatures]
        features <- t(features)

        iHet_method <- lapply(method, function(method_x){
          # modify the F1 based on the defined method
          if (method_x != "iHet") {

            feat_cor_tmp <- feat_cor_hifs
            rownames(feat_cor_tmp) <- feat_cor_tmp$feature
            feat_cor_tmp <- feat_cor_tmp[cfeatures,]
            ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<=0.05)

            if (method_x == "iHet_excl") {
              iHet <- F1_vec[ix_sign] %*% features[ix_sign, , drop=FALSE]

            } else if (method_x == "iHet_rev") {
              F1_vec[ix_sign] <- - F1_vec[ix_sign]
              iHet <- F1_vec %*% features

            }

            iHet <- t(iHet[1,,drop=FALSE])
            return(iHet)

            # just iHet method

          } else{
            iHet <- F1_vec %*% features
            #iHet <- qr.solve(F1_vec, features)
            iHet <- t(iHet[1,,drop=FALSE])
          }
          return(iHet)
        })
        names(iHet_method) <- method
        return(iHet_method)
      })
      names(iHet_bootstrap) <- 1:100
    }

    # to data.frame
    iHet_bootstrap_df <- reshape2::melt(iHet_bootstrap)
    iHet_bootstrap_df$Var2 <- NULL

    if(model == "NSCLC"){
      colnames(iHet_bootstrap_df) <- c("patient", "score", "method", "run", "model")
    }else{
      colnames(iHet_bootstrap_df) <- c("patient", "score", "method", "run")
    }

    return(iHet_bootstrap_df)

  }else{
    # Median
    if (model == "NSCLC"){

      iHet_median <- lapply(c("LUAD", "LUSC"), function(XX){

        # combine features
        features_cellfrac <- norm_boot_features %>%
          dplyr::filter(view=="cellfrac" & model == XX) %>%
          tidyr::pivot_wider(id_cols = !c(view, model), names_from = feature, values_from = median_value) %>%
          tibble::column_to_rownames(var = "patient")

        features_tf <- norm_boot_features %>%
          dplyr::filter(view=="tf" & model == XX) %>%
          tidyr::pivot_wider(id_cols = !c(view, model), names_from = feature, values_from = median_value) %>%
          tibble::column_to_rownames(var = "patient")

        features_path <- norm_boot_features %>%
          dplyr::filter(view=="path" & model == XX) %>%
          tidyr::pivot_wider(id_cols = !c(view, model), names_from = feature, values_from = median_value) %>%
          tibble::column_to_rownames(var = "patient")

        features <- cbind(features_cellfrac,
                          features_tf,
                          features_path)

        # match features
        cfeatures <- intersect(colnames(features), names(F1_weights_df))
        features <- features[, cfeatures]
        F1_vec <- F1_weights_df[cfeatures]
        features <- t(features)

        iHet_method <- lapply(method, function(method_x){
          # modify the F1 based on the defined method
          if (method_x != "iHet") {

            feat_cor_tmp <- feat_cor_hifs
            rownames(feat_cor_tmp) <- feat_cor_tmp$feature
            feat_cor_tmp <- feat_cor_tmp[cfeatures,]
            ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<0.05)

            if (method_x == "iHet_excl") {
              iHet <- F1_vec[ix_sign] %*% features[ix_sign,]

            } else if (method_x == "iHet_rev") {

              F1_vec[ix_sign] <- - F1_vec[ix_sign]
              iHet <- F1_vec %*% features
            }

          }
          # Compute iHet
          iHet <- F1_vec %*% features
          iHet <- t(iHet[1,,drop=TRUE])
          return(iHet)
        })
        names(iHet_method) <- method
        return(iHet_method)
      })
      names(iHet_median) <- c("LUAD", "LUSC")

    }else{

      # combine features
      features_cellfrac <- norm_boot_features %>%
        dplyr::filter(view=="cellfrac") %>%
        tidyr::pivot_wider(id_cols = !view, names_from = feature, values_from = median_value) %>%
        tibble::column_to_rownames(var = "patient")

      features_tf <- norm_boot_features %>%
        dplyr::filter(view=="tf") %>%
        tidyr::pivot_wider(id_cols = !view, names_from = feature, values_from = median_value) %>%
        tibble::column_to_rownames(var = "patient")

      features_path <- norm_boot_features %>%
        dplyr::filter(view=="path") %>%
        tidyr::pivot_wider(id_cols = !view, names_from = feature, values_from = median_value) %>%
        tibble::column_to_rownames(var = "patient")

      features <- cbind(features_cellfrac,
                        features_tf,
                        features_path)

      # match features
      cfeatures <- intersect(colnames(features), names(F1_weights_df))
      features <- features[, cfeatures]
      F1_vec <- F1_weights_df[cfeatures]
      features <- t(features)

      iHet_method <- lapply(method, function(method_x){
        # modify the F1 based on the defined method
        if (method_x != "iHet") {

          feat_cor_tmp <- feat_cor_hifs
          rownames(feat_cor_tmp) <- feat_cor_tmp$feature
          feat_cor_tmp <- feat_cor_tmp[cfeatures,]
          ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<0.05)

          if (method_x == "iHet_excl") {
            iHet <- F1_vec[ix_sign] %*% features[ix_sign,]

          } else if (method_x == "iHet_rev") {

            F1_vec[ix_sign] <- - F1_vec[ix_sign]
            iHet <- F1_vec %*% features
          }

        }
        # Compute iHet
        iHet <- F1_vec %*% features
        iHet <- t(iHet[1,,drop=TRUE])
        return(iHet)
      })
      names(iHet_method) <- method
      iHet_median <- iHet_method
    }

    # to data.frame
    iHet_median_df <- reshape2::melt(iHet_median)
    iHet_median_df$Var1 <- NULL

    if (model == "NSCLC"){
      colnames(iHet_median_df) <- c("patient", "score", "method", "model")
    }else{
      colnames(iHet_median_df) <- c("patient", "score", "method")
    }
    return(iHet_median_df)
  }
}

#' Assess iHet performance
assess_predictive_performance <- function(iHet,
                                          response,
                                          use_bootstrap_weights=FALSE,
                                          method,
                                          plot_roc=FALSE,
                                          metric="roc"){

  if(metric=="precision_recall"){
    yaxis <- "prec"; yaxis_name <- "Precision"
    xaxis <- "rec"; xaxis_name <- "Recall"
    auc_metric <- "aucpr"
  }else if (metric == "roc"){
    yaxis <- "tpr"; yaxis_name <- "True positive rate"
    xaxis <- "fpr"; xaxis_name <- "False positive rate"
    auc_metric <- "auc"
  }

  # get dataset name
  dataset <- unique(iHet$dataset)
  nR <- table(response)["R"]
  nNR <- table(response)["NR"]

  if (plot_roc==TRUE) {

    graphics::par(
      cex.axis = 1.3,
      #mar = c(5, 4, 2, 8),
      col.lab = "black",
      pty = "s", xpd = FALSE
    )

    ROCR::plot(c(0,1), c(0,1), col="white", lwd = 1,
               type = "S", cex.lab = 1.3, xlab=xaxis_name, ylab=yaxis_name)

    if (metric=="roc") abline(a=0, b=1, lty = 3, lwd = 1, col = "antiquewhite4")
    if (metric=="precision_recall") abline(a=nR/(nR+nNR), b=0, lty = 3, lwd = 1, col = "antiquewhite4")

    title(main = paste0(dataset, " (R=", nR, "; NR=", nNR, ")"))
  }

  # set colors
  col <- c("#7570B3","#D95F02", "#72a659", "#ca568b", "gray85")
  names(col) <- c("iHet", "iHet_rev", "PDL1", "TMB", "ESTIMATE")
  col <- col[c("iHet", "iHet_rev", "PDL1", "TMB", "ESTIMATE")]

  #col=colorRampPalette(brewer.pal(8, "Set2"))(6)
  #names(col) <- c("IR", "iHet", "iHet_inverted", "iHet_removed", "iHet_exclusion")
  #col <- col[c("iHet", "iHet_inverted", "iHet_removed", "iHet_exclusion")]

  rp_df <- list()

  ## categorize gold standards
  if(any(method %in% c("PDL1", "TMB", "ESTIMATE"))){
    iHet <- iHet %>%
      dplyr::group_by(method) %>%
      dplyr::mutate(score = ifelse(method == "TMB", categorize(score), score)) %>%
      dplyr::mutate(score = ifelse(method == "PDL1", categorize(score), score)) %>%
      dplyr::mutate(score = ifelse(method == "ESTIMATE", categorize(score), score))
  }

  ## sort iHet (for correct transformation into matrix, prior to predict response) ##
  iHet <- iHet %>%
    dplyr::group_by(patient) %>%
    dplyr::arrange(patient)

  # select common samples
  if(use_bootstrap_weights) {

    rp_df$iHet <- iHet
    # creating a response matrix to evaluate bootstrap performance
    rp_df$response <- matrix(response,
                             nrow = length(response),
                             ncol= length(unique(iHet$run))) # 100
    rownames(rp_df$response) <- names(response)
  }else{

    rp_df$iHet <- iHet
    rp_df$iHet$run <- 1

    # creating a response matrix to evaluate bootstrap performance
    rp_df$response <- matrix(response,
                             nrow = length(response),
                             ncol= 1)
    rownames(rp_df$response) <- names(response)
  }

  # predict with iHet
  ROC_data <- lapply(method, function(yy){

    iHet_sel <- subset(rp_df$iHet, method == yy)
    iHet_mat <- matrix(iHet_sel$score,
                       nrow=length(unique(iHet_sel$patient)),
                       ncol=length(unique(iHet_sel$run)),
                       byrow = TRUE)

    rownames(iHet_mat) <- unique(iHet_sel$patient)
    colnames(iHet_mat) <- unique(iHet_sel$run)

    if(yy %in% c("PDL1", "TMB", "ESTIMATE")){
      pred <- ROCR::prediction(iHet_mat[,1], rp_df$response[,1], label.ordering = c("NR", "R"))
    }else{
      pred <- ROCR::prediction(iHet_mat, rp_df$response, label.ordering = c("NR", "R"))
    }
    roc_perf <- ROCR::performance(pred, measure = yaxis, x.measure = xaxis)
    auc_perf <- unlist(ROCR::performance(pred, auc_metric)@y.values)
    auc_perf <- matrix(auc_perf)
    if(yy %in% c("PDL1", "TMB", "ESTIMATE")) {
      rownames(auc_perf) <- 1
    }else{
      rownames(auc_perf) <- 1:100
    }
    return(list(roc_perf=roc_perf, auc_perf=auc_perf))
  })
  names(ROC_data) <- method

  if (plot_roc==TRUE){

    if(use_bootstrap_weights) {

      legend_text <- sapply(method, function(xx){
        # confidence intervals
        tmp <- ROC_data[[xx]]$auc_perf
        mean <- mean(tmp)
        stdev <-  sqrt(var(tmp))
        n     <- 100
        ciw   <- qt(0.975, n) * stdev / sqrt(n) # 95% CI

        # using ROCR package
        if (xx %in% c("PDL1", "TMB", "ESTIMATE")){
          ROCR::plot(ROC_data[[xx]]$roc_perf, lwd=3, cex.lab = 1.3, col=col[xx], type = "s", add = TRUE)
        }else{
          ROCR::plot(ROC_data[[xx]]$roc_perf, lwd=3, cex.lab = 1.3, col=col[xx], avg = "threshold", type = "s", add = TRUE)
          # bootstrap ROC curves
          plot(ROC_data[[xx]]$roc_perf, col=col[xx], lty = 2, lwd = 0.25, add = TRUE)
        }

        if(metric=="precision_recall"){

          if (xx %in% c("PDL1", "TMB", "ESTIMATE")){

            legend_text <- paste0(xx,", AUCPR=", round(mean, 2))

          }else{

            legend_text <- paste0(xx,", AUCPR=", round(mean, 2),
                                  " (95% CI: ", round(mean-ciw,3), ", ",
                                  round(mean+ciw,3), ")")
          }

        }else{

          if (xx %in% c("PDL1", "TMB", "ESTIMATE")){

          legend_text <- paste0(xx,", AUC=", round(mean, 2))

          }else{

            legend_text <- paste0(xx,", AUC=", round(mean, 2),
                                  " (95% CI: ", round(mean-ciw,3), ", ",
                                  round(mean+ciw,3), ")")
          }
        }
        return(legend_text)
      })

    }else{

      legend_text <- sapply(method, function(xx){
        tmp <- ROC_data[[xx]]$auc_perf

        # using ROCR package
        if (xx %in% c("PDL1", "TMB", "ESTIMATE")){
          ROCR::plot(ROC_data[[xx]]$roc_perf, lwd=2, cex.lab = 1.3, col=col[xx], type = "s", add = TRUE)
        }else{
          ROCR::plot(ROC_data[[xx]]$roc_perf, lwd=2, cex.lab = 1.3, col=col[xx], avg = "threshold", type = "s", add = TRUE)
        }

        if(metric=="precision_recall"){
          legend_text <- paste0(xx,", AUCPR=", round(tmp, 3))
        }else{
          legend_text <- paste0(xx,", AUC=", round(tmp, 3))
        }
        return(legend_text)
      })

    }

    if(metric=="precision_recall"){
      legend(x=0.55, y=1.05, xjust = 0,
             legend = c(legend_text),
             fill = col[method], bty="n", cex=1)
    }else{
      legend(x=0.55, y=0.17, xjust = 0,
             legend = c(legend_text),
             fill = col[method], bty="n", cex=1)
    }

    plot <- grDevices::recordPlot()

  }

  # collect AUC bootstrap values

  auc_df <- do.call(rbind, lapply(method, function (zz){

    ROC_data[[zz]]$roc_perf <- NULL
    tmp_auc_df <- reshape2::melt(ROC_data[[zz]])
    tmp_auc_df$Var2 <- NULL; tmp_auc_df$L1 <- NULL;
    tmp_auc_df$method <- zz
    return(tmp_auc_df)

  }))
  colnames(auc_df) <- c("run", "auc", "method")

  if (plot_roc==TRUE) {
    return(plot)
  }else{
    return(auc_df)
  }
}

#' Integrate iHet, iHet_excl and TMB
integrated_score <- function(myscores_df, grid=seq(0, 1, length.out = 19),
                             with_TMB = FALSE,
                             normalize = c("01_first", "01_after")) {

  dataset_name <- unique(myscores_df$dataset)
  file_name <- sapply(strsplit(dataset_name, split = " (", fixed = TRUE), head, 1)

  # TMB #
  if (with_TMB){
    if (length(unique(myscores_df$TMB))>3){
      myscores_df$TMB_cat <- categorize_TMB(myscores_df$TMB)
    }else{
      if(file_name %in% c("oak.nsclc", "oak.luad", "oak.lusc", "poplar.nsclc", "poplar.luad", "poplar.lusc", "oak")) {
        myscores_df$TMB_cat <- gsub("<16", 1, myscores_df$TMB)
        myscores_df$TMB_cat <- gsub(">=16", 2, myscores_df$TMB_cat)
        myscores_df$TMB_cat <- as.numeric(myscores_df$TMB_cat)
      }else{
        myscores_df$TMB_cat <- as.numeric(myscores_df$TMB)
      }
    }

    if (!all(is.na(myscores_df$TMB_cat))){

      if(length(myscores_df$patient) > length(unique(myscores_df$patient))){
        ## bootstrap
        myscores_df_1 <- myscores_df %>% filter(run == 1)
        pred <- prediction(myscores_df_1$TMB_cat, myscores_df_1$response, label.ordering = c("NR", "R"))
        AUC_TMB <- performance(pred, measure = "auc")
        AUC_TMB <- mean(unlist(AUC_TMB@y.values))
      }else{
        ## median
        pred <- prediction(myscores_df$TMB_cat, myscores_df$response)
        AUC_TMB <- performance(pred, measure = "auc")
        AUC_TMB <- AUC_TMB@y.values[[1]]
      }
      #cat("AUC_TMB: ", AUC_TMB, "\n")
    }
  }

  if (with_TMB) myscores_df$TMB_lin <- linear_func(myscores_df$TMB_cat, range = "01")

  if (normalize == "01_first"){
    if(length(myscores_df$patient) > length(unique(myscores_df$patient))){
      ## bootstrap
      myscores_df <- myscores_df %>%
        dplyr::group_by(run) %>%
        dplyr::mutate(iHet = linear_func(iHet, range = "01"),
                      iHet_excl = linear_func(iHet_excl, range = "01"))

    }else{
      ## median
      myscores_df$iHet <- linear_func(myscores_df$iHet, range = "01")
      myscores_df$iHet_excl <- linear_func(myscores_df$iHet_excl, range = "01")
    }

  }

  # compute predictions
  if(length(myscores_df$patient) > length(unique(myscores_df$patient))){ ## bootstrap

    # response matrix
    myresponse_mat <- myscores_df
    if (with_TMB) myresponse_mat$TMB <- myresponse_mat$TMB_cat <- myresponse_mat$TMB_lin <- NULL
    myresponse_mat$iHet_rev <- myresponse_mat$iHet_excl <- myresponse_mat$iHet <- NULL
    myresponse_mat$dataset <- NULL
    myresponse_mat <- myresponse_mat %>%
      tidyr::pivot_wider(names_from = run, values_from = response) %>%
      as.data.frame()
    rownames(myresponse_mat) <- myresponse_mat$patient
    myresponse_mat$patient <- NULL

    # iHet_Excl pred matrix
    myscores_mat <- myscores_df
    if (with_TMB) myscores_mat$TMB <- myscores_mat$TMB_cat <- myscores_mat$TMB_lin <- NULL
    myscores_mat$iHet_rev <- myscores_mat$iHet <- NULL
    myscores_mat$response <- myscores_mat$dataset <- NULL
    myscores_mat <- myscores_mat %>%
      tidyr::pivot_wider(names_from = run, values_from = iHet_excl) %>%
      as.data.frame()
    rownames(myscores_mat) <- myscores_mat$patient
    myscores_mat$patient <- NULL

    pred <- ROCR::prediction(myscores_mat, myresponse_mat, label.ordering = c("NR", "R"))
    AUC_iHet_excl <- performance(pred, measure = "auc")
    AUC_iHet_excl <- unlist(AUC_iHet_excl@y.values)
    AUC_iHet_excl_CI <- qt(0.975, 100) * sqrt(var(AUC_iHet_excl)) / sqrt(100) # 95% CI
    #cat("Average (bootstrap) AUC_iHet_excl: ", AUC_iHet_excl, "\n")

    # iHet pred matrix
    myscores_mat <- myscores_df
    if (with_TMB) myscores_mat$TMB <- myscores_mat$TMB_cat <- myscores_mat$TMB_lin <- NULL
    myscores_mat$iHet_excl <- myscores_mat$iHet_rev <- NULL
    myscores_mat$response <- myscores_mat$dataset <- NULL
    myscores_mat <- myscores_mat %>%
      tidyr::pivot_wider(names_from = run, values_from = iHet) %>%
      as.data.frame()
    rownames(myscores_mat) <- myscores_mat$patient
    myscores_mat$patient <- NULL

    pred <- prediction(myscores_mat, myresponse_mat, label.ordering = c("NR", "R"))
    AUC_iHet <- performance(pred, measure = "auc")
    AUC_iHet <- unlist(AUC_iHet@y.values)
    AUC_iHet_CI <- qt(0.975, 100) * sqrt(var(AUC_iHet)) / sqrt(100) # 95% CI
    #cat("Average (bootstrap) AUC_iHet: ", AUC_iHet, "\n")

  }else{ ## median

    # iHet_Excl
    pred <- prediction(myscores_df$iHet_excl, myscores_df$response, label.ordering = c("NR", "R"))
    AUC_iHet_excl <- performance(pred, measure = "auc")
    AUC_iHet_excl <- AUC_iHet_excl@y.values[[1]]
    #cat("AUC_iHet_excl: ", AUC_iHet_excl, "\n")

    # iHet
    pred <- prediction(myscores_df$iHet, myscores_df$response, label.ordering = c("NR", "R"))
    AUC_iHet <- performance(pred, measure = "auc")
    AUC_iHet <- AUC_iHet@y.values[[1]]
    #cat("AUC_iHet: ", AUC_iHet, "\n")
  }

  # Integration of - iHet_excl, iHet and TMB

  if (with_TMB) {

    if(length(myscores_df$patient) > length(unique(myscores_df$patient))){ ## bootstrap
      data.tri <- do.call(rbind, lapply(1:100, function(bootstrap_sample){
        myscores_df_run <- subset(myscores_df, run == bootstrap_sample)
        data.tri <- do.call(rbind, lapply(grid, function(i.iHet){
          #cat("iHet: ", i.iHet, "\n")
          do.call(rbind, lapply(grid, function(i.excl){
            #cat("iHet_excluded: ", i.excl, "\n")
            do.call(rbind, lapply(grid, function(i.TMB){
              #cat("TMB: ", i.TMB, "\n")
              if (abs(1-(i.iHet+i.excl+i.TMB)) < 0.0001){

                ## duo
                pred_averaged_tmp <- i.iHet * myscores_df_run$iHet +
                  i.excl * myscores_df_run$iHet_excl

                if (normalize == "01_after"){
                  # avg duo to linear
                  if (all(pred_averaged_tmp == 0)) {
                    pred_averaged_lin <- 0
                  }else{
                    pred_averaged_lin <- linear_func(pred_averaged_tmp)
                  }
                }else{
                  pred_averaged_lin <- pred_averaged_tmp
                }

                ## trio
                pred_averaged = pred_averaged_lin + i.TMB*myscores_df_run$TMB_lin

                pred <- prediction(pred_averaged, myscores_df_run$response, label.ordering = c("NR", "R"))
                AUC_averaged <- performance(pred, measure = "auc")
                AUC_averaged <- AUC_averaged@y.values[[1]]
                return(data.frame(Top=i.TMB, Right=i.iHet, Left=i.excl, value=AUC_averaged,  run = bootstrap_sample))
              }else{
                return(NULL)
              }
            }))
          }))
        }))
      }))

      data.tri <- data.tri %>%
        dplyr::group_by(Top, Right, Left) %>%
        dplyr::summarise(mean_value = mean(value),
                         stdev_value = sqrt(var(value)),
                         n_value = length(value))

      data.tri$ciw_value <- qt(0.975, data.tri$n_value) * data.tri$stdev_value / sqrt(data.tri$n_value) # 95% CI
      data.tri$is.max <- FALSE
      data.tri$is.max[data.tri$mean_value == (max(data.tri[,"mean_value"]))] <- TRUE
      #cat("Average (bootstrap) AUC_max: ", max(data.tri[,"mean_value"]), "\n")

      # Add variable for the use of equal weights
      data.tri$equal.weight <- FALSE
      data.tri$equal.weight[data.tri$Right == 1/3 & data.tri$Right == 1/3 & data.tri$Left == 1/3] <- TRUE
      #cat("Average (bootstrap) AUC_equal.weight: ", data.tri$mean_value[data.tri$equal.weight], "\n")

    }else{ ## median
      data.tri <- do.call(rbind, lapply(grid, function(i.iHet){
        #cat("iHet: ", i.iHet, "\n")
        do.call(rbind, lapply(grid, function(i.excl){
          #cat("iHet_excluded: ", i.excl, "\n")
          do.call(rbind, lapply(grid, function(i.TMB){
            #cat("TMB: ", i.TMB, "\n")
            if (abs(1-(i.iHet+i.excl+i.TMB)) < 0.0001){

              ## duo
              pred_averaged_tmp = i.iHet * myscores_df$iHet +
                i.excl * myscores_df$iHet_excl

              if (normalize == "01_after"){
                # avg duo to linear
                if (all(pred_averaged_tmp == 0)) {
                  pred_averaged_lin <- 0
                }else{
                  pred_averaged_lin <- linear_func(pred_averaged_tmp)
                }
              }else{
                pred_averaged_lin <- pred_averaged_tmp
              }

              ## trio
              pred_averaged = pred_averaged_lin + i.TMB * myscores_df$TMB_lin

              pred <- prediction(pred_averaged, myscores_df$response, label.ordering = c("NR", "R"))
              AUC_averaged <- performance(pred, measure = "auc")
              AUC_averaged <- AUC_averaged@y.values[[1]]
              return(data.frame(Top=i.TMB, Right=i.iHet, Left=i.excl, value=AUC_averaged))
            }else{
              return(NULL)
            }
          }))
        }))
      }))

      data.tri$is.max <- FALSE
      data.tri$is.max[data.tri$value == (max(data.tri[,"value"]))] <- TRUE
      #cat("AUC_max: ", max(data.tri[,"value"]), "\n")

      # Add variable for the use of equal weights
      data.tri$equal.weight <- FALSE
      data.tri$equal.weight[data.tri$Right == 1/3 & data.tri$Right == 1/3 & data.tri$Left == 1/3] <- TRUE
      #cat("AUC_equal.weight: ", data.tri$value[data.tri$equal.weight], "\n")

    }

    values <- subset(data.tri, is.max == TRUE)[,1:3]
    palette <- c( "#FF9933", "#002C54", "#3375B2", "#CCDDEC", "#BFBFBF", "#000000")
    max_value <- subset(data.tri, is.max == TRUE)[,4:ncol(data.tri)]
    equal_value <- subset(data.tri, equal.weight == TRUE)[,4:ncol(data.tri)]

    if(length(myscores_df$patient) > length(unique(myscores_df$patient))){ ## bootstrap

      data.tri$check <- FALSE


      g <- ggtern(data=data.tri, ggtern::aes(x = Left, y = Top, z = Right)) +
        ggplot2::geom_point(shape=21,size=5, aes(fill=mean_value, col=is.max, alpha=is.max)) +
        ggtern::geom_crosshair_tern(data = subset(data.tri, is.max==TRUE), colour = "black") +
        ggtern::geom_crosshair_tern(data = subset(data.tri, equal.weight==TRUE), lty=2, colour = "black") +
        # ggtern::geom_crosshair_tern(data = subset(data.tri, Top == 1 & Right == 0 & Left == 0), colour = "black") +
        # ggtern::geom_crosshair_tern(data = subset(data.tri, Top == 0 & Right == 1 & Left == 0), colour = "black") +
        # ggtern::geom_crosshair_tern(data = subset(data.tri, Top == 0 & Right == 0 & Left == 1), colour = "black") +
        ggplot2::theme_bw() +
        ggtern::theme_hidetitles() +
        ggtern::theme_hidegrid() +
        ggtern::theme_showarrows() +
        # scale_fill_gradient(name="AUC", low="white", high=palette[2], limit=c(0.5, 0.9))+
        # scale_fill_gradient(name="AUC", trans = "log", low="#264653", high="#e76f51", limit=c(0.5, 0.9))+
        ggplot2::scale_alpha_manual(values = c(0.6,0.9))+
        viridis::scale_fill_viridis(name="AUC", option = "B", limit=c(0.5, 0.9)) +
        ggplot2::scale_colour_manual(values = c("white", "black")) +
        ggtern::scale_T_continuous(limits=c(0,1.0), breaks=seq(0,1,by=0.2), labels = seq(0,1,by=0.2)) +
        ggtern::scale_L_continuous(limits=c(0,1.0), breaks=seq(0,1,by=0.2), labels = seq(0,1,by=0.2)) +
        ggtern::scale_R_continuous(limits=c(0,1.0), breaks=seq(0,1,by=0.2), labels = seq(0,1,by=0.2)) +
        ggplot2::ggtitle(paste0(dataset_name,
                       "\nAUC(iHet) = ", paste0(round(mean(AUC_iHet),3),
                                                " (95% CI: ",
                                                round(mean(AUC_iHet)-AUC_iHet_CI,3), ", ",
                                                round(mean(AUC_iHet)+AUC_iHet_CI,3), ")"),
                       "\nAUC(iHet_Excl) = ", paste0(round(AUC_iHet_excl,3),
                                                     " (95% CI: ",
                                                     round(mean(AUC_iHet_excl)-AUC_iHet_excl_CI,3), ", ",
                                                     round(mean(AUC_iHet_excl)+AUC_iHet_excl_CI,3), ")"),
                       "\nAUC(TMB) = ", round(AUC_TMB,2),
                       "\nAUC(max) = ", paste0(round(max_value[,"mean_value"],2),
                                          " (95% CI: ",
                                          round(max_value[,"mean_value"]-max_value[,"ciw_value"],3), ", ",
                                          round(max_value[,"mean_value"]+max_value[,"ciw_value"],3), ")"),
                       "\nAUC(equal) = ", paste0(round(equal_value[,"mean_value"],2),
                                          " (95% CI: ",
                                          round(equal_value[,"mean_value"]-equal_value[,"ciw_value"],3), ", ",
                                          round(equal_value[,"mean_value"]+equal_value[,"ciw_value"],3), ")"),
                       "\n"))
                       #"\nweights:","\nTMB=", as.character(round(values[1],2)),
                       #",\niHet=", as.character(round(values[2],2)),
                       #",\niHet_excl=", as.character(round(values[3],2))))
      g = g +
        Tlab("TMB weight") +
        Rlab("iHet weight") +
        Llab("iHet_Excl weight") +
        ggplot2::theme(plot.title=element_text(vjust = -30),
                       axis.text = element_text(colour = "black", face = "bold", size = 12),
                       axis.title = element_text(face = "bold"))

    }else{ ## median

      g <- ggtern(data=data.tri,ggtern::aes(x = Left, y = Top, z = Right)) +
        ggplot2::geom_point(shape=21,size=5, aes(fill=value, col=is.max, alpha=is.max)) +
        ggtern::geom_crosshair_tern(data = subset(data.tri, is.max==TRUE), colour = "black") +
        ggplot2::theme_bw() +
        ggtern::theme_hidetitles() +
        ggtern::theme_hidegrid() +
        ggtern::theme_showarrows() +
        # scale_fill_gradient(name="AUC", low="white", high=palette[2], limit=c(0.5, 0.9))+
        # scale_fill_gradient(name="AUC", trans = "log", low="#264653", high="#e76f51", limit=c(0.5, 0.9))+
        ggplot2::scale_alpha_manual(values = c(0.6,0.9))+
        viridis::scale_fill_viridis(name="AUC", option = "B", limit=c(0, 1))+
        ggplot2::scale_colour_manual(values = c("white", "black")) +
        ggplot2::ggtitle(paste0(dataset_name,
                       "\niHet: ", round(AUC_iHet,2),
                       "\niHet_excl): ", round(AUC_iHet_excl,2),
                       "\nTMB: ", round(AUC_TMB,2),
                       "\nmax: ", round(max(data.tri[,"value"]),2), "\n",
                       "\nweights:","\nTMB=", as.character(round(values[1],2)),
                       ",\niHet=", as.character(round(values[2],2)),
                       ",\niHet_excl)=", as.character(round(values[3],2))))
      g = g +
        Tlab("TMB") +
        Rlab("iHet") +
        Llab("iHet_excl") +
        ggplot2::theme(plot.title=element_text(vjust = -30),
                       axis.text = element_text(colour = "black", face = "bold", size = 12),
                       axis.title = element_text(face = "bold"))
    }
    return(g)
    #return(list("data.duo"=data.duo, "data.duo_tmb"=data.duo_tmb,"data.tri"=data.tri))
  }else{
    return(list("data.duo"=data.duo))
  }
}

#' Apply hill function to numeric vector
hill_func <- function(x) {
  min_x <- min(x)
  max_x <- max(x)
  x01 = (x - min_x) / (max_x - min_x)
  k=0.5
  n=5
  y <- (k^n + 1)*(x01^n)/(k^n+x01^n)
  return(y)
}

#' Apply linear function to numeric vector
linear_func <- function(x, range = c("01", "11")) {
  min_x <- min(x)
  max_x <- max(x)
  x01 = (x - min_x) / (max_x - min_x)
  x11 = (2 * x01) - 1
  if (range == "01") return(x01)
  if (range == "11") return(x11)
}


loadRData <- function(fileName){ #to assign the loaded variable to the desired variable name
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
