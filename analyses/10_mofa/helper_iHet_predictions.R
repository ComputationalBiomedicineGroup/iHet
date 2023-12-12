
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
  library('dorothea')
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
  conflict_prefer("theme_bw", "ggplot2")
  library('ggforce')
  library('svglite')
  library('tidyr')
  library('plyr')
  library('purrr')
  library('survival')
  library('ggfortify')
  library('survminer')
  library('gtsummary')
  library('stringr')
  library('rstatix')
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
  gene_expr <- gene_expr[-which(genes_info$new_names == "entry withdrawn"), ]
  genes_info <- genes_info[-which(genes_info$new_names == "entry withdrawn"), ]
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

#' Import validation datasets
import_dataset <- function(dataset) {
  if (dataset %in% c("poplar", "oak")){
    myfilename <-list.files(path = paste0(folderValidation, "poplaroak"))
  }else{
    myfilename <-list.files(path = paste0(folderValidation, dataset))
  }
  myfilename <- myfilename[grep(".RData", myfilename)]
  myfilename <- myfilename[grep("all", myfilename)]
  
  #myfilename <-list.files(path = folderValidation)
  #myfilename <- myfilename[grep(dataset, myfilename, ignore.case=T)]
  
  if (length(grep("ontreatment", myfilename))>0){
    myfilename <- myfilename[-grep("ontreatment", myfilename)]
  }
  
  #all_dataset <- readRDS(paste0(folderValidation, myfilename))
  
  if (dataset %in% c("poplar", "oak")){
    all_dataset <- loadRData(paste0(folderValidation, "poplaroak", "/", myfilename))
    
    patients <- rownames(all_dataset$clinical_response %>% filter(STUDYNAME == toupper(dataset)))
    all_dataset$tpm <- all_dataset$tpm[, patients]
    all_dataset$counts <- all_dataset$counts[, patients]
    all_dataset$clinical_response <- all_dataset$clinical_response[patients, ]
    all_dataset$TMB <- all_dataset$TMB[patients, , drop=FALSE]
    all_dataset$pred_easier <- all_dataset$pred_easier[patients, ]
    all_dataset$pred_tasks <- all_dataset$pred_tasks[patients, ]
    all_dataset$computed_features <- lapply(all_dataset$computed_features, function(X) X[patients, ])
    
  }else{
    all_dataset <- loadRData(paste0(folderValidation, dataset, "/", myfilename))
    
    prediction_easier <- apply(all_dataset$pred_easier,1,median)
    prediction_tasks <- apply(scale(all_dataset$pred_tasks),1,median)
    #prediction_tasks <- all_dataset$pred_tasks$response
    names(prediction_tasks) <- rownames(all_dataset$pred_tasks)
    
    TGFb <- all_dataset$computed_features[["pathway_activity"]][,"TGFb"]
    names(TGFb) <- rownames(all_dataset$computed_features[["pathway_activity"]])
    
    EGFR <- all_dataset$computed_features[["pathway_activity"]][,"EGFR"]
    names(EGFR) <- rownames(all_dataset$computed_features[["pathway_activity"]])
  }
  
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

# categorize_TMB <- function(x, thresholds=NULL) {
#   vTert = quantile(x , c(0:3/3))
#   if (is.null(thresholds)){
#     x.out = cut(x, vTert, include.lowest = T, labels = c(1, 2, 3))
#   }else if ((is.numeric(thresholds)) & (length(thresholds)==2)){
#     # x.out = cut(x, c(min(x, na.rm = T), 100, 400, max(x, na.rm = T)), include.lowest = T, labels = c(1, 2, 3))
#     x.out = cut(x, c(min(x, na.rm = T), 100, 400, max(x, na.rm = T)), include.lowest = F, labels = c(1, 2, 3))
#   }
#   x.out <- as.numeric(x.out)
#   names(x.out) <- names(x)
#   return(x.out)
# }

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

get_iHet_associated_weights <- function(use_bootstrap_weights=FALSE, 
                                        model="NSCLC",
                                        use_fibroblasts_hifs_tissue = c("ESI", "tumor", "stroma", "epithelial"),
                                        cancer_type_hifs = c("NSCLC", "pancancer"),
                                        scale_mofa_weights = FALSE){
  
  # Load knowledge of features correlated with information about the density of fibroblasts in tumor region 
  feat_cor <- readRDS(corFibAllTissues)
  feat_cor_tmp <- feat_cor[[use_fibroblasts_hifs_tissue]][[cancer_type_hifs]]
  rownames(feat_cor_tmp) <- feat_cor_tmp$feature
  feat_cor_tmp <- feat_cor_tmp[!feat_cor_tmp$feature %in% c("Other", "DC", "Monocyte"),]
  feat_sign <- feat_cor_tmp$feature[which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<=0.05)]

  # 1. get MOFA weights
  # 2. Adjust sign of the factor 1 to have consistent positive 
  # association with immune response across bootstrap runs (get_rotation())
  
  ## Bootstrap
  if (use_bootstrap_weights){
    which_files <- grep(model, utils::unzip(zipfile=zipMOFAmodels, list=TRUE)[,1])
    files_to_extract <- utils::unzip(zipfile=zipMOFAmodels, list=TRUE)[which_files,1]
    files_to_extract <- files_to_extract[grep("hdf5", files_to_extract)]
    
    all_weights_run <- lapply(1:length(files_to_extract), function(run){
      # not load following 1:100 order
      model_bootstrap_run <- MOFA2::load_model(utils::unzip(zipfile=zipMOFAmodels, 
                                                            files=files_to_extract[run]))
      weights_bootstrap_run <- MOFA2::get_weights(model_bootstrap_run, scale = scale_mofa_weights)
      return(weights_bootstrap_run)
    })
    
    F1_weights_df <- lapply(1:length(all_weights_run), function(dataset_id) {
      weights <- all_weights_run[[dataset_id]]
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
      bind_rows() %>% select("feature", "factor 1", "view", "run")
    
    colnames(F1_weights_df)[2] <- "weight"
    
    ## Median  
  }else{
    ## FF provided file: (paste0(repo_dir, "/data/iHet_data/median_weights.rds"))
    ## r=0.92 between these computed weights and FF's
    median_weights <- lapply(1:length(all_weights_run), function(dataset_id) {
      weights <- all_weights_run[[dataset_id]]
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
    
    F1_weights_df <- median_weights %>% select(view, feature, "factor 1")
    
  }
  colnames(F1_weights_df)[3] <- "weight"
  
  # iHet_excl weights (set cor features to 0)
  weights_excl <- F1_weights_df %>%
    filter(feature %in% feat_sign) %>%
    mutate(weight=0, score = "iHet_excl")
  
  # iHet_rev (negate only cor features)
  weights_rev <- F1_weights_df %>%
    filter(feature %in% feat_sign) %>%
    mutate(weight=-weight, score = "iHet_rev")
  
  # iHet
  weights_iHet <- F1_weights_df %>%
    mutate(score = "iHet")
  
  all_weights <- rbind(weights_iHet, weights_excl, weights_rev)
  all_weights_df <- all_weights %>% pivot_wider(id_cols = c(feature, view), names_from = score, values_from = weight)
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
                         method=c("iHet", "iHet_removed", "iHet_inverted", "iHet_exclusion"),
                         use_fibroblasts_hifs_tissue = c("ESI", "tumor", "stroma", "epithelial"),
                         cancer_type_hifs = c("NSCLC", "pancancer"),
                         assess_cor_threshold = FALSE,
                         scale_mofa_weights = FALSE){
  
  # match arguments
  model <- match.arg(model, c("BLCA", "BRCA", "CESC","CRC", "GBM", "HNSC", "KIRC", "KIRP", 
                              "LIHC", "LUAD", "LUSC", "NSCLC", "OV", "PAAD", "PRAD",
                              "SKCM", "STAD", "THCA", "UCEC", "JiaSharma"))
  
  # use_fibroblasts_hifs_tissue <- match.arg(use_fibroblasts_hifs_tissue, c("ESI", "tumor", "stroma", "epithelial"))
  
  # 1. get MOFA weights
  # 2. Adjust sign of the factor 1 to have consistent positive 
  # association with immune response across bootstrap runs (get_rotation())
  
  ## Bootstrap
  if (use_bootstrap_weights){
    which_files <- grep(model, utils::unzip(zipfile=zipMOFAmodels, list=TRUE)[,1])
    files_to_extract <- utils::unzip(zipfile=zipMOFAmodels, list=TRUE)[which_files,1]
    files_to_extract <- files_to_extract[grep("hdf5", files_to_extract)]
    
    all_weights_run <- lapply(1:length(files_to_extract), function(run){
      # not load following 1:100 order
      model_bootstrap_run <- MOFA2::load_model(utils::unzip(zipfile=zipMOFAmodels, 
                                                            files=files_to_extract[run]))
      weights_bootstrap_run <- MOFA2::get_weights(model_bootstrap_run, scale = scale_mofa_weights)
      return(weights_bootstrap_run)
    })

    F1_weights_df <- lapply(1:length(all_weights_run), function(dataset_id) {
      weights <- all_weights_run[[dataset_id]]
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
      bind_rows() %>% select("feature", "factor 1", "view", "run")
  
    colnames(F1_weights_df)[2] <- "weight"

  ## Median  
  }else{
    ## FF provided file: (paste0(repo_dir, "/data/iHet_data/median_weights.rds"))
    ## r=0.92 between these computed weights and FF's
    median_weights <- lapply(1:length(all_weights_run), function(dataset_id) {
      weights <- all_weights_run[[dataset_id]]
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
  
  # immune cell fractions
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
  cellfrac <- log10(cellfrac*100+0.001)
  
  # ---
  # TF activity (using last version from dorothea)
  if (!any(packageVersion("dorothea") > '1.6.0')) stop("dorothea version 1.6.0 or 1.7.0 is required")
  tf <- easier::compute_TF_activity(RNA_tpm = dataset$tpm)
  
  # ---
  # pathway activity
  path <- easier::compute_pathway_activity(RNA_counts = dataset$counts,
                                           remove_sig_genes_immune_response = FALSE)
  # ---
  # scale features according to TCGA feature values from MOFA training
  boot_features_stats <- readRDS(TCGAbootstrapFeatureStats)
  boot_features_stats <- boot_features_stats[[model]]
  
  if (model == "NSCLC"){
    
    norm_boot_features <- lapply(c("LUAD", "LUSC"), function(XX){
      norm_boot_features <- lapply(1:length(boot_features_stats), function(run){
        mean_cellfrac <- boot_features_stats[[run]][[XX]]$cellfrac$mean[colnames(cellfrac)]
        cellfrac_centered <- sweep(cellfrac, 2, mean_cellfrac, FUN = "-")
        mean_tf <- boot_features_stats[[run]][[XX]]$tf$mean[colnames(tf)]
        tf_centered <- sweep(tf, 2, mean_tf, FUN = "-")
        path_centered <- sweep(path, 2, boot_features_stats[[run]][[XX]]$path$mean, FUN = "-")
        path_scaled <- sweep(path_centered, 2, boot_features_stats[[run]][[XX]]$path$sd, FUN = "/")
        return(list(path = as.matrix(path_scaled), cellfrac = as.matrix(cellfrac_centered), tf = as.matrix(tf_centered)))
      })
      return(norm_boot_features)
    })
    names(norm_boot_features) <- c("LUAD", "LUSC")
    
  }else{
    norm_boot_features <- lapply(1:length(boot_features_stats), function(run){
      mean_cellfrac <- boot_features_stats[[run]]$cellfrac$mean[colnames(cellfrac)]
      cellfrac_centered <- sweep(cellfrac, 2, mean_cellfrac, FUN = "-")
      mean_tf <- boot_features_stats[[run]]$tf$mean[colnames(tf)]
      tf_centered <- sweep(tf, 2,mean_tf, FUN = "-")
      path_centered <- sweep(path, 2, boot_features_stats[[run]]$path$mean, FUN = "-")
      path_scaled <- sweep(path_centered, 2, boot_features_stats[[run]]$path$sd, FUN = "/")
      return(list(path = as.matrix(path_scaled), cellfrac= as.matrix(cellfrac_centered), tf = as.matrix(tf_centered)))
    })
  }
  
  if (use_bootstrap_weights == FALSE){
    if (model == "NSCLC"){
      norm_boot_features <- reshape2::melt(norm_boot_features)
      colnames(norm_boot_features) <- c("patient", "feature", "value", "view", "run", "model")
      norm_boot_features <- norm_boot_features %>% 
        group_by(patient, feature, view, model) %>%
        dplyr::summarise(median_value = median(value)
        )
      
    }else{
      norm_boot_features <- reshape2::melt(norm_boot_features)
      colnames(norm_boot_features) <- c("patient", "feature", "value", "view", "run")
      norm_boot_features <- norm_boot_features %>% 
        group_by(patient, feature, view) %>%
        dplyr::summarise(median_value = median(value)
        )
    }
  }
  
  # Load knowledge of features correlated with information about the density of fibroblasts in tumor region 
  feat_cor <- readRDS(corFibAllTissues)

  cat("Computing iHet score... \n")
  if (use_bootstrap_weights){
    
    # Just for NSCLC mofa model
    if (model == "NSCLC"){
      
      iHet_bootstrap <- lapply(c("LUAD", "LUSC"), function(XX){
        
        # compute iHet score
        iHet_bootstrap <- lapply(1:length(files_to_extract), function(bot_sample){
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
              
              iHet_hifs <- lapply(use_fibroblasts_hifs_tissue, function(tissue){
                print(tissue)
                feat_cor_tmp <- feat_cor[[tissue]][[cancer_type_hifs]]
                rownames(feat_cor_tmp) <- feat_cor_tmp$feature
                feat_cor_tmp <- feat_cor_tmp[cfeatures,]
                
                # use correlation different thresholds
                if(assess_cor_threshold){
                  seq_cor <- seq(0,.3,.1)
                  iHet <- vector("list", length = length(seq_cor))
                  F1 <- vector("list", length = length(seq_cor))
                  ix_sign <- vector("list", length = length(seq_cor))
                  names(ix_sign) <- names(F1) <- names(iHet) <- seq_cor
                  
                  # Compute iHet for different corr thresholds
                  for(thr in seq_cor){
                    cat("Threshold", thr, ": \n")
                    F1[[as.character(thr)]] <- F1_vec
                    if (any(feat_cor_tmp$cor>=thr)){
                      ix_sign[[as.character(thr)]] <- which(feat_cor_tmp$cor>thr & feat_cor_tmp$FDR<=0.05)
                      
                      if (method_x == "iHet_exclusion") {
                        iHet[[as.character(thr)]] <- F1[[as.character(thr)]][ix_sign[[as.character(thr)]]] %*% features[ix_sign[[as.character(thr)]], , drop=FALSE]
                        
                      } else if (method_x == "iHet_inverted") {
                        F1[[as.character(thr)]][ix_sign[[as.character(thr)]]] <- - F1[[as.character(thr)]][ix_sign[[as.character(thr)]]]
                        iHet[[as.character(thr)]] <- F1[[as.character(thr)]] %*% features 
                        
                      } else if (method_x == "iHet_removed") {
                        F1[[as.character(thr)]][ix_sign[[as.character(thr)]]]<- 0
                        iHet[[as.character(thr)]] <- F1[[as.character(thr)]] %*% features 
                      }
                      
                      iHet[[as.character(thr)]] <- t(iHet[[as.character(thr)]][1,,drop=FALSE])
                      
                    }else{
                      ix_sign[[as.character(thr)]] <- NA
                      F1[[as.character(thr)]] <- NA
                      iHet[[as.character(thr)]] <- NA
                    }
                  }
                  return(iHet)
                  # assess_cor_threshold = FALSE, the p-value is used as threshold
                }else{
                  ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<=0.05)
                  #ix_sign <- which(feat_cor$cor>=0.2 & feat_cor$FDR<0.05)
                  
                  if (method_x == "iHet_exclusion") {
                    iHet <- F1[ix_sign] %*% features[ix_sign,]
                  } else if (method_x == "iHet_inverted") {
                    F1_vec[ix_sign] <- - F1_vec[ix_sign]
                    iHet <- F1_vec %*% features
                  } else if (method_x == "iHet_removed") {
                    F1_vec[ix_sign] <- 0
                    iHet <- F1_vec %*% features
                  }
                  
                  iHet <- t(iHet[1,,drop=FALSE])
                  return(iHet)
                }
              })
              names(iHet_hifs) <- use_fibroblasts_hifs_tissue
              return(iHet_hifs)
              # just iHet method
            }else{
              iHet <- F1_vec %*% features
              #iHet <- qr.solve(F1_vec, features)
              iHet <- t(iHet[1,,drop=FALSE])
              return(iHet)
            }
          })
          names(iHet_method) <- method
          return(iHet_method)
        })
        names(iHet_bootstrap) <- 1:length(files_to_extract)
        return(iHet_bootstrap)
      })
      names(iHet_bootstrap) <- c("LUAD", "LUSC")
      
      # For the rest of mofa models
    }else{
      
      # compute iHet score
      iHet_bootstrap <- lapply(1:length(files_to_extract), function(bot_sample){
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
            
            iHet_hifs <- lapply(use_fibroblasts_hifs_tissue, function(tissue){
              print(tissue)
              feat_cor_tmp <- feat_cor[[tissue]][[cancer_type_hifs]]
              rownames(feat_cor_tmp) <- feat_cor_tmp$feature
              feat_cor_tmp <- feat_cor_tmp[cfeatures,]
              
              # use correlation different thresholds
              if(assess_cor_threshold){
                seq_cor <- seq(0,.3,.1)
                iHet <- vector("list", length = length(seq_cor))
                iHet_exclusion <- vector("list", length = length(seq_cor))
                F1 <- vector("list", length = length(seq_cor))
                ix_sign <- vector("list", length = length(seq_cor))
                names(ix_sign) <- names(F1) <- names(iHet) <- names(iHet_exclusion) <- seq_cor
                
                # Compute iHet for different corr thresholds
                for(thr in seq_cor){
                  cat("Threshold", thr, ": \n")
                  F1[[as.character(thr)]] <- F1_vec
                  
                  if (any(feat_cor_tmp$cor>=thr)){
                      ix_sign[[as.character(thr)]] <- which(feat_cor_tmp$cor>thr & feat_cor_tmp$FDR<=0.05)
    
                    if (method_x == "iHet_exclusion") {
                      iHet[[as.character(thr)]] <- F1[[as.character(thr)]][ix_sign[[as.character(thr)]]] %*% features[ix_sign[[as.character(thr)]], , drop=FALSE]

                    } else if (method_x == "iHet_inverted") {
                      F1[[as.character(thr)]][ix_sign[[as.character(thr)]]] <- - F1[[as.character(thr)]][ix_sign[[as.character(thr)]]]
                      iHet[[as.character(thr)]] <- F1[[as.character(thr)]] %*% features 
                      
                    } else if (method_x == "iHet_removed") {
                      F1[[as.character(thr)]][ix_sign[[as.character(thr)]]]<- 0
                      iHet[[as.character(thr)]] <- F1[[as.character(thr)]] %*% features 
                    }
                    
                    iHet[[as.character(thr)]] <- t(iHet[[as.character(thr)]][1,,drop=FALSE])
                    
                  }else{
                    ix_sign[[as.character(thr)]] <- NA
                    F1[[as.character(thr)]] <- NA
                    iHet[[as.character(thr)]] <- NA
                  }
                }
                return(iHet)
                
                # assess_cor_threshold = FALSE, the p-value is used as threshold
              } else{
                
                ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<=0.05)
                #ix_sign <- which(feat_cor$cor>0.2 & feat_cor$FDR<=0.05)
                
                if (method_x == "iHet_exclusion") {
                  iHet <- F1[ix_sign] %*% features[ix_sign, , drop=FALSE]
                } else if (method_x == "iHet_inverted") {
                  F1_vec[ix_sign] <- - F1_vec[ix_sign]
                  iHet <- F1_vec %*% features
                } else if (method_x == "iHet_removed") {
                  F1_vec[ix_sign] <- 0
                  iHet <- F1_vec %*% features
                }
                
                iHet <- t(iHet[1,,drop=FALSE])
                return(iHet)
              }
            })
            names(iHet_hifs) <- use_fibroblasts_hifs_tissue
            return(iHet_hifs)
            # just iHet method
          }else{
            iHet <- F1_vec %*% features
            #iHet <- qr.solve(F1_vec, features)
            iHet <- t(iHet[1,,drop=FALSE])
            return(iHet)
          }
        })
        names(iHet_method) <- method
        return(iHet_method)
      })
      names(iHet_bootstrap) <- 1:length(files_to_extract)
    }
    
    # to data.frame
    iHet_bootstrap_df <- reshape2::melt(iHet_bootstrap)
    iHet_bootstrap_df$Var2 <- NULL
    if(assess_cor_threshold) {
      if ("iHet" %in% method){
        iHet_bootstrap_df[is.na(iHet_bootstrap_df[,3]), 3] <- "original"
        iHet_bootstrap_df[is.na(iHet_bootstrap_df[,4]), 4] <- "original"
      }
      if (model == "NSCLC"){
        colnames(iHet_bootstrap_df) <- c("patient", "score", "threshold", "hif_tissue", "method", "run", "model")
      }else{
        colnames(iHet_bootstrap_df) <- c("patient", "score", "threshold", "hif_tissue", "method", "run")
      }
    }else{
      if(length(method)==1 & "iHet" %in% method & model == "NSCLC"){
        colnames(iHet_bootstrap_df) <- c("patient", "score", "method", "run", "model")
      }
      if(length(method)>1 & "iHet" %in% method){
        iHet_bootstrap_df[is.na(iHet_bootstrap_df[,3]), 3] <- "original"
      }
      if (length(method)>1 & model == "NSCLC"){
        colnames(iHet_bootstrap_df) <- c("patient", "score", "hif_tissue", "method", "run", "model")
      }else if (length(method)>1 & model != "NSCLC"){
        iHet_bootstrap_df$patient <- names(iHet_bootstrap[[1]][[1]])
        colnames(iHet_bootstrap_df) <- c("patient", "score", "hif_tissue", "method", "run")
      }
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
          filter(view=="path"& model == XX) %>% 
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
            rownames(feat_cor_tmp) <- feat_cor_tmp$feature
            feat_cor_tmp <- feat_cor_tmp[cfeatures,]
            ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<0.05)
            if (method_x == "iHet_inverted") {
              F1_vec[ix_sign] <- - F1_vec[ix_sign]
            } else if (method_x == "iHet_removed") {
              F1_vec[ix_sign] <- 0
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
        cfeatures <- intersect(colnames(features), names(F1))
        features <- features[, cfeatures]
        F1_vec <- F1[cfeatures]
        features <- t(features)
        
        iHet_method <- lapply(method, function(method_x){
          # modify the F1 based on the defined method
          if (method_x != "iHet") {
            rownames(feat_cor_tmp) <- feat_cor_tmp$feature
            feat_cor_tmp <- feat_cor_tmp[cfeatures,]
            ix_sign <- which(feat_cor_tmp$cor>0 & feat_cor_tmp$FDR<0.05)
            if (method_x == "iHet_inverted") {
              F1_vec[ix_sign] <- - F1_vec[ix_sign]
            } else if (method_x == "iHet_removed") {
              F1_vec[ix_sign] <- 0
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
                                          use_fibroblasts_hifs_tissue,
                                          assess_cor_threshold=FALSE, 
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
  col <- c("#7570B3","#D95F02")
  names(col) <- c("iHet", "iHet_inverted")
  col <- col[c("iHet", "iHet_inverted")]
  
  #col=colorRampPalette(brewer.pal(8, "Set2"))(6)
  #names(col) <- c("IR", "iHet", "iHet_inverted", "iHet_removed", "iHet_exclusion")
  #col <- col[c("iHet", "iHet_inverted", "iHet_removed", "iHet_exclusion")]
  
  rp_df <- list()
  
  ## sort iHet (for correct transformation into matrix, prior to predict response) ##
  iHet <- iHet %>% group_by(patient) %>% arrange(patient)
  
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
    
    if(yy %in% c("iHet_inverted", "iHet_removed", "iHet_exclusion")){
      
      ROC_data <- lapply(use_fibroblasts_hifs_tissue, function(xx){ 
        
        iHet_sel <- subset(rp_df$iHet, method == yy & hif_tissue == xx)
        
        if(assess_cor_threshold){
          
          perf_thr <- lapply(unique(iHet_sel$threshold), function(zz){
            iHet_method_thr <- subset(iHet_sel, threshold == zz)
            iHet_mat <- matrix(iHet_method_thr$score, 
                               nrow=length(unique(iHet_method_thr$patient)),
                               ncol=length(unique(iHet_method_thr$run)),
                               byrow = TRUE)
            
            rownames(iHet_mat) <- unique(iHet_method_thr$patient)
            colnames(iHet_mat) <- unique(iHet_method_thr$run)
            
            # if (yy == "iHet_exclusion") pred <- ROCR::prediction(-iHet_mat, rp_df$response, label.ordering = c("NR", "R"))
            pred <- ROCR::prediction(iHet_mat, rp_df$response, label.ordering = c("NR", "R"))
            roc_perf <- ROCR::performance(pred,measure = yaxis, x.measure = xaxis)
            auc_perf <- unlist(ROCR::performance(pred, auc_metric)@y.values)
            auc_perf <- matrix(auc_perf); rownames(auc_perf) <- 1:100 
            return(list(auc_perf=auc_perf))
          })
          names(perf_thr) <- unique(iHet_sel$threshold)
          return(perf_thr)
          
        }else{
          
          # iHet_sel <- subset(iHet_sel, threshold == 0)
          iHet_mat <- matrix(iHet_sel$score, 
                             nrow=length(unique(iHet_sel$patient)),
                             ncol=length(unique(iHet_sel$run)),
                             byrow = TRUE)
          
          rownames(iHet_mat) <- unique(iHet_sel$patient)
          colnames(iHet_mat) <- unique(iHet_sel$run)
          
          pred <- ROCR::prediction(iHet_mat, rp_df$response, label.ordering = c("NR","R"))
          roc_perf <- ROCR::performance(pred, measure = yaxis, x.measure = xaxis)
          auc_perf <- unlist(ROCR::performance(pred, auc_metric)@y.values)
          auc_perf <- matrix(auc_perf); rownames(auc_perf) <- 1:100 
          return(list(roc_perf=roc_perf, auc_perf=auc_perf))
        }
        return(ROC_data)
      })
    }else{
      
      iHet_sel <- subset(rp_df$iHet, method == "iHet")
      iHet_mat <- matrix(iHet_sel$score, 
                         nrow=length(unique(iHet_sel$patient)),
                         ncol=length(unique(iHet_sel$run)),
                         byrow = TRUE)
      
      rownames(iHet_mat) <- unique(iHet_sel$patient)
      colnames(iHet_mat) <- unique(iHet_sel$run)
      pred <- ROCR::prediction(iHet_mat, rp_df$response, label.ordering = c("NR", "R"))
      roc_perf <- ROCR::performance(pred, measure = yaxis, x.measure = xaxis)
      auc_perf <- unlist(ROCR::performance(pred, auc_metric)@y.values)
      auc_perf <- matrix(auc_perf); rownames(auc_perf) <- 1:100 
      return(list(roc_perf=roc_perf, auc_perf=auc_perf))
    }
    names(ROC_data) <- use_fibroblasts_hifs_tissue
    return(ROC_data)
  })
  names(ROC_data) <- method
  
  if (plot_roc==TRUE){
    
    if(use_bootstrap_weights) {
      
      legend_text <- sapply(method, function(xx){
        # confidence intervals
        if(xx != "iHet") tmp <- as.vector(ROC_data[[xx]]$tumor$auc_perf)
        if(xx == "iHet") tmp <- ROC_data[[xx]]$auc_perf
        mean <- mean(tmp)
        stdev <-  sqrt(var(tmp))
        n     <- 100
        ciw   <- qt(0.975, n) * stdev / sqrt(n) # 95% CI
        
        # using ROCR package
        if(xx != "iHet") ROCR::plot(ROC_data[[xx]]$tumor$roc_perf, lwd=3, cex.lab = 1.3, col=col[xx], avg = "threshold", type = "S", add = TRUE)
        if(xx == "iHet") ROCR::plot(ROC_data[[xx]]$roc_perf, lwd=3, cex.lab = 1.3, col=col[xx], avg = "threshold", type = "S", add = TRUE)
        
        if(metric=="precision_recall"){
          legend_text <- paste0(xx,", AUCPR=", round(mean, 3), 
                                " (95% CI: ", round(mean-ciw,3), ", ", 
                                round(mean+ciw,3), ")")
        }else{
          legend_text <- paste0(xx,", AUC=", round(mean, 3), 
                                " (95% CI: ", round(mean-ciw,3), ", ", 
                                round(mean+ciw,3), ")")
        }
        return(legend_text)
      })
      
    }else{
      
      legend_text <- sapply(method, function(xx){
        if(xx != "iHet") tmp <- as.vector(ROC_data[[xx]]$tumor$auc_perf)
        if(xx == "iHet") tmp <- ROC_data[[xx]]$auc_perf
        
        # using ROCR package
        if(xx != "iHet") ROCR::plot(ROC_data[[xx]]$tumor$roc_perf, lwd=3, cex.lab = 1.3, col=col[xx], avg = "threshold", type = "S", add = TRUE)
        if(xx == "iHet") ROCR::plot(ROC_data[[xx]]$roc_perf, lwd=3, cex.lab = 1.3, col=col[xx], avg = "threshold", type = "S", add = TRUE)
        
        if(metric=="precision_recall"){
          legend_text <- paste0(xx,", AUCPR=", round(tmp, 3))
        }else{
          legend_text <- paste0(xx,", AUC=", round(tmp, 3))
        }
        return(legend_text)
      })
      
    }
    
    if(metric=="precision_recall"){
      legend(x=0.35, y=1.05, xjust = 0,
             legend = c(legend_text),
             fill = col[method], bty="n", cex=1)
    }else{
      legend(x=0.40, y=0.25, xjust = 0,
             legend = c(legend_text),
             fill = col[method], bty="n", cex = 0.7)
    }
    
    plot <- grDevices::recordPlot()
    
  }
  
  if(assess_cor_threshold){
    
    # collect AUC bootstrap values
    ROC_data$iHet$roc_perf <- NULL
    auc_df_iHet <- melt(ROC_data$iHet)
    ROC_data$iHet <- NULL
    
    auc_df <- melt(ROC_data)
    auc_df_iHet$Var2 <- auc_df$Var2 <- NULL
    auc_df_iHet[,3] <- auc_df[,3] <- NULL
    colnames(auc_df) <- c("run", "auc","threshold", "hif_tissue", "method")
    colnames(auc_df_iHet) <- c("run", "auc")
    auc_df_iHet$threshold <- "original"
    auc_df_iHet$hif_tissue <- "original"
    auc_df_iHet$method <- "iHet"
    auc_df <- rbind(auc_df_iHet, auc_df)
    
    AUC <- auc_df
    
  }else{
    
    # collect AUC bootstrap values
    ROC_data$iHet$roc_perf <- NULL
    auc_df_iHet <- melt(ROC_data$iHet)
    ROC_data$iHet <- NULL
    
    ROC_data$iHet_inverted$tumor$roc_perf <- NULL
    auc_df_iHetinverted <- melt(ROC_data)
    auc_df_iHet$Var2 <- auc_df_iHetinverted$Var2 <- NULL
    
    auc_df_iHet[,c(3,4)] <- NULL; auc_df_iHetinverted[,c(3,4)] <- NULL
    colnames(auc_df_iHet) <- c("run", "auc")
    auc_df_iHet$method <- "iHet"
    colnames(auc_df_iHetinverted) <- c("run", "auc", "method")
    
    auc_df <- rbind(auc_df_iHet, auc_df_iHetinverted)
    AUC <- auc_df
    
    # AUC <- NULL
    # AUC <- c(AUC, lapply(names(ROC_data), function(X) {
    #   if(X == "iHet") AUC <- ROC_data[[X]]$auc_perf
    #   if(X != "iHet") {
    #     AUC <- lapply(names(ROC_data[[2]]), function(Y){
    #       ROC_data[[X]][[Y]]$auc_perf
    #     })
    #   }
    #   return(AUC)
    # }))
    # names(AUC) <- names(ROC_data)
    
    if (plot_roc==TRUE) {
      return(plot)
    }else{
      return(AUC)
    }

  }
}

#' Integrate iHet, iHet_excl and TMB
integrated_score <- function(myscores_df, grid=seq(0, 1, length.out = 19),
                             with_TMB = FALSE, save_figure_format = c(".pdf", ".png"),
                             tissue, normalize = c("01_first", "01_after")) {

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
      cat("AUC_TMB: ", AUC_TMB, "\n")
    }
  }

  if (with_TMB) myscores_df$TMB_lin <- linear_func(myscores_df$TMB_cat, range = "01")
  
  if (normalize == "01_first"){
    if(length(myscores_df$patient) > length(unique(myscores_df$patient))){
      ## bootstrap
      myscores_df <- myscores_df %>% dplyr::group_by(run) %>%
        dplyr::mutate(iHet = linear_func(iHet, range = "01"),
                      iHet_exclusion = linear_func(iHet_exclusion, range = "01"))
      
    }else{
      ## median
      myscores_df$iHet <- linear_func(myscores_df$iHet, range = "01")
      myscores_df$iHet_exclusion <- linear_func(myscores_df$iHet_exclusion, range = "01")
    }

  }
  
  # compute predictions
  if(length(myscores_df$patient) > length(unique(myscores_df$patient))){
    ## bootstrap
    # response matrix
    myresponse_mat <- myscores_df
    if (with_TMB) myresponse_mat$TMB <- myresponse_mat$TMB_cat <- myresponse_mat$TMB_lin <- NULL
    myresponse_mat$iHet_removed <-  myresponse_mat$iHet_inverted <- myresponse_mat$iHet_exclusion <- myresponse_mat$iHet <- NULL
    myresponse_mat$dataset <- NULL
    myresponse_mat <- myresponse_mat %>% pivot_wider(names_from = run, values_from = response) %>% as.data.frame()
    rownames(myresponse_mat) <- myresponse_mat$patient
    myresponse_mat$patient <- NULL
    
    ### iHet_Exclusion
    # predictions matrix
    myscores_mat <- myscores_df
    if (with_TMB) myscores_mat$TMB <- myscores_mat$TMB_cat <- myscores_mat$TMB_lin <- NULL
    myscores_mat$iHet_removed <-  myscores_mat$iHet_inverted <- myscores_mat$iHet <- NULL
    myscores_mat$response <- myscores_mat$dataset <- NULL
    myscores_mat <- myscores_mat %>% pivot_wider(names_from = run, values_from = iHet_exclusion) %>% as.data.frame()
    rownames(myscores_mat) <- myscores_mat$patient
    myscores_mat$patient <- NULL
    
    pred <- ROCR::prediction(myscores_mat, myresponse_mat, label.ordering = c("NR", "R"))
    AUC_iHet_exclusion <- performance(pred, measure = "auc")
    AUC_iHet_exclusion <- mean(unlist(AUC_iHet_exclusion@y.values))
    cat("Average (bootstrap) AUC_iHet_exclusion: ", AUC_iHet_exclusion, "\n")
    
    ### iHet
    # predictions matrix
    myscores_mat <- myscores_df
    if (with_TMB) myscores_mat$TMB <- myscores_mat$TMB_cat <- myscores_mat$TMB_lin <- NULL
    myscores_mat$iHet_exclusion <- myscores_mat$iHet_removed <- myscores_mat$iHet_inverted <- NULL
    myscores_mat$response <- myscores_mat$dataset <- NULL
    myscores_mat <- myscores_mat %>% pivot_wider(names_from = run, values_from = iHet) %>% as.data.frame()
    rownames(myscores_mat) <- myscores_mat$patient
    myscores_mat$patient <- NULL
    
    pred <- prediction(myscores_mat, myresponse_mat, label.ordering = c("NR", "R"))
    AUC_iHet <- performance(pred, measure = "auc")
    AUC_iHet <- mean(unlist(AUC_iHet@y.values))
    cat("Average (bootstrap) AUC_iHet: ", AUC_iHet, "\n")
    
  }else{
    ## median
    ### iHet_Exclusion
    pred <- prediction(myscores_df$iHet_exclusion, myscores_df$response, label.ordering = c("NR", "R"))
    AUC_iHet_exclusion <- performance(pred, measure = "auc")
    AUC_iHet_exclusion <- AUC_iHet_exclusion@y.values[[1]]
    cat("AUC_iHet_exclusion: ", AUC_iHet_exclusion, "\n")

    ### iHet
    pred <- prediction(myscores_df$iHet, myscores_df$response, label.ordering = c("NR", "R"))
    AUC_iHet <- performance(pred, measure = "auc")
    AUC_iHet <- AUC_iHet@y.values[[1]]
    cat("AUC_iHet: ", AUC_iHet, "\n")
  }
  
  # Integration of - iHet_exclusion, iHet_removed and TMB

    if ("iHet" %in% names(myscores_df)) {
    myscores_df$iHet_removed <- myscores_df$iHet
    AUC_iHet_removed <- AUC_iHet
  }
  
  if (with_TMB) {
    
    if(length(myscores_df$patient) > length(unique(myscores_df$patient))){
      
      ## bootstrap
      data.tri <- do.call(rbind, lapply(1:100, function(bootstrap_sample){
        myscores_df_run <- subset(myscores_df, run == bootstrap_sample)
        data.tri <- do.call(rbind, lapply(grid, function(i.removed){
          #cat("iHet_removed: ", i.removed, "\n")
          do.call(rbind, lapply(grid, function(i.excluded){
            #cat("iHet_excluded: ", i.excluded, "\n")
            do.call(rbind, lapply(grid, function(i.TMB){
              #cat("TMB: ", i.TMB, "\n")
              if (abs(1-(i.removed+i.excluded+i.TMB)) < 0.0001){
                
                ## duo
                pred_averaged_tmp = i.removed*myscores_df_run$iHet_removed + 
                  i.excluded*myscores_df_run$iHet_exclusion
                
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
                return(data.frame(Top=i.TMB, Right=i.removed, Left=i.excluded, value=AUC_averaged,  run = bootstrap_sample))
              }else{
                return(NULL)
              }
            }))
          }))
        }))
      }))
      
      data.tri <- data.tri %>% group_by(Top, Right, Left) %>% 
        dplyr::summarise(mean_value = mean(value),
                         stdev_value = sqrt(var(value)),
                         n_value = length(value))
      
      data.tri$ciw_value <- qt(0.975, data.tri$n_value) * data.tri$stdev_value / sqrt(data.tri$n_value) # 95% CI
      data.tri$is.max <- FALSE
      data.tri$is.max[data.tri$mean_value == (max(data.tri[,"mean_value"]))] <- TRUE
      cat("Average (bootstrap) AUC_max: ", max(data.tri[,"mean_value"]), "\n")
      
      # Add variable for the use of equal weights
      data.tri$equal.weight <- FALSE
      data.tri$equal.weight[data.tri$Right == 1/3 & data.tri$Right == 1/3 & data.tri$Left == 1/3] <- TRUE
      cat("Average (bootstrap) AUC_equal.weight: ", data.tri$mean_value[data.tri$equal.weight], "\n")
      
    }else{
      
      ## median
      data.tri <- do.call(rbind, lapply(grid, function(i.removed){
        #cat("iHet_removed: ", i.removed, "\n")
        do.call(rbind, lapply(grid, function(i.excluded){
          #cat("iHet_excluded: ", i.excluded, "\n")
          do.call(rbind, lapply(grid, function(i.TMB){
            #cat("TMB: ", i.TMB, "\n")
            if (abs(1-(i.removed+i.excluded+i.TMB)) < 0.0001){
              
              ## duo
              pred_averaged_tmp = i.removed*myscores_df$iHet_removed + 
                i.excluded*myscores_df$iHet_exclusion
              
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
              pred_averaged = pred_averaged_lin + i.TMB*myscores_df$TMB_lin
              
              pred <- prediction(pred_averaged, myscores_df$response, label.ordering = c("NR", "R"))
              AUC_averaged <- performance(pred, measure = "auc")
              AUC_averaged <- AUC_averaged@y.values[[1]]
              return(data.frame(Top=i.TMB, Right=i.removed, Left=i.excluded, value=AUC_averaged))
            }else{
              return(NULL)
            }
          }))
        }))
      }))

      data.tri$is.max <- FALSE
      data.tri$is.max[data.tri$value == (max(data.tri[,"value"]))] <- TRUE
      cat("AUC_max: ", max(data.tri[,"value"]), "\n")
      
      # Add variable for the use of equal weights
      data.tri$equal.weight <- FALSE
      data.tri$equal.weight[data.tri$Right == 1/3 & data.tri$Right == 1/3 & data.tri$Left == 1/3] <- TRUE
      cat("AUC_equal.weight: ", data.tri$value[data.tri$equal.weight], "\n")
      
    }
      
    values <- subset(data.tri, is.max == TRUE)[,1:3]
    palette <- c( "#FF9933", "#002C54", "#3375B2", "#CCDDEC", "#BFBFBF", "#000000")
    max_value <- subset(data.tri, is.max == TRUE)[,4:ncol(data.tri)]
    equal_value <- subset(data.tri, equal.weight == TRUE)[,4:ncol(data.tri)]
    
    if(length(myscores_df$patient) > length(unique(myscores_df$patient))){
      ## bootstrap
      
      g <- ggtern(data=data.tri, ggtern::aes(x = Left, y = Top, z = Right)) +
        theme_bw() + theme_hidetitles() + theme_hidearrows() +
        geom_point(shape=21,size=5, aes(fill=mean_value, col=is.max, alpha=is.max)) +
        geom_crosshair_tern(data = subset(data.tri, is.max==TRUE), colour = "black") + 
        geom_crosshair_tern(data = subset(data.tri, equal.weight==TRUE), lty=2, colour = "black") + 
        # scale_fill_gradient(name="AUC", low="white", high=palette[2], limit=c(0.5, 0.9))+
        # scale_fill_gradient(name="AUC", trans = "log", low="#264653", high="#e76f51", limit=c(0.5, 0.9))+
        scale_alpha_manual(values = c(0.6,0.9))+
        scale_fill_viridis(name="AUC", option = "A", limit=c(0.5, 0.9))+
        scale_colour_manual(values = c("white", "black")) + 
        ggtitle(paste0(dataset_name, 
                       "\niHet: ", round(AUC_iHet_removed,2), 
                       "\niHet_IExcl: ", round(AUC_iHet_exclusion,2), 
                       "\nTMB: ", round(AUC_TMB,2), 
                       "\nmax: ", paste0(round(max_value[,"mean_value"],2),
                                         " (95% CI: ", 
                                         round(max_value[,"mean_value"]-max_value[,"ciw_value"],3), ", ", 
                                         round(max_value[,"mean_value"]+max_value[,"ciw_value"],3), ")"),
                       "\nequal: ", paste0(round(equal_value[,"mean_value"],2),
                                         " (95% CI: ", 
                                         round(equal_value[,"mean_value"]-equal_value[,"ciw_value"],3), ", ", 
                                         round(equal_value[,"mean_value"]+equal_value[,"ciw_value"],3), ")"), "\n",
                       "\nweights:","\nTMB=", as.character(round(values[1],2)), 
                       ",\niHet=", as.character(round(values[2],2)),
                       ",\nIExcl=", as.character(round(values[3],2))))
      g = g +
        Tlab("TMB") +
        Rlab("iHet") +
        Llab("iHet_IExcl") +
        theme_bw() +
        theme(plot.title=element_text(vjust = -30), axis.text = element_text(colour = "black"),
              axis.title = element_text(face = "bold"))
      
      print(g)
      if (save_figure_format == ".pdf"){
        ggsave(filename = paste0("./iHet_ms/iHet_iHetExcl_TMB_", normalize,  "_" ,file_name, "_", tissue, "_performance_results_bootstrap.pdf"),
               width = 8, height = 6)
      }else if (save_figure_format == ".png"){
        ggsave(filename = paste0("./iHet_ms/iHet_iHetExcl_TMB_", normalize,  "_" ,file_name, "_", tissue, "_performance_results_bootstrap.png"),
               width = 8, height = 6)
      }
      
  
    }else{
      ## median
      
      g <- ggtern(data=data.tri,ggtern::aes(x = Left, y = Top, z = Right)) +
        theme_bw() + theme_hidetitles() + theme_hidearrows() +
        geom_point(shape=21,size=5, aes(fill=value, col=is.max, alpha=is.max)) +
        geom_crosshair_tern(data = subset(data.tri, is.max==TRUE), colour = "black") + 
        # scale_fill_gradient(name="AUC", low="white", high=palette[2], limit=c(0.5, 0.9))+
        # scale_fill_gradient(name="AUC", trans = "log", low="#264653", high="#e76f51", limit=c(0.5, 0.9))+
        scale_alpha_manual(values = c(0.6,0.9))+
        scale_fill_viridis(name="AUC", option = "A", limit=c(0.5, 0.9))+
        scale_colour_manual(values = c("white", "black")) + 
        ggtitle(paste0(dataset_name, 
                       "\niHet (IR): ", round(AUC_iHet_removed,2), 
                       "\niHet (noExcl): ", round(AUC_iHet_exclusion,2), 
                       "\nTMB: ", round(AUC_TMB,2), 
                       "\nmax: ", round(max(data.tri[,"value"]),2), "\n",
                       "\nweights:","\nTMB=", as.character(round(values[1],2)), 
                       ",\niHet(IR)=", as.character(round(values[2],2)),
                       ",\niHet(noExcl)=", as.character(round(values[3],2))))
      g = g +
        Tlab("TMB") +
        Rlab("iHet\n(IR)") +
        Llab("iHet\n(noExcl)") +
        theme_bw() +
        theme(plot.title=element_text(vjust = -30), axis.text = element_text(colour = "black"),
              axis.title = element_text(face = "bold"))
      
      print(g)
      
      if (save_figure_format == ".pdf"){
        ggsave(filename = paste0("./iHet_ms/iHet_iHetExcl_TMB_", normalize,  "_" ,file_name, "_", tissue, "_performance_results_median.pdf"),
               width = 8, height = 6)
      }else if (save_figure_format == ".png"){
        ggsave(filename = paste0("./iHet_ms/iHet_iHetExcl_TMB_", normalize,  "_" ,file_name, "_", tissue, "_performance_results_median.png"),
               width = 8, height = 6)
      }
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
