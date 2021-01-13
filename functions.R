# Load libraries
library(MOFA2)
library(reticulate)
library(dplyr)
library(DESeq2)
library(stringr)
library(operators)
library(gtools)
library(matrixStats)
library(doBy)
library(reshape2)
library(sjmisc)
library(ComplexHeatmap)
library(Cairo)
library(RColorBrewer)
library(UpSetR)
library(corrplot)
library(ggpubr)
library(data.table)
library(lares)
library(extrafont)
library(ggrepel)
library(ggplot2)
library(EnhancedVolcano)
library(rstatix)
library(cowplot)
library(beeswarm)
library(ggbeeswarm)
library(ggExtra)
library(tikzDevice)
library(heatmaply)
library(plotly)


getData <- function(type){
  # Import data
  
  datafilepath = paste0('./Data/',type,"_data.RData")
  load(datafilepath)
  
  colnames(cellfrac) <- gsub(".","_",colnames(cellfrac),fixed=TRUE)
  colnames(cellfrac)[9] <- "T_regulatory"
  
  # Load proxy scores separate
  
  if (type!="PBMC"){
    load(paste0('./Data/Proxys_scores_',type,"_data.RData"))
    proxy <- Proxy_measures.matrix
    rownames(proxy) <- rownames(pathways)}
  
  # Delete normal samples from Sharma RNA dataset and rename
  if (type=="SHARMA"){
    
    counts<-counts[,-grep("_N_RNA",colnames(counts))]
    cellfrac<-cellfrac[-grep("_N_RNA",rownames(cellfrac)),]
    lrpairs<-lrpairs[-grep("_N_RNA",rownames(lrpairs)),]
    tfs <- tf_act[-grep("_N_RNA",rownames(tf_act)),]
    pathways<-pathways[-grep("_N_RNA",rownames(pathways)),]
    proxy<-proxy[-grep("_N_RNA",rownames(proxy)),]
    
    colnames(counts) <- substr(colnames(counts), 1, 5)}
  
  # Rename Jia RNA dataset
  if (type=="JIA"){
    tfs <- tf_act
    colnames(counts) <- substr(colnames(counts), 1, 7)}
  
  
  
  # Collect data into single list
  
  if (type!="PBMC"){
    data =  mget(c("counts","cellfrac","pathways","tfs","proxy"))
  } else{
    data =  mget(c("counts","cellfrac","pathways","tfs"))}
  
  print(paste0(type," imported")) 
  return(data)}

processRNAseq <- function(counts,regression=1,group=NULL) {
  
  ### Normalize RNAseq data ###
  RNAseqNorm<-function(counts, mincounts=3, minsamples=0.10) {
    
    # Keep only genes with at least 'mincounts' counts
    # in at least 'minsamples' % of samples
    
    ind.RNAseq<-which(rowSums(counts>=mincounts)/ncol(counts)>=minsamples)
    counts<-counts[ind.RNAseq, ]
    
    # To integers
    counts<-as.matrix(counts)
    mode(counts) <- "integer"
    
    # Apply VST normalization
    dset<-DESeqDataSetFromMatrix(counts, 
                                 colData=data.frame(id=colnames(counts)), 
                                 design=~1)
    dset<-estimateSizeFactors(dset)
    dset<-estimateDispersions(dset)
    normdata<-getVarianceStabilizedData(dset)
    
    return(normdata)
  }
  RNAseq.norm <- RNAseqNorm(counts)
  
  ### Select genes based on variability rank ###
  PatientData <- RNAseq.norm
  
  # Calculate between-patient variability (TCGA) #
  
  if (str_contains(colnames(PatientData),"TCGA") || str_contains(colnames(PatientData),"pbmc")  ){
    print("Processing TCGA RNA-seq")
    
    vars <- rank(rowVars(PatientData))
    
    selection <- which.maxn(vars, 5000)
    
    # Select the RNAseq.norm data based on the ranks
    RNAseq.norm <- RNAseq.norm[selection,]
    
    # Order the samples
    RNAseq.norm<- RNAseq.norm[, mixedsort(colnames(RNAseq.norm))]
    
    # Convert to dataframe for MOFA
    RNAseq.norm.df <- as.data.frame(as.table(RNAseq.norm))
    
    # Make column names
    colnames(RNAseq.norm.df) <- c("feature", "sample", "value")
    
    # Add view as column
    RNAseq.norm.df$view <- "RNAseq"
    
    RNAseq.norm.df$group <- group
    
    # Calculate within-patient variability (Jia & Sharma) #
    
  } else {
    print("Processing JIA/SHARMA RNA-seq")
    
    
    colnames(PatientData) <- sapply(colnames(PatientData), function(x){if(str_length(x) ==5){substr(x,1,2)} else{substr(x,1,4)}})
    Patients <- unique(colnames(PatientData))
    
    # Create an empty matrix with Patients as columns and rows as features (genes)
    var.matrix <- matrix(NA, ncol = length(Patients), nrow= nrow(PatientData))
    colnames(var.matrix) <- Patients
    rownames(var.matrix) <- rownames(PatientData)
    
    # Calculate the variability within patients, consequently rank the genes/features
    for(i in 1:length(Patients)){
      mask <- colnames(PatientData) %in% Patients[i]
      subset <- PatientData[, mask]
      var.matrix[,i] <- rank(rowVars(subset))
    }
    
    # Calculate the average ranks across the patients. 
    average.ranks <-rowMeans(var.matrix)
    
    # Only select the 5000 genes/features with the highest rank
    selection <- which.maxn(average.ranks, 5000)
    features <- names(average.ranks)[selection]
    
    # Select the RNAseq.norm data based on the ranks
    RNAseq.norm <- RNAseq.norm[features,]
    
    
    # Order the samples
    RNAseq.norm<- RNAseq.norm[, mixedsort(colnames(RNAseq.norm))]
    
    # Convert to dataframe for MOFA
    RNAseq.norm.df <- as.data.frame(as.table(RNAseq.norm))
    
    # Make column names
    colnames(RNAseq.norm.df) <- c("feature", "sample", "value")
    
    # Add view as column
    RNAseq.norm.df$view <- "RNAseq"
    
    # Add group
    RNAseq.norm.df$group <- sapply((RNAseq.norm.df$sample), function(x){if(str_length(x) ==5){"Sharma"} else{"Jia"}})
    
    if(regression == 1){
      print("Processing regression")
      dataLM <- RNAseq.norm.df
      
      # Get the patients from the ID number
      dataLM$patient <- as.factor(sapply(dataLM$sample, function(x){if(str_length(x) ==5){substr(x, 1, 2)} else{substr(x, 1, 4)}}))
      
      # Fit a linear model
      lm.out <- lm(value~patient, data=dataLM)
      
      # Extract the residuals from the linear model and assign to the RNAseq data
      residuals <- lm.out[["residuals"]]
      RNAseq.norm.df$value <- residuals
      
    }
  }
  return(RNAseq.norm.df) }

processFeatures<-function(data,group=NULL){
  
  # Transpose data
  cellfrac <- data.matrix(t(data$cellfrac))
  #lrpairs <- data.matrix(t(data$lrpairs))
  pathways <- data.matrix(t(data$pathways))
  tfs <- data.matrix(t(data$tfs))
  
  # Processing features
  
  # Remove missing entries
  #lrpairs <- lrpairs[complete.cases(lrpairs),]
  
  # Remove cell type 'Others'
  cellfrac <- cellfrac[-which(rownames(cellfrac)=="Other"),]
  # Log transformation
  cellfrac <- log10(cellfrac*100+0.001)
  
  # Scale pathways
  pathways <- t(scale(t(pathways)))
  
  
  # Convert to dataframe
  cellfrac.df <- as.data.frame(as.table(cellfrac))
  #lrpairs.df <- as.data.frame(as.table(lrpairs))
  pathways.df <- as.data.frame(as.table(pathways))
  tfs.df <- as.data.frame(as.table(tfs))
  
  # Add 'view' as column
  cellfrac.df["view"] <- "Immune cells quantification"
  #lrpairs.df['view'] <- "Ligand-receptor pairs"
  pathways.df["view"] <- "Pathway scores"
  tfs.df["view"] <- "Transcription factors"
  
  # Rename columns
  colnames(cellfrac.df) <- c("feature", "sample", "value", "view" )
  #colnames(lrpairs.df) <- c("feature", "sample", "value", "view" )
  colnames(tfs.df) <- c("feature", "sample", "value", "view" )
  colnames(pathways.df) <- c("feature", "sample", "value", "view" )
  
  
  # Check if multiple samples exist per patient
  
  Patients <- sapply(colnames(data$counts), function(x){if(str_length(x) ==5){substr(x,1,2)} else{substr(x,1,4)}})
  
  if (!str_contains(Patients,"TCGA") && !str_contains(Patients,"pbmc")){
    print("Processing regression")
    
    linRegression <- function(dataframe.df){
      dataLM <- dataframe.df
      
      # Add the patients as a new column
      dataLM$patient <- as.factor(substr(dataLM$sample, 1, 4))
      
      # Fit a linear model
      lm.out <- lm(value~patient, data=dataLM)
      residuals <- lm.out[["residuals"]]
      
      # Assign the residuals as new values
      dataframe.df$value <- residuals
      
      return(dataframe.df)
    }
    
    cellfrac.df <- linRegression(cellfrac.df)
    #lrpairs.df <- linRegression(lrpairs.df)
    tfs.df <- linRegression(tfs.df)
    pathways.df <- linRegression(pathways.df)
  }
  
  MOFAdata <- do.call("rbind", list(cellfrac.df, pathways.df, tfs.df))
  MOFAdata$group <- group
  
  return(MOFAdata)
}

generateHeatmap <- function(data,title,row_names = F){
  
  data <- acast(data, feature ~sample, value.var = "value")
  
  if (sum(grepl("TCGA",colnames(data)))>1){
    patients <- NULL
    top_ann <- NULL
  }else{
    patients <- substr(colnames(data), 1,4)
    colors = colorRampPalette(rev(brewer.pal(n =7, name = "Paired")))(length(unique(patients)))
    
    ann_colors = list(patient = setNames(colors, unique(patients)))
    
    top_ann <- HeatmapAnnotation(patient = patients,  col = ann_colors)}
  
  
  if (title=="Transcription factors"){
    row_names = F}
  
  
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  col_fun <- circlize::colorRamp2(seq(-4,4, by =8/99), color)
  
  
  
  data <-t(scale(t(data)))
  
  
  
  data[is.nan(data)] <- 0
  
  set.seed(30)
  
  
  min.lim =round(min(data, na.rm = T)-1)
  max.lim = round(max(data, na.rm = T)+1)
  col_fun <- circlize::colorRamp2(seq(min.lim, max.lim, by = (max.lim-min.lim)/99), color)
  
  
  hm <- Heatmap(data, column_title = title, col = col_fun,top_annotation = top_ann, show_row_dend = F, column_title_gp = gpar(fontsize = 18), row_names_gp = gpar(fontsize = 17),column_names_gp = gpar(fontsize = 17),
                show_column_dend = T,  name = "z-score", cluster_columns  = F, cluster_rows = T, show_column_names =  F, show_row_names = row_names)
  
  draw(hm)
  
}

plotResults <- function(model){
  # Plot results
  
  # Explore explained variance per view and per factor
  ## Values
  print(model@cache$variance_explained$r2_total)
  print(model@cache$variance_explained$r2_per_factor)
  
  ## Plotting
  # Total explained variance per view
  
  print(plot_variance_explained(model, x= "view", y ="factor", plot_total = TRUE)[[1]])
  
  print(plot_variance_explained(model, x= "view", y ="factor", plot_total = TRUE)[[2]])
  
  if (model@dimensions$G>1){
    
    ## Plot group results
    
    # Total explained variance per view
    plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
    
    # Variance per feature per view
    print(plot_variance_explained(model, x="group", y="factor"))}
}

buildHeatmap <- function(model){
  
  
  
  col_fun = circlize::colorRamp2(seq(-1, 1, by = 2/10), rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                              "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                              "#4393C3", "#2166AC", "#053061")))
  
  
  
  setup <- function(weights){
    
    
    
    # Define vectors to split the heatmap
    row.split.vector <- c(rep("Immune cell\n Quantification",10), rep("Pathway\n Scores",14), rep("Transcription Factor Activity", nrow(weights)-24))
    
    
    if (paste(names(model@dimensions$N),collapse="") == "LUADLUSC"){
      title = "NSCLC"}
    else if ( paste(names(model@dimensions$N),collapse="") == "JiaSharma"){
      title = "J+S"}
    else{
      title = names(model@dimensions$N)}
    
    
    
    # Size of the each cell
    setsize = 0.30 #cm
    hm <- Heatmap(weights, 
                  name = "weights", 
                  rect_gp = gpar(col = "grey", lwd = 0.5), 
                  col = col_fun, border= T,  
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  show_row_names = T, show_column_names = T, row_names_side = "left",
                  
                  # Annotation
                  # bottom_annotation  = HeatmapAnnotation(
                  #`correlated factors` = combiFactorsFeatures,
                  # col = col.list.Features, 
                  # simple_anno_size = unit(1, "mm"),
                  # annotation_height = unit(3, "mm"), 
                  # annotation_name_gp = gpar(fontsize = 8, fontface = "bold"), 
                  # annotation_legend_param = list(title_gp = gpar(fontsize = 10, fontface = 2))), 
                  
                  #labels_gp = gpar(fontsize = 7)),
                  #grid_width = unit(5, "mm"),
                  #legend_direction = "horizontal"),
                  
                  
                  row_names_gp = gpar(fontsize = 10),
                  row_split = row.split.vector,
                  row_title_gp = gpar(fontsize =15, fontface = 2),
                  
                  column_names_gp  = gpar(fontsize = 9),
                  column_title_gp  = gpar(fontsize =9, fontface = 2),
                  column_title = title,
                  
                  width = unit(setsize, "cm")*ncol(weights) , height  = unit(setsize, "cm")*nrow(weights))
    
    return(hm) }
  
  
  data <- cbind(do.call( "rbind",get_weights(model,  scale = T)))[,1:3]
  
  sign_vector <- getRotation(corr.feat)
  
  sign_vector <- sign_vector[grepl(paste0("^",paste0(names(model@dimensions$N),collapse=""),"$"),gsub('.{2}$', '', rownames(sign_vector))),]
  
  if (names(model@dimensions$N) == "Sharma"){
    sign_vector <- c(1,1,1)}
  
  hm1 <- setup(t(t(data)*sign_vector))
  
  return(hm1)
}

runGSEA <- function(model, feature_set, factors, p = 0.1){
  
  enrichment.pos <- run_enrichment(model, view = "RNAseq", factors = factors, feature.set = feature_set, sign = "positive", statistical.test = "parametric" )
  enrichment.neg <- run_enrichment(model, view = "RNAseq", factors = factors, feature.set = feature_set, sign = "negative", statistical.test = "parametric")
  
  enrichment.pos.pval.adj <- data.frame(enrichment.pos$pval.adj)
  enrichment.neg.pval.adj <- data.frame(enrichment.neg$pval.adj)
  
  # Merge the results of the enrichment with positive and negative weights.
  enrichment.merged <- cbind(enrichment.pos.pval.adj, enrichment.neg.pval.adj)
  rownames(enrichment.merged) <- rownames(enrichment.pos.pval.adj)
  
  
  enrichment.merged[is.na(enrichment.merged)] <- 1
  
  enrichment.merged <- apply(enrichment.merged, c(1,2), function(x){if(x < p){x =x} else{x = NA}})
  
  # Convert p-value to 'score'
  enrichment.merged <- -log10(enrichment.merged)
  
  # Remove rows with only NAs
  ind <- apply(enrichment.merged, 1, function(x) all(is.na(x)))
  enrichment.merged <- enrichment.merged[ !ind, ]
  
  return(enrichment.merged)
}

getGSEAdata <- function(model,gene_annotation){
  out <- runGSEA(model, get(gene_annotation), 1:3, p = 0.1)
  out <- as.data.frame(out)
  return(out)}


plotGSEA <- function(list_of_GSEA_results,type){
  
  getTopGenes <- function(single_model){
    return(unlist(lapply(seq(6), function(x) {rownames(single_model[x])[which.maxn(single_model[,x],10)]} )))}
  
  
  if (type == "GO"){
    genes <- table(unlist(lapply(list_of_GSEA_results,function(x) {getTopGenes(x)})))
    genes <- genes[genes>1]
    genes <- as.matrix(genes[order(-genes), drop=FALSE])
    genes <- rownames(genes)
    
  }else{
    genes <- unique(unlist(lapply(list_of_GSEA_results,function(x) {rownames(x)})))}
  
  
  getGSEAHeatmap <- function(model,title){
    
    data <- matrix(nrow=length(genes),ncol=1)
    rownames(data) <- genes
    
    
    data <- merge(data,model,by=0,all=TRUE)
    
    data[is.na(data)] <- 0
    
    
    rownames(data) <- data[,1]
    
    data[,1:2] <- NULL
    
    data <- data[genes,]
    
    
    colnames(data) <- rep(c("Factor1","Factor2","Factor3"),2)
    
    
    rownames(data) <- sapply(rownames(data), function(x){str_remove(x, paste0(type,"_"))})
    rownames(data) <- sapply(rownames(data), function(x){str_replace_all(x, "_", " ")})
    
    
    
    sign_vector <- getRotation(corr.RNA)
    
    
    
    sign_vector <- sign_vector[grepl(title,rownames(sign_vector)),]
    
    
    
    sign.vector <- c(sign_vector,sign_vector*-1)
    
    
    if (title == "Jia-Sharma"){
      sign.vector <- c(rep(-1,3),rep(1,3))}
    
    if (title == "Jia"){
      sign.vector <- c(-1,-1,1,1,1,-1)}
    
    if (title == "Sharma"){
      sign.vector <- c(rep(1,3),rep(-1,3))}
    
    
    #sign.vector <- c(rep("positive", ncol(data)/2),rep("negative", ncol(data)/2)) 
    
    
    sign.color<- setNames((c("#B2182B", "#2166AC")) ,(c(1, -1)))
    col.list <- list(sign = sign.color)
    
    
    max.lim <-  round(max(data, na.rm = T))+1
    # Define color function for the heatmaps. Colors chosen based on corrplot
    col_fun = circlize::colorRamp2(seq(0,max.lim, by =(max.lim)/5), c("#FFFFFF", "#D1E5F0", "#92C5DE",
                                                                      "#4393C3", "#2166AC", "#053061"))
    setsize=0.5
    
    hm <- Heatmap(as.matrix(data), 
                  name = "-log10(p-value)", 
                  rect_gp = gpar(col = "grey", lwd = 0.5), 
                  col = col_fun ,border= TRUE,  na_col = "white",
                  cluster_rows = F, cluster_columns = FALSE,
                  show_row_names = TRUE, show_column_names = TRUE, 
                  row_names_side = "left", 
                  heatmap_legend_param = list(
                    title_gp = gpar(fontsize = 15, fontface = 2),
                    title_position = "leftcenter-rot",
                    labels_gp = gpar(fontsize =15),
                    grid_width = unit(2, "mm"),
                    legend_direction = "vertical"
                  ),
                  # Annotation
                  bottom_annotation  = HeatmapAnnotation(
                    sign = sign.vector,
                    col = col.list,
                    simple_anno_size = unit(1, "mm"),
                    annotation_height = unit(3, "mm"),
                    annotation_name_gp = gpar(fontsize = 15, fontface = "bold"), show_annotation_name = F,
                    annotation_legend_param = list(title_gp = gpar(fontsize = 15, fontface = 2))),
                  
                  # Rows
                  row_names_gp = gpar(fontsize =12), 
                  row_title_gp = gpar(fontsize =12, fontface = 2),
                  row_gap = unit(2, "mm"),
                  
                  # Columns
                  column_gap = unit(1, "mm"),
                  column_names_gp  = gpar(fontsize = 15),
                  column_title_gp  = gpar(fontsize =15, fontface = 2),
                  column_title = title,
                  cluster_row_slices = FALSE, 
                  cluster_column_slices = FALSE,
                  # Size heatmap
                  width = unit(setsize, "cm")*ncol(data), height  = unit(setsize, "cm")*nrow(data), show_heatmap_legend = T)
    
    return(hm)}
  
  
  
  draw(  getGSEAHeatmap(list_of_GSEA_results[["RNA.jia"]],"Jia") +  getGSEAHeatmap(list_of_GSEA_results[["RNA.sharma"]],"Sharma") + 
           getGSEAHeatmap(list_of_GSEA_results[["RNA.comb.jiasharma"]],"Jia-Sharma") + getGSEAHeatmap(list_of_GSEA_results[["RNA.lusc"]],"LUSC") + 
           getGSEAHeatmap(list_of_GSEA_results[["RNA.luad"]],"LUAD") + getGSEAHeatmap(list_of_GSEA_results[["RNA.nsclc"]],"NSCLC") +
           getGSEAHeatmap(list_of_GSEA_results[["RNA.crc"]],"CRC") +    getGSEAHeatmap(list_of_GSEA_results[["RNA.skcm"]],"SKCM") + getGSEAHeatmap(list_of_GSEA_results[["RNA.pbmc"]],"PBMC")) }




runMOFA <- function(MOFAdata){
  
  MOFAobject <- create_mofa(get(MOFAdata))
  
  
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  train_opts <- get_default_training_options(MOFAobject)
  model_opts$num_factors <-7# number of samples is very small for learning 10 factors, it should not exceed ~7
  train_opts$convergence_mode <- "medium"
  train_opts$maxiter <- 7000
  #train_opts$drop_factor_threshold <- 0
  use_condaenv("C:/Users/s149072/Anaconda3", required = TRUE)
  
  set.seed(1234)
  
  # Run MOFA
  MOFAobject <- prepare_mofa(object = MOFAobject, data_options = data_opts , model_options = model_opts, training_options = train_opts)
  
  model <- run_mofa(MOFAobject, outfile = paste0("./Models/",MOFAdata,".hdf5"))
  
  set.seed(NULL)
  
  return(model) }


runMOFA_bootstrap <- function(MOFAdata,type,N){
  
  MOFAobject <- create_mofa(MOFAdata)
  
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  train_opts <- get_default_training_options(MOFAobject)
  model_opts$num_factors <-7# number of samples is very small for learning 10 factors, it should not exceed ~7
  train_opts$convergence_mode <- "medium"
  train_opts$maxiter <- 7000
  #train_opts$drop_factor_threshold <- 0
  use_condaenv("C:/Users/s149072/Anaconda3", required = TRUE)
  
  set.seed(1234)
  
  # Run MOFA
  MOFAobject <- prepare_mofa(object = MOFAobject, data_options = data_opts , model_options = model_opts, training_options = train_opts)
  
  model <- run_mofa(MOFAobject, outfile = paste0("./Models/Bootstrap/",type,N,".hdf5"))
  
  set.seed(NULL)
  
  return(model) }

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

transformFactors <- function(model){
  
  
  data <- data.frame(do.call("rbind",get_weights(model)))[1:3] 
  colSums(data)
  
  if (paste(names(model@dimensions$N),collapse="") == "CRCLUADLUSCSKCM"){
    colnames(data) <- paste("TCGA",substr(colnames(data),start=7,stop=7),sep=" ")
  } else if (paste(names(model@dimensions$N),collapse="") == "LUADLUSC"){
    colnames(data) <- paste("NSCLC",substr(colnames(data),start=7,stop=7),sep=" ")
  }else{
    colnames(data) <- paste(paste(names(model@dimensions$N),collapse=""),substr(colnames(data),start=7,stop=7),sep=" ")}
  return(data)}

getCorr <- function(list_of_models,RNA_corr=FALSE){
  factors <- lapply(list_of_models,transformFactors)
  
  cnames <- Reduce(intersect,lapply(factors,rownames))
  
  factors <- lapply(factors,function(i) i[cnames,])
  
  data <- do.call(cbind,factors)
  
  colnames(data) <- unlist(lapply(factors,colnames))
  
  cor_matrix <- cor(data)

  rotation <- getRotation(cor_matrix)
  
  
  data[rownames(rotation)] <-  t(t(data[rownames(rotation)])*rotation$Sign)
  
  if (RNA_corr == FALSE){
    data['JS B1'] <- get_median_bootstrap_data(jiasharma)
    data['NSCLC B1'] <- get_median_bootstrap_data(nsclc)}
  
  cor_matrix_plot <- cor(data)
  
  p.mat = cor.mtest(cor_matrix_plot)
  
  
  cor_matrix_plot[p.mat>0.05] <- 1
  
  p.mat[p.mat>0.05] <- 0.00

  plot <- heatmaply_cor(
    cor_matrix_plot,
    node_type = "scatter",
    column_text_angle = 90,
    fontsize_row = 13,
    fontsize_col = 13,
    point_size_mat = -log10(p.mat),
    point_size_name = "-log10",
    label_names = c("x", "y", "Pearson \ncorrelation") ) 
  
  ggplotly(plot)
  
  cor_matrix_plot[p.mat>0.05] <- 0

  
  corrplot(cor_matrix_plot, tl.col="black", tl.srt=75, tl.cex=1.2, p.mat = p.mat, sig.level = 0.001, insig = "blank")
  
  
  
  return(cor_matrix)
  
  
}

getRotation <- function(corr_matrix){

  
  corr_matrix <- corr_matrix[,22]
  
  
  
  dir_factors <- as.data.frame(sign(corr_matrix))
  
  dir_factors[dir_factors==0] = 1
  
  colnames(dir_factors) <- "Sign"
  
  
  return(dir_factors)}


plotResponseCorr <- function(model,type,corr){
  
  
  # Get first three factors explaining largest variance
  factors <- data.frame(do.call("rbind" ,get_factors(model, scale = T))[,1:3])
  
  # Get rotation and apply to factors
  sign_vector <- getRotation(corr)
  
  sign_vector <- sign_vector[grepl(paste0("^",paste0(names(model@dimensions$N),collapse=""),"$"),gsub('.{2}$', '', rownames(sign_vector))),]
  
  if (names(model@dimensions$N) == "Sharma"){
    sign_vector <- c(1,1,1)}
  factors <- t(t(factors)*sign_vector)
  
  
  getProxyData <- function(name){
    # Get proxy data and calculate the median
    proxy <- get(name)$proxy
    proxy <- getProxyMedian(proxy)
    return(proxy)}
  
  
  if (names(model@dimensions$N) == "NSCLC"){
    proxy_luad <- getProxyData("LUAD")
    proxy_lusc <- getProxyData("LUSC")
    proxy <- c(proxy_luad,proxy_lusc)
    # Combine proxy and factor weights
    data <- as.data.frame(cbind(factors,proxy))
    
  } else if (length(names(model@dimensions$N)) > 1){
    proxy_jia <- getProxyData("JIA")
    proxy_sharma <- getProxyData("SHARMA")
    proxy <- c(proxy_jia,proxy_sharma)
    # Combine proxy and factor weights
    data <- as.data.frame(cbind(factors,proxy))
    
  }else{
    proxy <- getProxyData(toupper(names(model@dimensions$N)))
    data <- as.data.frame(cbind(factors,proxy)) }
  
  title <- paste0(names(model@dimensions$N),collapse="")
  
  # Prepare data for plot
  table <- as.data.frame(melt(setDT(data), id.vars = c("proxy"), variable.name = "Factor"))
  
  # Generate figure
  fig <- (ggscatter(table, y = "proxy", x= "value", add = "reg.line", conf.int = TRUE, facet.by = "Factor",
                    cor.coef = TRUE, cor.method = "pearson", formula = "x ~ y",
                    xlab = paste0("Factor values"), ylab ="Response score", size = 1, height = 1) )+ ggtitle(paste0(title," ",type))
  
  return(fig)}


getProxyMedian <- function(proxy){
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  
  patterns <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
  
  #filter_scores <- c("CYT", "IS_Davoli", "RohIS", "IFny", "ExpandedImmune", "T_cell_inflamed")
  
  #filter_scores <- grepl(paste(patterns, collapse="|"), colnames(proxy))
  
  #proxy <- proxy[,filter_scores]
  
  proxy <- proxy[,patterns]
  proxy <- apply(proxy, 2, range01)
  proxy <- apply(proxy, 1, median)
  
  return(proxy)}

responseHeatmap <- function(type){
  
  dataset <- get(paste0(type))
  data <- dataset$proxy
  data <- getProxyMedian(data)
  data <- as.data.frame(data)
  
  rownames(data) <- rownames(dataset$pathways)
  
  if (type == "JIA"){
    patients <- unique(substr(rownames(data),1,4))
    title <- paste0("Within-patient response Jia")
    k=7
  }else{
    patients <- unique(substr(rownames(data),1,2))
    title <- paste0("Within-patient response Sharma")
    k=5 }
  
  
  df <- data.frame(matrix(ncol=6,nrow=length(patients)),row.names = patients)
  colnames(df) <- c(seq(1:6))
  
  for (i in seq_along(patients)){
    for (j in seq(1:6)){
      samples <- substr(rownames(data)[grepl(patients[i],rownames(data))],k,k)
      if (length(which(samples == j)>0)) {
        df[i,j] <- data[grepl(patients[i],rownames(data)),][which(samples == j)]}
    }
  }
  
  weights <- as.matrix(df)
  setsize= 0.6
  col_fun = circlize::colorRamp2(seq(0, 1, by = 1/10), c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                         "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                         "#4393C3", "#2166AC", "#053061"))
  row.split.vector <- c(rep("Patient",nrow(weights)))
  colnames(weights) <- lapply(colnames(weights),function(x) paste0("Sample ",x))
  
  hm <- Heatmap(weights, 
                name = "weights", 
                rect_gp = gpar(col = "grey", lwd = 0.5), 
                #col = col_fun,
                border= T,  
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = T, show_column_names = T, row_names_side = "left",
                row_names_gp = gpar(fontsize = 12),
                row_split = row.split.vector,
                row_title_gp = gpar(fontsize =15, fontface = 1),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", weights[i, j]), x, y, gp = gpar(fontsize = 10))},
                column_names_gp  = gpar(fontsize = 12),
                column_title_gp  = gpar(fontsize =15, fontface = 2),
                column_title = title,)
  
  return(hm) }

get_median_bootstrap_data <- function(data,variance_cutoff=0){
  data <- data[data$Variance>variance_cutoff,]
  data$Variance <- NULL
  getMedianperFactor <- function(data,i){
    
    #row <- as.data.frame(colMeans(data[grepl(paste0("Factor",i),rownames(data)),]))
    
    row <- as.data.frame(apply(data[grepl(paste0("Factor",i),rownames(data)),],2,median))
    
    colnames(row) <- paste0("Factor",i)
    return(row)}
  
  data <- as.data.frame(do.call(cbind,lapply(seq(1),function(x) {getMedianperFactor(data,x)})))
  data <- Filter(function(x)!all(is.na(x)), data)
  
  return(data) }

plotHMs_bootstrap <- function(data,title){
  
  
  data <- get_median_bootstrap_data(data)
  
  
  
  setsize= 1.5
  col_fun = circlize::colorRamp2(seq(-1, 1, by = 2/10), rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                              "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                              "#4393C3", "#2166AC", "#053061")))
  
  
  row.split.vector <- c(rep("Immune cells\n Quantification",10), rep("Pathways Scores",14), rep("Transcription Factors Activity", nrow(data)-24))
  
  
  data <- as.matrix(data)
  
  if (title == "NSCLC"){
    corr <- t(cor(get_median_bootstrap_data(jiasharma),get_median_bootstrap_data(nsclc)))
    sign_vector <- as.matrix(sign(corr[cbind(1:nrow(corr), max.col(abs(corr)))]))
    sign_vector <- sign_vector[1:ncol(data)]
    data <- t(t(data)*sign_vector)}
  
  
  hm <- Heatmap(data, 
                name = "weights", 
                rect_gp = gpar(col = "grey", lwd = 0.5), 
                col = col_fun,
                border= T,  
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = T, show_column_names = T, row_names_side = "left",
                
                column_title = paste0(title," bootstrap"),
                column_title_gp  = gpar(fontsize =9, fontface = 2),
                column_title_rot = 90,
                row_names_gp = gpar(fontsize = 10),
                row_split = row.split.vector,
                row_title_gp = gpar(fontsize =15, fontface = 2),
                
                column_names_gp  = gpar(fontsize = 9),
                
                width = unit(setsize*0.2, "cm")*ncol(data) , height  = unit(setsize*0.2, "cm")*nrow(data) )
  
  return(hm)}

Bootstrap_MOFA <- function(type,N){
  
  print(paste0("Iteration N=",N," of dataset: ",type))
  
  model <- processed_features[[type]]
  
  samples = sample(unique(model$sample),replace=TRUE)
  
  getSampleFromData <- function(counter){
    extract <- model[model$sample==samples[counter],]
    extract$sample <- paste0("P",counter)
    return(extract)}
  
  mofa_data <- as.data.frame(do.call(rbind, lapply(seq(length(samples))  , function(x) {getSampleFromData(x)})))
  
  model <- runMOFA_bootstrap(mofa_data,type,N)}

getWeightsBootstrap <- function(type,N){
  
  model <- load_model(paste0("./Models/Bootstrap/",type,N,".hdf5"))
  
  
  factors <- which(rowSums(Reduce('+',model@cache$variance_explained$r2_per_factor))>1)
  
  if (length(factors)>0){
    weights <- do.call( "rbind",get_weights(model,  scale = T))[,factors]
    Variance <- rowSums(Reduce('+',model@cache$variance_explained$r2_per_factor)) [factors]
    weights <- rbind(Variance,weights)
    colnames(weights) <- paste0(substr(type,6,str_length(type)),"_",colnames(weights),"_",N)}
  
  else{
    rows <- rownames(do.call( "rbind",get_weights(model,  scale = T)))
    rows <- c("Variance",rows)
    weights <- data.frame(matrix(,nrow=length(rows),ncol=0))
    rownames(weights) <- rows}
  
  weights <- t(weights)
  return(weights)
}


plotWeightsFactors <- function(model1,model2,xlabel,ylabel){
  
  NSCLC <- model1
  JiaSharma <- model2
  
  hist(NSCLC,100,col=rgb(1,0,0,0.5),main="Histogram of Factor weights",xlab="weights")
  hist(JiaSharma,add=T,100,col=rgb(0,0,1,0.5))
  legend("topright", legend=c(xlabel,ylabel), col=c(rgb(1,0,0,0.5), 
                                                    rgb(0,0,1,0.5)), pt.cex=2, pch=15 )
  
  
  data <- as.data.frame(cbind(NSCLC,JiaSharma))
  colnames(data) <- c("v1","v2")
  
  
  ggscatter(data,x = "v1", y = "v2", add = "reg.line",xlab=xlabel,ylab=ylabel) +
    stat_cor(label.x = 0.0, label.y = 0.-0.5) +
    stat_regline_equation(label.x = 0.0, label.y = -1)
}


plotBootstrapWeights <- function(results1,results2){
  
  js <- results1[grepl("Factor1",rownames(results1)),]
  js$Variance <- NULL
  
  nsclc <- results2[grepl("Factor1",rownames(results2)),]
  nsclc$Variance <- NULL
  
  plotWeightsFactors(colMeans(nsclc),colMeans(js),"NSCLC-boot F1","JS-boot F1")
}


plotHorBarDiff <- function(y,x){
  diff <- (y-x)^2
  diff <- sort(diff,decreasing=TRUE)
  diff <- as.data.frame(diff)
  colnames(diff) <- c("value")
  diff$key <- rownames(diff)
  
  
  ggplot(diff,aes(x=key,y=value)) + geom_bar(stat='identity') + 
    coord_flip() + scale_y_continuous(name="Difference in weight") +
    scale_x_discrete(name="Feature") }


plotRawFeatures <- function(dataset,type){
  
  data <- cbind(dataset$cellfrac,dataset$pathways,dataset$tfs)
  
  data$Other <- NULL
  
  X <- colSds(as.matrix(data))
  
  Y <- abs(results$effsize)
  
  df <- as.data.frame(cbind(X,Y))
  
  rownames(df) <- colnames(data)
  
  df <- df[rownames(df) %in% features_js,]
  
  
  ggplot(data=df,aes(x=X,y=Y)) + geom_point() + scale_y_log10(limits = c(1e-02,1e3)) + scale_x_log10(limits = c(1e-02,1e+01))
  
  
  
  ggplot(data=df,aes(x=X,y=V2)) + geom_point() + scale_y_log10(limits = c(1e-02,1e3)) + scale_x_log10(limits = c(1e-02,1e+03)) + xlab("Between") + ylab("Within") +
    stat_cor(label.x = 0.5, label.y = 0.-0.80) + 
    stat_regline_equation(label.x = 0.5, label.y = -1) +
    geom_smooth(method = "lm", se = FALSE ,colour="black", size=0.8) + annotation_logticks()   +
    
    geom_hline(yintercept=mean(Y),linetype = 'dotted') + geom_vline(xintercept=mean(X),linetype = 'dotted') + theme_bw() + ggtitle(paste0("Within- vs between-patient variability ",type))}



plotOverview <- function(jiasharma,nsclc){
  
  #Extract first factors of bootstrapped data
  js <- jiasharma[grepl("Factor1",rownames(jiasharma)),]
  js$Variance <- NULL
  js <- t(as.data.frame(lapply(js,median)))
  
  nsclc <- nsclc[grepl("Factor1",rownames(nsclc)),]
  nsclc$Variance <- NULL
  nsclc <- t(as.data.frame(lapply(nsclc,median)))
  
  
  data <- as.data.frame(cbind(js,nsclc))
  colnames(data) <- c("JS","NSCLC")
  
  data$Type <- c(rep("Immune cells\n Quantification",10), rep("Pathway Scores",14), rep("Transcription Factor \nActivity", nrow(data)-24))
  
  
  
  getOutliers <- function(dataset){
    a <- sort(rownames(dataset)[dataset<(mean(dataset)-2*sd(dataset))])
    b <- sort(rownames(dataset)[dataset>(mean(dataset)+2*sd(dataset))])
    return(c(a,b))}
  
  
  outliers <- c(getOutliers(js),getOutliers(nsclc))
  
  labels <- rownames(data)
  labels[!(labels %in% outliers)] <- ""
  
  both <- intersect(getOutliers(js),getOutliers(nsclc))
  
  data$Color <- "Non-outlier"
  data$Color[rownames(data) %in% setdiff(getOutliers(js),both)] <- "Unique Jia-Sharma"
  data$Color[rownames(data) %in% setdiff(getOutliers(nsclc),both)] <- "Unique NSCLC"
  data$Color[rownames(data) %in% both] <- "Intersection"
  data$names <- labels
  
  plot <- ggplot(data,aes(JS,NSCLC), add = "reg.line",conf.int=TRUE) + # ggtitle("Median of factor 1 weights Jia-Sharma and NSCLC dataset") + 
    stat_cor(label.x = 0.3, label.y = 0.-0.80, size=6) + 
    stat_regline_equation(label.x = 0.3, label.y = -0.95, size=6) +
    geom_smooth(method = "lm", se = FALSE ,colour="black", size=0.8) +
    geom_point(aes(shape = Type,col=Color), size = 3) +
    scale_shape_manual(values = c(15, 16, 17)) +
    scale_color_manual(values = c("purple", "black", "blue", "red")) + geom_text_repel(aes(label = labels),size=6,face="bold")
  
  
  p <- plot + geom_hline(yintercept=0,linetype = 'dotted') + 
    geom_vline(xintercept=0,linetype = 'dotted') + 
    geom_segment(aes(x = -1, y = -1, xend = 1, yend = 1),linetype='dotted') +  theme(panel.border = element_blank()) + theme_bw() +
    theme(legend.position = c(1.12,0.2),legend.text=element_text(size=13), axis.title = element_text(size = 20), axis.text = element_text(size = 15),plot.title = element_text(size = 20, face = "bold"))
  
  plot(ggMarginal(p, type = "histogram")) 
  
  return(data)
}    

plotVolcano <- function(jiasharma,nsclc,labels){
  
  # Extract first factors of bootstrapped data
  js <- jiasharma[grepl("Factor1",rownames(jiasharma)),]
  js$Variance <- NULL
  
  nsclc <- nsclc[grepl("Factor1",rownames(nsclc)),]
  nsclc$Variance <- NULL
  
  js$group <- "js"
  nsclc$group <- "nsclc"
  
  variables <- js
  variables$group <- NULL
  variables <- colnames(variables)
  
  data <- rbind(js,nsclc)
  
  
  
  getWilcoxon <- function(variable){
    df <- cbind(data$group,data[variable])
    colnames(df) <- c("group","x")
    wilcox <- wilcox_test(df,x~group,p.adjust.method = "BH",detailed = TRUE)

    mean <- mean(abs(df[(df$group=="js"),]$x)) - mean(abs(df[(df$group=="nsclc"),]$x))
    
    matrix <- as.data.frame(cbind(wilcox$p,wilcox_effsize(df,x~group)$effsize,wilcox$estimate,wilcox$statistic,mean))
    rownames(matrix) <- variable
    colnames(matrix) <- c("pvalue","effsize","estimate","statistic","delta_means")
    return(matrix)}
  
  
  
  results <- as.data.frame(do.call( "rbind",lapply(variables,getWilcoxon)))
  
  results$effsize <- results$effsize*sign(results$delta_means)
  
  
  labels <- labels["Color"]
  
  
  keyvals <- c(rep('grey80',nrow(results)))
  
  keyvals[(results$pvalue < 0.01) & (abs(results$effsize)>0.5)]  <- 'grey20'
  
  keyvals[labels$Color == "Unique Jia-Sharma"] <- 'royalblue'
  keyvals[labels$Color == "Unique NSCLC"] <- 'red'
  keyvals[labels$Color == "Intersection"] <- 'purple'
  
  
  
  names(keyvals)[keyvals == 'royalblue'] <- 'Outlier Jia-Sharma'
  names(keyvals)[keyvals == 'red'] <- 'Outlier NSCLC'
  names(keyvals)[keyvals == 'purple'] <- 'Intersection'
  names(keyvals)[keyvals == 'grey80'] <- 'Not significant'
  names(keyvals)[keyvals == 'grey20'] <- 'Significant'
  
  
  
  keyvals.shape <- c(rep(1,nrow(results)))
  
  
  keyvals.shape[labels$Color != "Non-outlier"] <- 19
  
  
  names(keyvals.shape)[keyvals.shape == 1] <- '-'
  names(keyvals.shape)[keyvals.shape == 19] <- '--'
  
  
  
  
  results$Color <- labels$Color
  results$corrsign <- as.factor(sign(results$estimate))
  
  
  results$threshold = as.numeric(as.factor(results$pvalue < 0.05))
  results$threshold<-factor(results$threshold, levels = c(1,2), labels = c("notSign", "Sign"))
  
  
  
  results$names <- rownames(results)
  
  # Calculate dimensions of plot
  xminmax <- max(abs(results$effsize))
  xminmax <- xminmax + xminmax*0.01
  ymax <- max(-log10(results$pvalue))
  ymax <- ymax + ymax*0.01
  
  # Specify important features
  subset_results <- results[results$Color!="Non-outlier",]
  subset_results$Dataset <- subset_results$Color
  subset_results$Specificity <- subset_results$Color
  subset_results$Specificity = "Intra- / inter-variability"
  subset_results$Specificity[abs(subset_results$effsize)<0.5] = "Intra- + inter-variability"
  
  
  g = ggplot(data=results, aes(x=effsize, y=-log10(pvalue)), environment = environment()) +
    geom_point(alpha=0.05) + geom_point(data=subset_results, aes(fill=Dataset, size=Specificity), color="black", alpha=0.8, shape=21) +
    scale_fill_manual(values = c("purple","blue","red")) +
    scale_size_manual(values=c(5,3)) +
    xlim(c(-xminmax, xminmax)) + ylim(c(0, ymax)) +
    xlab("larger in NSCLC                 effect size             larger in Jia-Sharma") + ylab("-log10 p-value") + ggtitle("") +
    theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour="#9e9e9e") + 
    geom_vline(xintercept = 0, linetype = "solid", colour="#9e9e9e") + geom_vline(xintercept = -0.5, linetype = "longdash", colour="#9e9e9e") +
    geom_vline(xintercept = 0.5, linetype = "longdash", colour="#9e9e9e")
  
  g=g+geom_text_repel(data=subset(results, Color!="Non-outlier"), aes(x=effsize, y=-log10(pvalue), label=names), show.legend = NA, inherit.aes = F)
  
  plot(g)
  
}




plotBeeswarm <- function(labels,dataset){
  
  # Import data modalities
  data <- cbind(dataset$cellfrac,dataset$pathways,dataset$tfs)
  
  # Separate patients
  patients <- sapply(rownames(data), function(x){if(str_length(x) == 9){substr(x,1,2)} else{substr(x,1,4)}})
  
  # Calculate predicted immune response score
  data$patient <- patients
  data$proxy <- getProxyMedian(dataset$proxy)
  
  
  # Extract important features
  features <- rownames(labels[labels$Color!="Non-outlier",])
  features <- c(features,'proxy')
  

  
  # Rename JAK-STAT feature
  features[grepl('JAK',features)] <- "JAKSTAT"
  names(data)[names(data)=='JAK-STAT'] <- "JAKSTAT"
  
  
  # Rename proxy feature
  features[grepl('proxy',features)] <- "Response_score"
  names(data)[names(data)=='proxy'] <- "Response_score"
  
  
  getBeeswarm <- function(feature){
    return( ggplot(data, aes_string(x="patient", y=feature,col = "patient")) +  
              geom_beeswarm(size=3,groupOnX=FALSE,priority='density') + theme(text = element_text(size=20),axis.text.x=element_text(angle = 90, vjust = 0.5)) + theme(legend.position = "none" ) )  }
  
  
  # Plot beeswarm
  plot_list <- lapply(features,getBeeswarm)
  
  
  cowplot::plot_grid(plotlist = plot_list,ncol= 3)
  
}
