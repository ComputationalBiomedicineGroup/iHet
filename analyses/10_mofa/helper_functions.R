library(MOFA2)
library(reticulate)
library(dplyr)
library(DESeq2)
library(stringr)
# library(operators)
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
# library(lares)
library(extrafont)
library(ggrepel)
library(ggplot2)
library(EnhancedVolcano)
library(rstatix)
library(cowplot)
library(beeswarm)
library(ggbeeswarm)
library(ggExtra)
# library(tikzDevice)
library(heatmaply)
library(plotly)


#' Bring features into 'tidy' format required by MOFA
processFeatures <- function(data, dataset_id) {
  # Transpose data
  cellfrac <- data.matrix(t(data$cellfrac[[dataset_id]]))
  #lrpairs <- data.matrix(t(data$lrpairs))
  pathways <- data.matrix(t(data$pathway[[dataset_id]]))
  tfs <- data.matrix(t(data$tf[[dataset_id]]))

  # Processing features

  # Remove missing entries
  #lrpairs <- lrpairs[complete.cases(lrpairs),]

  # Remove cell type 'Others'
  cellfrac <- cellfrac[-which(rownames(cellfrac) == "Other"), ]
  # Log transformation
  cellfrac <- log10(cellfrac * 100 + 0.001)

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
  colnames(cellfrac.df) <- c("feature", "sample", "value", "view")
  #colnames(lrpairs.df) <- c("feature", "sample", "value", "view" )
  colnames(tfs.df) <- c("feature", "sample", "value", "view")
  colnames(pathways.df) <- c("feature", "sample", "value", "view")


  # These are datasets with multiple biospies per patient. Regressing
  # out effects related to the inter-patient variability (we want to focus
  # on the intra-patient variability).
  if (dataset_id %in% c("Jia2018", "Sharma2019")) {
    print("Processing regression")

    linRegression <- function(dataframe.df) {
      dataLM <- dataframe.df

      # Add the patients as a new column
      # TODO split by `_` - there's likely an issue with the Sharma data (patient id only 2 chars)
      dataLM$patient <- as.factor(substr(dataLM$sample, 1, 4))

      # Fit a linear model
      lm.out <- lm(value ~ patient, data = dataLM)
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

  MOFAdata <-
    do.call("rbind", list(cellfrac.df, pathways.df, tfs.df))
  MOFAdata$group <- dataset_id

  return(MOFAdata)
}


#' Plot heatmap
generateHeatmap <- function(data, title, row_names = F) {
  data <- acast(data, feature ~ sample, value.var = "value")

  if (sum(grepl("TCGA", colnames(data))) > 1) {
    patients <- NULL
    top_ann <- NULL
  } else {
    patients <- substr(colnames(data), 1, 4)
    colors = colorRampPalette(rev(brewer.pal(n = 7, name = "Paired")))(length(unique(patients)))

    ann_colors = list(patient = setNames(colors, unique(patients)))

    top_ann <-
      HeatmapAnnotation(patient = patients,  col = ann_colors)
  }

  if (title == "Transcription factors") {
    row_names = F
  }

  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  col_fun <- circlize::colorRamp2(seq(-4, 4, by = 8 / 99), color)

  data <- t(scale(t(data)))

  data[is.nan(data)] <- 0

  set.seed(30)

  min.lim = round(min(data, na.rm = T) - 1)
  max.lim = round(max(data, na.rm = T) + 1)
  col_fun <-
    circlize::colorRamp2(seq(min.lim, max.lim, by = (max.lim - min.lim) / 99), color)

  hm <-
    Heatmap(
      data,
      column_title = title,
      col = col_fun,
      top_annotation = top_ann,
      show_row_dend = F,
      column_title_gp = gpar(fontsize = 18),
      row_names_gp = gpar(fontsize = 17),
      column_names_gp = gpar(fontsize = 17),
      show_column_dend = T,
      name = "z-score",
      cluster_columns  = F,
      cluster_rows = T,
      show_column_names =  F,
      show_row_names = row_names
    )

  draw(hm)

}

#' train the MOFA model
runMOFA <- function(MOFAdata){

  MOFAobject <- create_mofa(MOFAdata)

  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  train_opts <- get_default_training_options(MOFAobject)
  model_opts$num_factors <-7# number of samples is very small for learning 10 factors, it should not exceed ~7
  train_opts$convergence_mode <- "medium"
  train_opts$maxiter <- 7000
  #train_opts$drop_factor_threshold <- 0

  set.seed(1234)

  # Run MOFA
  MOFAobject <- prepare_mofa(object = MOFAobject, data_options = data_opts , model_options = model_opts, training_options = train_opts)

  model <- run_mofa(MOFAobject, outfile = paste0("./Models/",MOFAdata,".hdf5"), use_basilisk=FALSE)

  set.seed(NULL)

  return(model)

}