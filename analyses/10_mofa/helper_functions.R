library(conflicted)
library(MOFA2)
library(reticulate)
library(dplyr)
conflict_prefer("filter", "dplyr")
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(tibble)
library(corrplot)
library(stringr)


#' Bring features into 'tidy' format required by MOFA
processFeatures <- function(data, dataset_id) {

  add_patient_col = function(df) {
    if (dataset_id %in% c("Jia2018", "Sharma2019")) {
      df %>% separate("sample", c("patient", "replicate"), remove=FALSE) %>%
        # exclude normal samples
        filter(replicate != "N") %>%
        select(-replicate)
    } else {
      df %>% mutate(patient = sample)
    }
  }

  do_pivot = function(df, view) {
    df %>% as_tibble(rownames="sample") %>%
      pivot_longer(cols=-sample, names_to="feature") %>%
      mutate(view = !!view) %>%
      add_patient_col
  }

  cellfrac.df = do_pivot(data$cellfrac[[dataset_id]], "Immune cells quantification") %>%
    filter(feature != "Other") %>%
    # Log transform
    mutate(value = log10(value * 100 + 0.001))
  pathways.df = do_pivot(scale(data$pathway[[dataset_id]]), "Pathway scores")
  tfs.df = do_pivot(data$tf[[dataset_id]], "Transcription factors")

  # These are datasets with multiple biospies per patient. Regressing
  # out effects related to the inter-patient variability (we want to focus
  # on the intra-patient variability).
  if (dataset_id %in% c("Jia2018", "Sharma2019")) {
    print("Processing regression")

    linRegression <- function(dataframe.df) {
      # Fit a linear model
      lm.out <- lm(value ~ patient, data = dataframe.df)
      residuals <- lm.out[["residuals"]]

      # Assign the residuals as new values
      dataframe.df$value <- residuals

      return(dataframe.df)
    }

    cellfrac.df <- linRegression(cellfrac.df)
    tfs.df <- linRegression(tfs.df)
    pathways.df <- linRegression(pathways.df)
  }

  MOFAdata = bind_rows(cellfrac.df, pathways.df, tfs.df) %>%
     mutate(group=dataset_id)

  return(MOFAdata)
}


#' Plot heatmap
generateHeatmap <- function(data, title, row_names = F) {
  data <- pivot_wider(data, names_from = "sample", values_from="value", id_cols="feature") %>%
    column_to_rownames("feature")

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
runMOFA <- function(MOFAdata, out_file){

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

  model <- run_mofa(MOFAobject, outfile = out_file, use_basilisk=FALSE)

  set.seed(NULL)

  return(model)

}

#' Plot MOFA results
plotResults <- function(model, title){
  # Plot results

  # Explore explained variance per view and per factor
  ## Values
  print(model@cache$variance_explained$r2_total)
  print(model@cache$variance_explained$r2_per_factor)

  ## Plotting
  # Total explained variance per view

  print(plot_variance_explained(model, x= "view", y ="factor", plot_total = TRUE)[[1]] + ggtitle(title))

  print(plot_variance_explained(model, x= "view", y ="factor", plot_total = TRUE)[[2]] + ggtitle(title))

  if (model@dimensions$G>1){

    ## Plot group results

    # Total explained variance per view
    plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]] + ggtitle(title)

    # Variance per feature per view
    print(plot_variance_explained(model, x="group", y="factor")) + ggtitle(title)}
}


#' Make a complex heatmap showing the different factor weitghts for different datasets.
buildHeatmap <- function(model, title){
  col_fun = circlize::colorRamp2(seq(-1, 1, by = 2/10), rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                              "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                              "#4393C3", "#2166AC", "#053061")))

  setup <- function(weights){
    # Define vectors to split the heatmap
    row.split.vector <- c(rep("Immune cell\n Quantification",10), rep("Pathway\n Scores",14), rep("Transcription Factor Activity", nrow(weights)-24))

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

    return(hm)
  }

  data <- cbind(do.call( "rbind",get_weights(model,  scale = T)))[,1:3]

  # The signs of MOFA factors are meaningless. "Rotate" factors among
  # data sources to make sure the sign is consistent among data sources.
  #
  # We want all factors to have positive correlation with CD8+ T cells,
  # which makes the factors more intuitive to interpret.
  hm1 <- setup(t(t(data)*sign(data["CD8+ T", ])))

  return(hm1)
}

# transformFactors <- function(model){

#   data <- data.frame(do.call("rbind",get_weights(model)))[1:3]
#   colSums(data)

#   if (paste(names(model@dimensions$N),collapse="") == "CRCLUADLUSCSKCM"){
#     colnames(data) <- paste("TCGA",substr(colnames(data),start=7,stop=7),sep=" ")
#   } else if (paste(names(model@dimensions$N),collapse="") == "LUADLUSC"){
#     colnames(data) <- paste("NSCLC",substr(colnames(data),start=7,stop=7),sep=" ")
#   }else{
#     colnames(data) <- paste(paste(names(model@dimensions$N),collapse=""),substr(colnames(data),start=7,stop=7),sep=" ")}
#   return(data)
# }

# #' Plot Correlation matrix
# getCorr <- function(list_of_models,RNA_corr=FALSE){
#   factors <- lapply(list_of_models,transformFactors)

#   cnames <- Reduce(intersect,lapply(factors,rownames))

#   factors <- lapply(factors,function(i) i[cnames,])

#   data <- do.call(cbind,factors)

#   colnames(data) <- unlist(lapply(factors,colnames))

#   cor_matrix <- cor(data)

#   # TODO rotate!
#   # rotation <- getRotation(cor_matrix)
#   # data[rownames(rotation)] <-  t(t(data[rownames(rotation)])*rotation$Sign)

#   # if (RNA_corr == FALSE){
#   #   data['JS B1'] <- get_median_bootstrap_data(jiasharma)
#   #   data['NSCLC B1'] <- get_median_bootstrap_data(nsclc)}

#   cor_matrix_plot <- cor(data)

#   p.mat = cor.mtest(cor_matrix_plot)


#   # cor_matrix_plot[p.mat>0.05] <- 1

#   # p.mat[p.mat>0.05] <- 0.00

#   # # plot <- heatmaply_cor(
#   # #   cor_matrix_plot,
#   # #   node_type = "scatter",
#   # #   column_text_angle = 90,
#   # #   fontsize_row = 13,
#   # #   fontsize_col = 13,
#   # #   point_size_mat = -log10(p.mat),
#   # #   point_size_name = "-log10",
#   # #   label_names = c("x", "y", "Pearson \ncorrelation") )

#   # # print(plot)

#   # cor_matrix_plot[p.mat>0.05] <- 0


#   corrplot(cor_matrix_plot, tl.col="black", tl.srt=75, tl.cex=1.2, p.mat = p.mat, sig.level = 0.001, insig = "blank")


#   return(cor_matrix)

# }
