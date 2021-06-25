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
library(rstatix)
library(stringr)
library(ggbeeswarm)
library(BayesFactor)
library(ggrepel)
library(readr)


rdbu10_palette <- c(
  "#67001F", "#B2182B", "#D6604D", "#F4A582",
  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
  "#4393C3", "#2166AC", "#053061"
)


#' Function to compute the iHet score starting from immune cell, pathway, and TF features
#' using either the NSCLC or JiaSharma feature weights
#' TODO: fix path to weightfile
getiHet <- function(featobj, dataset, weightfile = "/Users/francescafinotello/Dropbox/Research_projects/MODiSH/RData/median_weights.rds", model = c("NSCLC", "JiaSharma")) {

  # Select model weights
  model <- match.arg(model, c("NSCLC", "JiaSharma"))

  # Derive the median weights for factor 1
  weights <- readRDS(weightfile)
  if (model == "JiaSharma") {
    F1 <- weights$JiaSharma$`factor 1`
    names(F1) <- weights$JiaSharma$feature
  } else if (model == "NSCLC") {
    F1 <- weights$NSCLC$`factor 1`
    names(F1) <- weights$NSCLC$feature
  }

  # Extract, aggregate, and normalize the input features
  immcell <- featobj$cellfrac[[dataset]]
  immcell <- immcell[which(rownames(immcell) != "Other"), ]
  immcell[, "CD4 T"] <- immcell[, "CD4 T"] + immcell[, "Treg"]
  immcell <- log10(immcell * 100 + 0.001)
  tf <- featobj$tf[[dataset]]
  path <- featobj$pathway[[dataset]]
  path <- t(scale(t(path)))
  features <- cbind(immcell, tf, path)

  # Select the common features
  # TODO: put some warnings for missing/exceeding ones?
  cfeatures <- intersect(colnames(features), names(F1))
  features <- features[, cfeatures]
  F1 <- F1[cfeatures]
  features <- t(features)

  # Compute iHet on the input features
  iHet <- F1 %*% features
  iHet <- iHet[1, , drop = TRUE]

  return(iHet)
}

#' Plot MOFA results
plotResults <- function(model, title) {
  # Plot results

  # Explore explained variance per view and per factor
  ## Values
  print(model@cache$variance_explained$r2_total)
  print(model@cache$variance_explained$r2_per_factor)

  ## Plotting
  # Total explained variance per view

  print(plot_variance_explained(model, x = "view", y = "factor", plot_total = TRUE)[[1]] + ggtitle(title))

  print(plot_variance_explained(model, x = "view", y = "factor", plot_total = TRUE)[[2]] + ggtitle(title))

  if (model@dimensions$G > 1) {

    ## Plot group results

    # Total explained variance per view
    plot_variance_explained(model, x = "group", y = "factor", plot_total = T)[[2]] + ggtitle(title)

    # Variance per feature per view
    print(plot_variance_explained(model, x = "group", y = "factor")) + ggtitle(title)
  }
}


#' Make a complex heatmap showing the different factor weitghts for different datasets.
buildHeatmap <- function(weight_df, title) {
  col_fun <- circlize::colorRamp2(seq(-1, 1, by = 2 / 10), rev(rdbu10_palette))
  # Define vectors to split the heatmap
  row.split.vector <- recode(weight_df$view, "Immune cells quantification" = "Immune\ncells", "Pathway scores" = "Pathways", "Transcription factors" = "Transcription\nfactors")
  weights <- weight_df %>%
    select(-view) %>%
    column_to_rownames("feature") %>%
    as.matrix() %>%
    .[, 1:3]

  # Size of the each cell
  setsize <- 0.30 # cm
  hm <- Heatmap(weights,
    name = "weights",
    rect_gp = gpar(col = "grey", lwd = 0.5),
    col = col_fun, border = T,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = T, show_column_names = T, row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    row_split = row.split.vector,
    row_title_gp = gpar(fontsize = 15, fontface = 2),
    row_title_rot = 0,
    column_names_gp = gpar(fontsize = 9),
    column_title_gp = gpar(fontsize = 9, fontface = 2),
    column_title = title,
    width = unit(setsize, "cm") * ncol(weights), height = unit(setsize, "cm") * nrow(weights)
  )

  return(hm)
}
