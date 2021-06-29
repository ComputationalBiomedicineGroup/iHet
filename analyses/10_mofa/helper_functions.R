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

#' Rotate factors such that they are positively associated with CD8+ T-cell infiltration.
get_rotation <- function(weights) {
  diag(sign(weights[["Immune cells quantification"]]["CD8 T", ]))
}

#' extract and preprocess weights from a MOFA model
get_weight_df <- function(model) {
  weights <- get_weights(model, scale = TRUE)
  lapply(names(weights), function(modality) {
    weights[[modality]] %*% get_rotation(weights) %>%
      as_tibble(rownames = "feature") %>%
      mutate(view = modality) %>%
      rename_with(function(x) {
        str_replace(x, "V", "factor ")
      }, starts_with("V"))
  }) %>%
    bind_rows() %>%
    arrange(view, feature)
}

#' extract factors from a MOFA model
get_factor_df <- function(model) {
  weights <- get_weights(model, scale = TRUE)
  get_factors(model) %>%
    do.call("rbind", .) %*% get_rotation(weights) %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample") %>%
    rename_with(function(x) {
      str_replace(x, "V", "factor ")
    }, starts_with("V"))
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
