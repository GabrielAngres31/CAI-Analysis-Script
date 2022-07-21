#' @name GRAPHER.tukeyPlotter
#' @description 
#' @param tukey_results
#' @example


GRAPHER.tukeyPlotter <- function(tukey_results) {
  tukeyplot <- plot(tukey_results, las=1 , col="gray29")
  show(tukeyplot)
  #return(tukeyplot)
}

#' @name GRAPHER.tukeyHypergraph
#' @description
#' @param tukey_edgelist
#' @param measure
#' @example

GRAPHER.tukeyHypergraph <- function(tukey_edgelist, measure) {
  plot(hypergraph_from_edgelist(tukey_edgelist))
  text(0, 1.25, cex = 1.75, paste0("Hypergraph of ", measure))  
}
#------------

#]]DEBUG
#GRAPHER.tukeyBoxplot(SWITCHBOARD.csvAUGMFILE, "height")

#-----------------

#' @name GRAPHER.qqGen
#' @description
#' @param df_in
#' @param accession
#' @param measure
#' @example

GRAPHER.qqGen <- function(df_in, accession, measure) {
  target_data <- df_in[[measure]][df_in$accession == accession]
  UTIL.quickPNG(paste0(measure, "--", accession, "_LIN"), "QQPlot_Images")
  qqPlot(target_data, main = paste0(accession, ", ", measure), envelope = 0.95)
  dev.off()
}

#' @name GRAPHER.intermeasure
#' @description
#' @param fits_df
#' @param accession
#' @param modelname
#' @example

GRAPHER.intermeasure <- function(fits_df, accession, modelname) {
  #GET x_measure and y_measure from model object
  model_in <- CALC.ModelFetch(fits_df, accession, modelname)
  
  model_in <- UTIL.AugmentModelWithLeverageTagging(model_in)

  model_in <- UTIL.AugmentModelWithID(model_in, accession)
  
  colnames_data <- names(model_in$model)
  model_df <- select(model_in$model, c(colnames_data))
  y_measure <- colnames_data[1]
  x_measure <- colnames_data[2]
  slope <- model_in$coefficients[2]
  intercept <- model_in$coefficients[1]
 
  adj_rsq <- round(summary(model_in)[["adj.r.squared"]], 3)

  intermeasure_plot_out <- ggplot(model_df, aes(!!as.symbol(x_measure), !!as.symbol(y_measure),
                                                label = ifelse(tagged == "N/A", "", ID))) +
    geom_point(aes(color = tagged)) +
    labs(title = (paste0("INTERMEASURE\nMODEL:", modelname, "\nACCESSION:", accession, "\nADJ.R^2:", adj_rsq))) +
    geom_text(hjust=1.05, vjust=0) +
    geom_text(hjust=-.20, vjust=0) +
    geom_abline(intercept = intercept, slope = slope)
  plot(intermeasure_plot_out)
}

#]]DEBUG
#GRAPHER.intermeasure(test_aug_filt, DATA_OBJECTS.LinearRegressions_DF, 242, "height~width")


GRAPHER.R2heatmap <- function(R2_df_in, abbreviations) {
  
  R2_heat_palette <- colorRampPalette(brewer.pal(8, "RdYlGn"))(25)
  
  R2_df_in %>%
    setNames(c("accession", abbreviations)) %>%
    column_to_rownames("accession") %>%
    data.matrix() %>%
    heatmap(Rowv = NA, Colv = NA, col = R2_heat_palette)
}

#]]DEBUG
#GRAPHER.R2heatmap(test_heatmap_in, "Test_Heatmap_R2", "Heatmaps")

#' @name GRAPHER.SBCheatmap
#' @description
#' @param 
#' @example

  
GRAPHER.SBCheatmap <- function(SBC_df_in, abbreviations) {

  SBC_heat_palette <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(25))
  
  SBC_df_in %>%
    setNames(c("accession", abbreviations)) %>%
    column_to_rownames("accession") %>%
    data.matrix() %>%
    heatmap(Rowv = NA, Colv = NA, scale = "row", col = SBC_heat_palette)
}

#]]DEBUG
#GRAPHER.SBCheatmap(test_heatmap_in_SBC, "Test_Heatmap_SBC", "Heatmaps")