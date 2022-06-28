#' @name GRAPHER.tukeyPlotter
#' @description
#' @param tukey_results
#' @example


GRAPHER.tukeyPlotter <- function(tukey_results) {
  tukeyplot <- plot(tukey_results, las=1 , col="gray29")
  show(tukeyplot)
  #return(tukeyplot)
}

#' @name 
#' @description
#' @param 
#' @example

GRAPHER.tukeyHypergraph <- function(tukey_edgelist, measure) {
  plot(hypergraph_from_edgelist(tukey_edgelist))
  text(0, 1.25, cex = 1.75, paste0("Hypergraph of ", measure))  
}
#------------

#]]DEBUG
#GRAPHER.tukeyBoxplot(SWITCHBOARD.csvAUGMFILE, "height")

#-----------------

#' @name 
#' @description
#' @param 
#' @example

GRAPHER.qqGen <- function(df_in, accession, measure) {
  target_data <- df_in[[measure]][df_in$accession == accession]
  UTIL.quickPNG(paste0(measure, "--", accession, "_LIN"), "QQPlot_Images")
  qqPlot(target_data, main = paste0(accession, ", ", measure), envelope = 0.95)
  dev.off()
}

#' @name 
#' @description
#' @param 
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
  
  #print(model_df)
  #print(CALC.PadIDfromModelValuePairs(model_df, data_df))
  #print(model_df)
  
  # model_df <- mutate(model_df, tagged = "N/A")
  # model_df$tagged[indices_MID] <- ">2p/n Threshold"
  # model_df$tagged[indices_HI] <- ">3p/n Threshold"

  intermeasure_plot_out <- ggplot(model_df, aes(!!as.symbol(x_measure), !!as.symbol(y_measure),
                                                label = ifelse(tagged == "N/A", "", ID))) +
    geom_point(aes(color = tagged)) +
    labs(title = (paste0("INTERMEASURE\nMODEL:", modelname, "\nACCESSION:", accession))) +
    geom_text(hjust=1.05, vjust=0) +
    geom_text(hjust=-.20, vjust=0) +
    geom_abline(intercept = intercept, slope = slope)
  plot(intermeasure_plot_out)
}

#]]DEBUG
#GRAPHER.intermeasure(test_aug_filt, DATA_OBJECTS.LinearRegressions_DF, 242, "height~width")


#' @name 
#' @description
#' @param 
#' @example

GRAPHER.multiplicative <- function(fits_df, accession, modelname) {
  model_in <- CALC.ModelFetch(fits_df, accession, modelname)
  
}

GRAPHER.elliptical <- function() {
  NULL
}

#' @name 
#' @description
#' @param 
#' @example

GRAPHER.R2heatmap <- function(R2_df_in, filenamestring, subfolder) {
  UTIL.quickPNG(filenamestring, subfolder)
  R2_heat_palette <- colorRampPalette(brewer.pal(8, "RdYlGn"))(25)
  
  R2_df_in %>%
    column_to_rownames("accession") %>%
    data.matrix() %>%
    heatmap(Rowv = NA, Colv = NA, col = R2_heat_palette)
  
  dev.off()
}

#]]DEBUG
#GRAPHER.R2heatmap(test_heatmap_in, "Test_Heatmap_R2", "Heatmaps")

#' @name 
#' @description
#' @param 
#' @example

GRAPHER.SBCheatmap <- function(SBC_df_in, filenamestring, subfolder) {
  UTIL.quickPNG(filenamestring, subfolder)
  
  SBC_heat_palette <- colorRampPalette(brewer.pal(8, "RdBu"))(25)
  
  SBC_df_in %>%
    column_to_rownames("accession") %>%
    data.matrix() %>%
    heatmap(Rowv = NA, Colv = NA, scale = "column", col = SBC_heat_palette)
  
  dev.off()
}

#]]DEBUG
#GRAPHER.SBCheatmap(test_heatmap_in_SBC, "Test_Heatmap_SBC", "Heatmaps")