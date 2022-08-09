# qqGen and Intermeasure Plots --------------------------------------------

#' @name COMPILE.qqGen
#' @description Generates and stores all Quantile-Quantile plot for measures in a dataframe as .png files.
#' @param data_df The table of measure values by individual pads.

COMPILE.qqGen <- function(data_df) {

  for(accession in SWITCHBOARD.ACCESSIONLIST) {
    for(measure in SWITCHBOARD.qqMeasures_VEC) {
      GRAPHER.qqGen(data_df, accession, measure)
    }
  }
}

#' @name COMPILE.intermeasure
#' @description Generates and stores all intermeasure plots for all measures and accessions in a dataframe as .png files.
#' @param model_df The list of models to plot, by name.

COMPILE.intermeasure <- function(model_df) {
  
  modelnames = names(model_df)[-1]

  
  for(accession in SWITCHBOARD.ACCESSIONLIST) {
    for (modelname in modelnames) {
      
      UTIL.quickPNG(paste0(accession, "_", modelname, "_INTERMEASURE"), "INTERMEASURE")
      GRAPHER.intermeasure(model_df, accession, modelname)
      dev.off()
    }
  }
}

# Tukey Mean Test Graphs --------------------------------------------------

#' @name COMPILE.tukeygraphs
#' @description Generates all mean intercomparison Tukey tests for a given dataset, comparing between accessions.
#' @param data_df The table of measure values by individual pads.

COMPILE.tukeygraphs <- function(data_df) {
  measures <- c("width", "height", "diameter", "thickness", "fresh_weight")
  
  for (measure in measures) {
    UTIL.quickPNG(measure, "Tukey_Graphs\\Tukey_Intercomparisons")
    CALC.tukeyTester(data_df, measure) %T>%
      GRAPHER.tukeyPlotter() ->  tukey_results
    dev.off()
    
    UTIL.quickPNG(measure, "Tukey_Graphs\\Tukey_Hypergraphs")
    tukey_edgelist <- tukey_results %>%
      CALC.tukeyLabel_DF() %>%
      CALC.tukeyRelation() %>%
      CALC.tukeyEdgelist()
    GRAPHER.tukeyHypergraph(tukey_edgelist, measure)
    dev.off()
  }
} 


# Summary Statistic Graphs ------------------------------------------------

#' @name COMPILE.R2heatmap
#' @description Generates and stores an adj. R^2 Heatmap for a dataframe of values, by model and accession, as a .png file.
#' @param R2_df_in The dataframe of SBC values, with model columns and accession rows.
#' @param filenamestring The name given for the .png file.
#' @param subfolder The folder that the file should be stored in.
#' @param abbreviations A vector of string abbreviations for the models for intelligible graphing.

COMPILE.R2heatmap <- function(R2_df_in, filenamestring, subfolder, heatmap_title, abbreviations) {
  UTIL.quickPNG(filenamestring, subfolder)
  GRAPHER.R2heatmap(R2_df_in, heatmap_title, abbreviations)
  dev.off()
}

#' @name COMPILE.SBCheatmap
#' @description Generates and stores a SBC Heatmap for a dataframe of values, by model and accession, as a .png file.
#' @param SBC_df_in The dataframe of SBC values, with model columns and accession rows.
#' @param filenamestring The name given for the .png file.
#' @param subfolder The folder that the file should be stored in.
#' @param abbreviations A vector of string abbreviations for the models for intelligible graphing.

COMPILE.SBCheatmap <- function(SBC_df_in, filenamestring, subfolder, heatmap_title, abbreviations) {
  UTIL.quickPNG(filenamestring, subfolder)
  GRAPHER.SBCheatmap(SBC_df_in, heatmap_title, abbreviations)
  dev.off()
}