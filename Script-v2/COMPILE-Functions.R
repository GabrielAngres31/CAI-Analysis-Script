# qqGen and Intermeasure Plots --------------------------------------------

#' @name COMPILE.qqGen
#' @description Generates and stores 
#' all Quantile-Quantile plot for measures in a dataframe as .png files.
#' @param data_df The table of measure values by individual pads.
#' @param accessions A vector of accessions IDs.
#' @param measures A vector of measures to generate plots for.

COMPILE.qqGen <- function(data_df, accessions, measures) {

  for(accession in accessions) {
    for(measure in measures) {
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
    cat(paste0("For Accession ", accession, "\n"))
    for (modelname in modelnames) {
      cat(paste0(modelname, "..."))
      UTIL.quickPNG(paste0(accession, "_", modelname, "_INTERMEASURE"), "INTERMEASURE")
      GRAPHER.intermeasure(model_df, accession, modelname)
      dev.off()
    }
    cat("\n")
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

#' @name COMPILE.ViolinPlot
#' @description Generates and stores a Violin plot of the range of values for a particular set of measures across multiple accessions.
#' @param data_df The table of measure values by individual pads.
#' @param measurelist The list of measures to iterate over, as a character vector.
#' @param data_title The title of the heatmap plot, as a string.

COMPILE.ViolinPlot <- function(data_df, measurelist, data_title) {
  for (measure in measurelist) {
    UTIL.quickPNG(paste0(measure, " -  ", data_title), "Violin_Graphs\\Graphs")
    GRAPHER.ViolinPlot(data_df, measure)
    dev.off()
  }
}