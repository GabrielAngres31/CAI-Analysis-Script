
#' @name COMPILE.qqGen
#' @description
#' @param df_in
#' @example

COMPILE.qqGen <- function(df_in) {
  progress <- 0
  newprogress <- 0
  numdots <- 50
  barsize <- length(SWITCHBOARD.ACCESSIONLIST)*length(SWITCHBOARD.qqMeasures_VEC)/numdots

  
  for(accession in SWITCHBOARD.ACCESSIONLIST) {
    for(measure in SWITCHBOARD.qqMeasures_VEC) {
      GRAPHER.qqGen(df_in, accession, measure)
      
      newprogress <- newprogress + 1
      if (newprogress %/% barsize > progress) {
        rep_num <- newprogress %/% barsize - progress
        cat(rep(".", rep_num))
        progress <- progress + rep_num
      }
    }
  }
  cat("\n")
}

#' @name COMPILE.intermeasure
#' @description
#' @param model_df_in
#' @example

COMPILE.intermeasure <- function(model_df_in) {
  
  modelnames = names(model_df_in)[-1]
  
  progress <- 0
  newprogress <- 0
  numdots <- 50
  barsize <- length(SWITCHBOARD.ACCESSIONLIST)*length(modelnames)/numdots
  
  for(accession in SWITCHBOARD.ACCESSIONLIST) {
    for (modelname in modelnames) {
      
      UTIL.quickPNG(paste0(accession, "_", modelname, "_INTERMEASURE"), "INTERMEASURE")
      GRAPHER.intermeasure(model_df_in, accession, modelname)
      dev.off()
      
      newprogress <- newprogress + 1
      if (newprogress %/% barsize > progress %/% numdots) {
        rep_num <- newprogress %/% barsize - progress
        cat(rep(".", rep_num))
        progress <- progress + rep_num
      }
    }
  }
  cat("\n")
}

# COMPILE.multiplicative <- function() {
#   NULL
# }
# COMPILE.elliptical <- function() {
#   NULL
# }

COMPILE.boxplots <- function() {
  NULL
}

#' @name COMPILE.tukeygraphs
#' @description
#' @param dataset
#' @example

COMPILE.tukeygraphs <- function(dataset) {
  measures <- c("width", "height", "diameter", "thickness", "fresh_weight")
  
  for (measure in measures) {
    UTIL.quickPNG(measure, "Tukey_Graphs\\Tukey_Intercomparisons")
    CALC.tukeyTester(dataset, measure) %T>%
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

COMPILE.R2heatmap <- function(R2_df_in, filenamestring, subfolder, abbreviations) {
  UTIL.quickPNG(filenamestring, subfolder)
  GRAPHER.R2heatmap(R2_df_in, abbreviations)
  dev.off()
}

COMPILE.SBCheatmap <- function(SBC_df_in, filenamestring, subfolder, abbreviations) {
  UTIL.quickPNG(filenamestring, subfolder)
  GRAPHER.SBCheatmap(SBC_df_in, abbreviations)
  dev.off()
}