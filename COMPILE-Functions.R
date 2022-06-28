
#' @name COMPILE.qqGen
#' @description
#' @param df_in
#' @example

COMPILE.qqGen <- function(df_in) {
  for(accession in SWITCHBOARD.ACCESSIONLIST) {
    for(measure in SWITCHBOARD.qqMeasures_VEC) {
      GRAPHER.qqGen(df_in, accession, measure)
    }
  }
}

#' @name COMPILE.intermeasure
#' @description
#' @param model_df_in
#' @example

COMPILE.intermeasure <- function(model_df_in) {
  modelnames = names(model_df_in)[-1]
  for(accession in SWITCHBOARD.ACCESSIONLIST) {
    for (modelname in modelnames) {
      UTIL.quickPNG(paste0(accession, "_", modelname, "_INTERMEASURE"), "INTERMEASURE")
      GRAPHER.intermeasure(model_df_in, accession, modelname)
      dev.off()
    }
  }
}

COMPILE.multiplicative <- function() {
  NULL
}
COMPILE.elliptical <- function() {
  NULL
}

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