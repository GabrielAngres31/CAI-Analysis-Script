cat("GENERATING PLOTS AND TABLES\n")

# Quantile-Quantile Plot Generation ---------------------------------------

cat("Quantile-Quantile Plots\n")
COMPILE.qqGen(SWITCHBOARD.csvMAINFILE, SWITCHBOARD.ACCESSIONLIST, SWITCHBOARD.initialMeasures_VEC)

# Intermeasure Analysis Figure Generation ---------------------------------

cat("Intermeasure Plots\n")
COMPILE.intermeasure(DATA_OBJECTS.LinearRegressions_DF)

cat("Violin Plots\n")
COMPILE.ViolinPlot(SWITCHBOARD.csvMAINFILE, SWITCHBOARD.initialMeasures_VEC, "Unfiltered")
COMPILE.ViolinPlot(DATA_OBJECTS.csvAUGM_3PN, SWITCHBOARD.initialMeasures_VEC, "Filtered")

cat("Tukey Intercomparisons and Grouping Hypergraphs\n")
COMPILE.tukeygraphs(SWITCHBOARD.csvAUGMFILE)


# Heatmap Generation ------------------------------------------------------

cat("HEATMAP: Intermeasure Comparisons\n")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_intermeasure, "Lin_UNF",   "Heatmaps", "Intermeasure Models\n", SWITCHBOARD.models_intermeasure_abbr)


cat("HEATMAP: Multiplicative Model Comparisons - Unfiltered\n")
cat("With Interactions...")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_multiplicative.splitinteraction[[2]], "Mult_UNF_mult",  "Heatmaps", "Multiplicative Models\nUnfiltered Dataset with Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[2]]$models)
cat("Without Interactions...")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_multiplicative.splitinteraction[[1]], "Mult_UNF_plus",  "Heatmaps", "Multiplicative Models\nUnfiltered Dataset with No Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[1]]$models)
cat("Altogether...")
cat("\n")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative, "Mult_UNF",  "Heatmaps", "Multiplicative Models\nUnfiltered Dataset",SWITCHBOARD.models_multiplicative_abbr)


cat("HEATMAP: Multiplicative Model Comparisons - FILTERED on 3P/N Leverage\n")
cat("With Interactions...")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_multiplicative.splitinteraction[[2]], "Mult_3PN_mult",  "Heatmaps", "Multiplicative Models\n3P/N Dataset with Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[2]]$models)
cat("Without Interactions...")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_multiplicative.splitinteraction[[1]], "Mult_3PN_plus",  "Heatmaps", "Multiplicative Models\n3P/N Dataset with No Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[1]]$models)
cat("Altogether...")
cat("\n")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_multiplicative, "Mult_3PN",  "Heatmaps", "Multiplicative Models\n3P/N Dataset", SWITCHBOARD.models_multiplicative_abbr)

cat("HEATMAP: Multiplicative Model Comparisons - Unfiltered - Alternate Formula\n")
cat("With Interactions...")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_multiplicative.splitinteraction_RATIO[[2]], "Mult_UNF_mult_ALT",  "Heatmaps", "Multiplicative Models ALT\nUnfiltered Dataset with Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction_RATIO[[2]]$models)
cat("Without Interactions...")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_multiplicative.splitinteraction_RATIO[[1]], "Mult_UNF_plus_ALT",  "Heatmaps", "Multiplicative Models ALT\nUnfiltered Dataset with No Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction_RATIO[[1]]$models)
cat("Altogether...")
cat("\n")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative_RATIO, "Mult_UNF_ALT",  "Heatmaps", "Multiplicative Models ALT\nUnfiltered Dataset",SWITCHBOARD.models_multiplicative_abbr_RATIO)


cat("HEATMAP: Multiplicative Model Comparisons - FILTERED on 3P/N Leverage - Alternate Formula\n")
cat("With Interactions...")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_multiplicative.splitinteraction_RATIO[[2]], "Mult_3PN_mult_ALT",  "Heatmaps", "Multiplicative Model ALTs\n3P/N Dataset with Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction_RATIO[[2]]$models)
cat("Without Interactions...")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_multiplicative.splitinteraction_RATIO[[1]], "Mult_3PN_plus_ALT",  "Heatmaps", "Multiplicative Models ALT\n3P/N Dataset with No Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction_RATIO[[1]]$models)
cat("Altogether...")
cat("\n")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_multiplicative_RATIO, "Mult_3PN_ALT",  "Heatmaps", "Multiplicative Models ALT\n3P/N Dataset", SWITCHBOARD.models_multiplicative_abbr_RATIO)


cat("HEATMAP: Elliptical Model Comparisons - Unfiltered\n")
cat("With Interactions...")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_elliptical.splitinteraction[[2]], "Ellip_UNF_mult",  "Heatmaps", "Elliptical Models\nUnfiltered Dataset with Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[2]]$models)
cat("Without Interactions...")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_elliptical.splitinteraction[[1]], "Ellip_UNF_plus",  "Heatmaps", "Elliptical Models\nUnfiltered Dataset with No Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[1]]$models)
cat("Altogether...")
cat("\n")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_elliptical, "Ellip_UNF", "Heatmaps", "Elliptical Models\nUnfiltered Dataset", SWITCHBOARD.models_elliptical_abbr)


cat("HEATMAP: Elliptical Model Comparisons - FILTERED on 3P/N Leverage\n")
cat("With Interactions...")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_elliptical.splitinteraction[[2]], "Ellip_3PN_mult",  "Heatmaps", "Elliptical Models\n3P/N Dataset with Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[2]]$models)
cat("Without Interactions...")
  UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_elliptical.splitinteraction[[1]], "Ellip_3PN_plus",  "Heatmaps", "Elliptical Models\n3P/N Dataset with No Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[1]]$models)
cat("Altogether...")
cat("\n")
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_elliptical, "Ellip_3PN", "Heatmaps", "Elliptical Models\n3P/N Data", SWITCHBOARD.models_elliptical_abbr)


cat("HEATMAP: Reis et. al. Result Comparisons\n")
cat("Fresh Weight...")
CALC.CompareREIS_FRESHWEIGHT()
cat("Area...")
CALC.CompareREIS_AREA()
cat("Dry Weight...")
cat("\n")
CALC.CompareREIS_DRYWEIGHT()
cat("\n")


cat("--------------------------------\n")
cat("--------------------------------\n")
