#Load relevant packages into library

cat("Loading Packages...\n")
library("dplyr")
library("tibble")
library("rsq")
library("car")
library("Rcmdr")
library("hash")
library("multcompView")
library("HyperG")
library("plotly")
library("htmlwidgets")
library("magrittr")
library("tidyr")
library("multiApply")
library("stringr")
library("viridis")
library("pheatmap")
library("RColorBrewer")
library("data.table")

# SWITCHBOARD Setup -------------------------------------------------------

cat("Setting up SWITCHBOARD: \n")

#NOTE: ALL DATAFILES MUST BE LOADED INTO THE WORKING DIRECTORY PRIOR TO THE PROGRAM START
SWITCHBOARD.DIRECTORY <- getwd()

cat("::| Loading Dataset of Interest...\n")
SWITCHBOARD.csvMAINFILE <- read.csv(file = DATA_FILE)

# This section hard-codes the dry_weight data collected in the experiment for addition into the data set.

### USE THIS one as agument file 
# https://community.rstudio.com/t/recognize-usr-bin-env-rscript-header-as-r-script/59237
SWITCHBOARD.FRESH_DRY_RAWS <- tribble(
  ~accession, ~mean_freshweight, ~SEM_mean_fresh, ~mean_dryweight, ~SEM_mean_dry,
  242, 9.17, 0.94, 0.6, 0.07,
  246, 8.48, 0.94, 0.49, 0.07,
  319, 8.73, 0.94, 0.6, 0.07,
  325, 10.23, 1, 0.68, 0.07,
  326, 5.76, 0.94, 0.49, 0.07,
  390, 10.65, 0.97, 0.89, 0.07,
  572, 10.37, 0.94, 0.65, 0.07,
  580, 5.47, 0.91, 0.4, 0.07,
  582, 8.11, 0.94, 0.57, 0.07,
  584, 9.71, 1.13, 0.66, 0.08,
  585, 4.56, 0.95, 0.32, 0.07,
  839, 3.97, 0.94, 0.33, 0.07,
  845, 7.59, 0.94, 0.71, 0.07,
  854, 2.31, 0.94, 0.21, 0.07
) 

cat("::| Setting Dry Weights...\n")
SWITCHBOARD.DRY_PROPORTIONS <- SWITCHBOARD.FRESH_DRY_RAWS %>%
  mutate(proportion = mean_dryweight/mean_freshweight) %>%
  select(accession, proportion)

cat("::| Setting Accession List...\n")
## This one could be read from the file load
SWITCHBOARD.ACCESSIONLIST <- c(242, 246, 319, 325, 326, 390, 572, 580, 582, 584, 585, 839, 845, 854)

## Save the CSV file
cat("::| Setting Derived Measures...\n")
SWITCHBOARD.csvAUGMFILE <- SWITCHBOARD.csvMAINFILE %>%
  mutate(H_div_W = height / width,
         FW_div_W = fresh_weight / width,
         FW_div_D = fresh_weight / diameter,
         FW_div_H = fresh_weight / height,
         FW_div_T = fresh_weight / thickness,
         D_div_W = diameter / width) %>% 
  merge(read.csv('Pad_Area_Estimations.csv', fileEncoding = 'UTF-8-BOM'), by = "completeID") %>%
  mutate(Theo_Area = pi*height*width*0.25) %>%
  mutate(PartRatio = ((height-width)/(height+width))^2) %>%
  mutate(Pade_Peri = pi*(height+width)*(64-3*PartRatio^2)/(64-16*PartRatio)) %>%
  mutate(Pade_Derived_Diam = Pade_Peri/pi) %>%
  full_join(SWITCHBOARD.DRY_PROPORTIONS, by = "accession") %>%
  mutate(dry_weight = proportion*fresh_weight)

cat("::| Setting Model Guides...\n")
DATA_OBJECTS.MODELS_intermeasure <- tribble(
  ~ xname, ~ yname,
  'width',     'fresh_weight',
  'height',    'fresh_weight',
  'diameter',  'fresh_weight',
  'thickness', 'fresh_weight',
  'width',     'height',
  'width',     'diameter',
  'width',     'thickness',
  'height',    'diameter',
  'height',    'thickness',
  'diameter',  'thickness',
)

cat("::| Setting Intermeasure Models\n")
SWITCHBOARD.models_intermeasure_abbr <- c("FW~W", "FW~H", "FW~D", "FW~T", "H~W", "D~W", "T~W", "D~H", "T~H", "T~D")

cat("::| Setting Initial Measures\n")
SWITCHBOARD.initialMeasures_VEC <- c("width", "diameter", "height", "thickness", "fresh_weight")

cat("::| Setting Multiplicative Models for Rectangular Model Calculation\n")
SWITCHBOARD.models_multiplicative <- c(
  "fresh_weight~width", "fresh_weight~height", "fresh_weight~diameter","fresh_weight~thickness",
  "fresh_weight~D_div_W", 
  "fresh_weight~width+height", "fresh_weight~width*height",
  "fresh_weight~width+diameter", "fresh_weight~width*diameter",
  "fresh_weight~width+thickness", "fresh_weight~width*thickness",
  "fresh_weight~height+diameter", "fresh_weight~height*diameter",
  "fresh_weight~height+thickness", "fresh_weight~height*thickness",
  "fresh_weight~diameter+thickness", "fresh_weight~diameter*thickness",
  "fresh_weight~width+height+diameter", "fresh_weight~width*height*diameter", 
  "fresh_weight~width+height+thickness", "fresh_weight~width*height*thickness",
  "fresh_weight~width+diameter+thickness", "fresh_weight~width*diameter*thickness",
  "fresh_weight~height+diameter+thickness", "fresh_weight~height*diameter*thickness",
  "fresh_weight~width+height+diameter+thickness", "fresh_weight~width*height*diameter*thickness"
)

cat("::| Setting Multiplicative Models for Elliptical Model Calculation\n")
SWITCHBOARD.models_multiplicative_RATIO <- c(
  "fresh_weight~width", "fresh_weight~height", "fresh_weight~diameter","fresh_weight~thickness",
  "fresh_weight~D_div_W", 
  "fresh_weight~width+height", "fresh_weight~width*height",
  "fresh_weight~width+D_div_W", "fresh_weight~width*D_div_W",
  "fresh_weight~width+thickness", "fresh_weight~width*thickness",
  "fresh_weight~height+D_div_W", "fresh_weight~height*D_div_W",
  "fresh_weight~height+thickness", "fresh_weight~height*thickness",
  "fresh_weight~D_div_W+thickness", "fresh_weight~D_div_W*thickness",
  "fresh_weight~width+height+D_div_W", "fresh_weight~width*height*D_div_W", 
  "fresh_weight~width+height+thickness", "fresh_weight~width*height*thickness",
  "fresh_weight~width+D_div_W+thickness", "fresh_weight~width*D_div_W*thickness",
  "fresh_weight~height+D_div_W+thickness", "fresh_weight~height*D_div_W*thickness",
  "fresh_weight~width+height+D_div_W+thickness", "fresh_weight~width*height*D_div_W*thickness"
)

cat("::| Setting Multiplicative Abbreviations\n")
SWITCHBOARD.models_multiplicative_abbr <- c(
  "W", "H", "D", "T",
  "D/W",
  "W+H", "W*H",
  "W+D", "W*D",
  "W+T", "W*T",
  "H+D", "H*D",
  "H+T", "H*T",
  "D+T", "D*T",
  "W+H+D", "W*H*D",
  "W+H+T", "W*H*T",
  "W+D+T", "W*D*T",
  "H+D+T", "H*D*T",
  "W+H+D+T", "W*H*D*T"
  #"I(width*height*diameter*thickness)"
)

SWITCHBOARD.models_multiplicative_abbr_RATIO <- c(
  "W", "H", "D", "T",
  "D/W (r)",
  "W+H", "W*H",
  "W+r", "W*r",
  "W+T", "W*T",
  "H+r", "H*r",
  "H+T", "H*T",
  "T+r", "T*r",
  "W+H+r", "W*H*r",
  "W+H+T", "W*H*T",
  "W+T+r", "W*T*r",
  "H+T+r", "H*T*r",
  "W+H+T+r", "W*H*T*r"
)

cat("::| Setting Elliptical Abbreviations\n")
SWITCHBOARD.models_elliptical <- c(
  "fresh_weight~Area",
  "fresh_weight~Theo_Area",
  "fresh_weight~Area+thickness",
  "fresh_weight~Area*thickness",
  "fresh_weight~Theo_Area+thickness",
  "fresh_weight~Theo_Area*thickness",
  "fresh_weight~Theo_Area+thickness+Pade_Derived_Diam",
  "fresh_weight~Theo_Area*thickness*Pade_Derived_Diam"
)

SWITCHBOARD.models_elliptical_abbr <- c(
  "A", "Th_A",
  "A+T", "A*T",
  "Th_A+T", "Th_A*T",
  "Th_A+T+PD", "Th_A*T*PD"
)

#---
cat("::| Setting Reis et. al. Comparison Models\n")
SWITCHBOARD.models_to_compare.REIS.AREA <- c(
  "Area~height*width",
  "Area~Theo_Area",
  "Area~Theo_Area*diameter"
)

SWITCHBOARD.models_to_compare.REIS.AREA_abbr <- c(
  "A~HW",
  "A~A_th",
  "A~A_th*D"
)

SWITCHBOARD.models_to_compare.REIS.DRYWEIGHT <- c(
  "dry_weight~width*thickness",
  "dry_weight~width*thickness*Pade_Derived_Diam",
  "dry_weight~Theo_Area*thickness*Pade_Derived_Diam"
)

SWITCHBOARD.models_to_compare.REIS.DRYWEIGHT_abbr <- c(
  "DW~W*T",
  "DW~W*T*pdD",
  "DW~A_th*T*pdD"
)

#SWITCHBOARD.models_to_compare.REIS.DRY_WEIGHT

# DATA_OBJECT Generation: Unfiltered --------------------------------------------------
### https://stackoverflow.com/questions/16714020/loop-through-data-frame-and-variable-names


cat("::| Polishing Model Guides...\n")
DATA_OBJECTS.MODELS_intermeasure$models <- apply(DATA_OBJECTS.MODELS_intermeasure[, c("yname", "xname")], 1, paste, collapse = "~") 

DATA_OBJECTS.MODELS_multiplicative <- as.data.frame(SWITCHBOARD.models_multiplicative) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_multiplicative_abbr <- as.data.frame(SWITCHBOARD.models_multiplicative_abbr) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_multiplicative_RATIO <- as.data.frame(SWITCHBOARD.models_multiplicative_RATIO) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_multiplicative_abbr_RATIO <- as.data.frame(SWITCHBOARD.models_multiplicative_abbr_RATIO) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_elliptical <- as.data.frame(SWITCHBOARD.models_elliptical) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_elliptical_abbr <- as.data.frame(SWITCHBOARD.models_elliptical_abbr) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_REIS_AREA <- as.data.frame(SWITCHBOARD.models_to_compare.REIS.AREA) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_REIS_AREA_abbr <- as.data.frame(SWITCHBOARD.models_to_compare.REIS.AREA_abbr) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_REIS_DRYWEIGHT <- as.data.frame(SWITCHBOARD.models_to_compare.REIS.DRYWEIGHT) %>%
  set_colnames(c("models"))
  
DATA_OBJECTS.MODELS_REIS_DRYWEIGHT_abbr <- as.data.frame(SWITCHBOARD.models_to_compare.REIS.DRYWEIGHT_abbr) %>%
  set_colnames(c("models"))

DATA_OBJECTS.models_multiplicative.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative, 4)
DATA_OBJECTS.models_multiplicative.splitinteraction_RATIO <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_RATIO, 4)
DATA_OBJECTS.models_elliptical.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical, 2)

DATA_OBJECTS.models_multiplicative_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_abbr, 4)
DATA_OBJECTS.models_multiplicative_abbr.splitinteraction_RATIO <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_abbr_RATIO, 4)
DATA_OBJECTS.models_elliptical_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical_abbr, 2)

cat("Generating Model Dataframes\n")

### loop 
DATA_OBJECTS.LinearRegressions_DF <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_intermeasure)
DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative)
DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT_RATIO <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative_RATIO)
DATA_OBJECTS.EllipticalRegressions_DF_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_elliptical)

DATA_OBJECTS.ComparativeRegressions_REIS_AREA_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_REIS_AREA)

DATA_OBJECTS.ComparativeRegressions_REIS_DRYWEIGHT_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_REIS_DRYWEIGHT)

# DATA_OBJECT Generation: Filtered ---------------------------------------------
cat("Generating Filtered Datasets on 3p/n Threshold Value by Model\n")

cat("::| Writing Blacklists...\n")
DATA_OBJECTS.blacklist_3PN <- UTIL.Assemble3PNDroplist(DATA_OBJECTS.LinearRegressions_DF[,-which(names(DATA_OBJECTS.LinearRegressions_DF) %in% c("thickness~width", "thickness~height", "thickness~diameter"))], SWITCHBOARD.csvMAINFILE)


cat("::| Writing Filtered Datasets...\n")
DATA_OBJECTS.csvAUGM_3PN <- UTIL.DropOnNPN(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.blacklist_3PN)

DATA_OBJECTS.ComparativeRegressions_REIS_AREA_3PN <- CALC.ModelGenerator(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_REIS_AREA)

DATA_OBJECTS.ComparativeRegressions_REIS_DRYWEIGHT_3PN <- CALC.ModelGenerator(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_REIS_DRYWEIGHT)

# Summary Statistic Heatmap Generation ------------------------------------------------------

cat("Generating Heatmaps and Writing CSV Data\n")
cat("::| Heatmaps...\n")
SWITCHBOARD.models_multiplicative.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative, 5)
SWITCHBOARD.models_multiplicative.splitinteraction_RATIO <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_RATIO, 5)
SWITCHBOARD.models_elliptical.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical, 2)

SWITCHBOARD.models_multiplicative_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_abbr, 5)
SWITCHBOARD.models_multiplicative_abbr.splitinteraction_RATIO <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_abbr_RATIO, 5)
SWITCHBOARD.models_elliptical_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical_abbr, 2)

cat("::| CSVs...\n")
DATA_OBJECTS.comparison_mult_R2_FILT <- DATA_OBJECTS.csvAUGM_3PN %>%
  CALC.ModelGenerator(DATA_OBJECTS.MODELS_multiplicative) %>%
  CALC.df_R2() %>%
  UTIL.AllAccessionToMain(CALC.ModelGenerator_AllAccessions_R2(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_multiplicative))

DATA_OBJECTS.comparison_mult_R2_FILT_RATIO <- DATA_OBJECTS.csvAUGM_3PN %>%
  CALC.ModelGenerator(DATA_OBJECTS.MODELS_multiplicative_RATIO) %>%
  CALC.df_R2() %>%
  UTIL.AllAccessionToMain(CALC.ModelGenerator_AllAccessions_R2(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_multiplicative_RATIO))


DATA_OBJECTS.comparison_ellip_R2_FILT <- DATA_OBJECTS.csvAUGM_3PN %>%
  CALC.ModelGenerator(DATA_OBJECTS.MODELS_elliptical) %>%
  CALC.df_R2() %>%
  UTIL.AllAccessionToMain(CALC.ModelGenerator_AllAccessions_R2(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_elliptical))



DATA_OBJECTS.Comparison.REIS.area <- as.data.frame(SWITCHBOARD.models_to_compare.REIS.AREA) %>%
  set_colnames(c("models"))

DATA_OBJECTS.Comparison.REIS.dryweight <- as.data.frame(SWITCHBOARD.models_to_compare.REIS.DRYWEIGHT) %>%
  set_colnames(c("models"))
