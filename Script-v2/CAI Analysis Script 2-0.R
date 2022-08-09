# This program was written to automate on a large scale the statistical analyses
# done in R, in order to save time - so that minor changes to analytical processes
# may be quickly reflected across the entire analysis. This allows hours of work
# to be repeated, with or without changes, in less than two minutes while analyzing 
# thousands of datapoints.
#
# The program is interspersed with diagnostic output messages ("cat(...)") that 
# output to the console and indicate the progression of the program.
# 
# Contact gangres@nevada.unr.edu for questions regarding the code's function and 
# behavior.

#Establish working directory with target files.


# Workstation Setup -------------------------------------------------------

cat("Setting Working Directory...\n")
setwd("C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\Script-v2")

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
library(RColorBrewer)
library("data.table")



# SWITCHBOARD Setup -------------------------------------------------------

cat("Setting up SWITCHBOARD: \n")

SWITCHBOARD.DIRECTORY <-
  "C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\Script-v2"

SWITCHBOARD.subdirectories <- c("INTERMEASURE", 
                                "QQPlot_Images", 
                                "Tukey_Graphs",
                                "Tukey_Graphs\\Tukey_Hypergraphs",
                                "Tukey_Graphs\\Tukey_Intercomparisons",
                                "Heatmaps",
                                "Heatmaps\\Graphs", 
                                "Heatmaps\\Tables")

for (subdir in SWITCHBOARD.subdirectories) {
  dirpath <- file.path(SWITCHBOARD.DIRECTORY, subdir)
    if (!(dir.exists(dirpath))) {
      dir.create(dirpath)
    }
}

cat("::| Loading Dataset of Interest...\n")
SWITCHBOARD.csvMAINFILE <- read.csv(file = 'PARL0_06202022.csv')

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
  mutate(Pade_Derived_Diam = Pade_Peri/pi)

cat("::| Setting Accession List...\n")
SWITCHBOARD.ACCESSIONLIST <- c(242, 246, 319, 325, 326, 390, 572, 580, 582, 584, 585, 839, 845, 854)

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

SWITCHBOARD.models_intermeasure_abbr <- c("FW~W", "FW~H", "FW~D", "FW~T", "H~W", "D~W", "T~W", "D~H", "T~H", "T~D")

SWITCHBOARD.qqMeasures_VEC <- c("width", "diameter", "height", "thickness", "fresh_weight")

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


# DATA_OBJECT Generation: Unfiltered --------------------------------------------------

cat("::| Polishing Model Guides...\n")
DATA_OBJECTS.MODELS_intermeasure$models <- apply(DATA_OBJECTS.MODELS_intermeasure[, c("yname", "xname")], 1, paste, collapse = "~") 

DATA_OBJECTS.MODELS_multiplicative <- as.data.frame(SWITCHBOARD.models_multiplicative) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_elliptical <- as.data.frame(SWITCHBOARD.models_elliptical) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_multiplicative_abbr <- as.data.frame(SWITCHBOARD.models_multiplicative_abbr) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_elliptical_abbr <- as.data.frame(SWITCHBOARD.models_elliptical_abbr) %>%
  set_colnames(c("models"))

DATA_OBJECTS.models_multiplicative.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative, 4)
DATA_OBJECTS.models_elliptical.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical, 2)

DATA_OBJECTS.models_multiplicative_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_abbr, 4)
DATA_OBJECTS.models_elliptical_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical_abbr, 2)

cat("Generating Model Dataframes\n")
DATA_OBJECTS.LinearRegressions_DF <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_intermeasure)
DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative)
DATA_OBJECTS.EllipticalRegressions_DF_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_elliptical)



# DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT.split_plus <- DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT[SWITCHBOARD.models_multiplicative.splitinteraction[[1]]]
# DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT.split_mult <- DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT[SWITCHBOARD.models_multiplicative.splitinteraction[[2]]]
# 
# DATA_OBJECTS.EllipticalRegressions_DF_UNFILT.split_plus <- DATA_OBJECTS.EllipticalRegressions_DF_UNFILT[SWITCHBOARD.models_elliptical.splitinteraction[[1]]]
# DATA_OBJECTS.EllipticalRegressions_DF_UNFILT.split_mult <- DATA_OBJECTS.EllipticalRegressions_DF_UNFILT[SWITCHBOARD.models_elliptical.splitinteraction[[2]]]



# DATA_OBJECT Generation: Filtered ---------------------------------------------
cat("Generating Filtered Datasets on 3p/n Threshold Value by Model\n")

cat("::| Writing Blacklists...\n")
DATA_OBJECTS.blacklist_3PN <- UTIL.Assemble3PNDroplist(DATA_OBJECTS.LinearRegressions_DF[,-which(names(DATA_OBJECTS.LinearRegressions_DF) %in% c("thickness~width", "thickness~height", "thickness~diameter"))], SWITCHBOARD.csvMAINFILE)


cat("::| Writing Filtered Datasets...\n")
DATA_OBJECTS.csvAUGM_3PN <- UTIL.DropOnNPN(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.blacklist_3PN)


# Quantile-Quantile Plot Generation ---------------------------------------

cat("Generating Quantile-Quantile Plots\n")
COMPILE.qqGen(SWITCHBOARD.csvMAINFILE)
# Intermeasure and Tukeygraph Generation ---------------------------------

cat("Generating Intermeasure Plots\n")
COMPILE.intermeasure(DATA_OBJECTS.LinearRegressions_DF)

cat("Generating Tukey Intercomparisons and Grouping Hypergraphs\n")
COMPILE.tukeygraphs(SWITCHBOARD.csvAUGMFILE)

# Summary Statistic Heatmap Generation ------------------------------------------------------

DATA_OBJECTS.csvAUGM_3PN
DATA_OBJECTS.EllipticalRegressions_DF_UNFILT

cat("Generating Heatmaps and Writing CSV Data\n")

SWITCHBOARD.models_multiplicative.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative, 5)
SWITCHBOARD.models_elliptical.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical, 2)

SWITCHBOARD.models_multiplicative_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_multiplicative_abbr, 5)
SWITCHBOARD.models_elliptical_abbr.splitinteraction <- UTIL.SplitModelListOnInteractions(DATA_OBJECTS.MODELS_elliptical_abbr, 2)

UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_intermeasure, "Lin_UNF",   "Heatmaps", "Intermeasure Models", SWITCHBOARD.models_intermeasure_abbr)

UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative, "Mult_UNF",  "Heatmaps", "Multiplicative Models - Unfiltered Dataset",SWITCHBOARD.models_multiplicative_abbr)
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_multiplicative, "Mult_3PN",  "Heatmaps", "Multiplicative Models - 3P/N Dataset", SWITCHBOARD.models_multiplicative_abbr)

UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_elliptical, "Ellip_UNF", "Heatmaps", "Elliptical Models - Unfiltered Dataset", SWITCHBOARD.models_elliptical_abbr)
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_elliptical, "Ellip_3PN", "Heatmaps", "Elliptical Models - 3P/N Dataset", SWITCHBOARD.models_elliptical_abbr)

#---

UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_multiplicative.splitinteraction[[1]], "Mult_UNF_plus",  "Heatmaps", "Multiplicative  Models - Unfiltered Dataset with No Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[1]]$models)
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_multiplicative.splitinteraction[[2]], "Mult_UNF_mult",  "Heatmaps", "Multiplicative  Models- Unfiltered Dataset with Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[2]]$models)

UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_multiplicative.splitinteraction[[1]], "Mult_3PN_plus",  "Heatmaps", "Multiplicative Models - 3P/N Dataset with No Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[1]]$models)
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_multiplicative.splitinteraction[[2]], "Mult_3PN_mult",  "Heatmaps", "Multiplicative Models - 3P/N Dataset with Interactions", DATA_OBJECTS.models_multiplicative_abbr.splitinteraction[[2]]$models)

UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_elliptical.splitinteraction[[1]], "Ellip_UNF_plus",  "Heatmaps", "Elliptical Models - Unfiltered Dataset with No Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[1]]$models)
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.models_elliptical.splitinteraction[[2]], "Ellip_UNF_mult",  "Heatmaps", "Elliptical Models- Unfiltered Dataset with Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[2]]$models)

UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_elliptical.splitinteraction[[1]], "Ellip_3PN_plus",  "Heatmaps", "Elliptical Models - 3P/N Dataset with No Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[1]]$models)
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.models_elliptical.splitinteraction[[2]], "Ellip_3PN_mult",  "Heatmaps", "Elliptical Models - 3P/N Dataset with Interactions", DATA_OBJECTS.models_elliptical_abbr.splitinteraction[[2]]$models)

