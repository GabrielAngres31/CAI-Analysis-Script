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
setwd("C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script")

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
  #"I(width*height*diameter*thickness)"
)

SWITCHBOARD.models_multiplicative_abbr <- c(
  "FW~W", "FW~H", "FW~D", "FW~T",
  "FW~D/W",
  "FW~W+H", "FW~W*H",
  "FW~W+D", "FW~W*D",
  "FW~W+T", "FW~W*T",
  "FW~H+D", "FW~H*D",
  "FW~H+T", "FW~H*T",
  "FW~D+T", "FW~D*T",
  "FW~W+H+D", "FW~W*H*D",
  "FW~W+H+T", "FW~W*H*T",
  "FW~W+D+T", "FW~W*D*T",
  "FW~H+D+T", "FW~H*D*T",
  "FW~W+H+D+T", "FW~W*H*D*T"
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
  "FW~A", "FW~Th_A",
  "FW~A+T", "FW~A*T",
  "FW~Th_A+T", "FW~Th_A*T",
  "FW~Th_A+T+PDD", "FW~Th_A*T*PDD"
)


# DATA_OBJECT Generation: Unfiltered --------------------------------------------------

cat("::| Polishing Model Guides...\n")
DATA_OBJECTS.MODELS_intermeasure$models <- apply(DATA_OBJECTS.MODELS_intermeasure[, c("yname", "xname")], 1, paste, collapse = "~") 

DATA_OBJECTS.MODELS_multiplicative <- as.data.frame(SWITCHBOARD.models_multiplicative) %>%
  set_colnames(c("models"))

DATA_OBJECTS.MODELS_elliptical <- as.data.frame(SWITCHBOARD.models_elliptical) %>%
  set_colnames(c("models"))

cat("Generating Model Dataframes\n")
DATA_OBJECTS.LinearRegressions_DF <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_intermeasure)
DATA_OBJECTS.MultiplicativeRegressions_DF_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative)
DATA_OBJECTS.EllipticalRegressions_DF_UNFILT <- CALC.ModelGenerator(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_elliptical)


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

cat("Generating Heatmaps and Writing CSV Data\n")
UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_intermeasure, "Lin_UNF",   "Heatmaps", SWITCHBOARD.models_intermeasure_abbr)

UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_multiplicative, "Mult_UNF",  "Heatmaps", SWITCHBOARD.models_multiplicative_abbr)
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_multiplicative, "Mult_3PN",  "Heatmaps", SWITCHBOARD.models_multiplicative_abbr)

UTIL.StoreHeatmapsAndFile(SWITCHBOARD.csvAUGMFILE, DATA_OBJECTS.MODELS_elliptical, "Ellip_UNF", "Heatmaps", SWITCHBOARD.models_elliptical_abbr)
UTIL.StoreHeatmapsAndFile(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.MODELS_elliptical, "Ellip_3PN", "Heatmaps", SWITCHBOARD.models_elliptical_abbr)