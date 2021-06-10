#Establish working directory with target files.
setwd("C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script")

#Load relevant packages into library
library("xlsx")
library("dplyr")
library("tibble")
library("rsq")
library("qwraps2")
library("ggplot2")
library("car")
library("dunn.test")

#SWITCHBOARD
#   Hacky way of getting consistent legend generation, also controls value rounding over the program
#     and some contant definitions

SWITCHBOARD.CLEAN_DATA <- FALSE
SWITCHBOARD.strALLACCESSIONS <- "All Accessions"
SWITCHBOARD.strALLDATA <- "Entire"
SWITCHBOARD.roundto <- 3
SWITCHBOARD.STDDEV <- 3
SWITCHBOARD.strACCESSIONLIST <-
  c(
    "242",
    "246",
    "319",
    "325",
    "326",
    "390",
    "572",
    "580",
    "582",
    "584",
    "585",
    "839",
    "845",
    "854"
  )
SWITCHBOARD.strQQMEASURESLIST <-
  c(
    "fresh_weight", 
    "width",
    "height",
    "diameter",
    "thickness",
    "D_div_W",
    "FW_div_H",
    "FW_div_W",
    "FW_div_D",
    "FW_div_T"
  )

#List of all plots to be generated.
SWITCHBOARD.GRAPHS_LIST <- tribble(
  ~ xname, ~ yname,
  "width",  "fresh_weight",
  'height', 'fresh_weight',
  'diameter', 'fresh_weight',
  'thickness', 'fresh_weight',
  'D_div_W', 'fresh_weight',
  'width', 'height',
  'width', 'diameter',
  'width', 'thickness',
  'height', 'diameter',
  'height', 'thickness',
  'diameter', 'thickness',
  #  'add_HW', 'thickness',
  #  'add_HW', 'fresh_weight',
  #  'add_HWT', 'fresh_weight',
  #  'add_DT', 'fresh_weight',
  #  'mul_HWDT', 'fresh_weight'
)

#R^2 Output csv dataframes
csvR2Frame <- data.frame(Graph = character(), R2 = double())

#List of anomalous datapads--TODO: Establish statistical, heuristic basis for removal
errorPads <-
  vector() #c("PARL 585-1-8", "PARL 319-1-1", "PARL 854-2-3", "PARL 580-1-2", "PARL 584-1-1", "PARL 242-4-3")

#FUNCTION DEFINITIONS

'
speciesIdent <- function(ID) {
  species_pass_value = switch(
    substr(ID, 6, 8),
    "580" = "O. ficus-indica",
    "854" = "O. ficus-indica",
    "584" = "O. ficus-indica",
    "585" = "O. ficus-indica",
    "572" = "O. ficus-indica",
    "325" = "O. crassa",
    "845" = "O. ficus-indica",
    "246" = "O. ficus-indica",
    "326" = "O. ficus-indica",
    "390" = "O. robusta",
    "242" = "N. cochinillifera",
    "582" = "N. cochinillifera",
    "319" = "N. cochinillifera",
    "839" = "O. undulata"
  )
  return(species_pass_value)
}
'

#For AIC, etc. parameter calculation
'
paraNum <- function(colname) {
  if (colname == "add_HW" | colname == "add_DT") {
    return(2)
  }
  else if (colname == ("add_HWT")) {
    return(3)
  }
  else if (colname == ("mul_HWDT")) {
    return(4)
  }
  else {
    return(1)
  }
}
'
#Just for easy parameterization in the Generator function
ACCESSORY.colnameToLegend <- function(column_text) {
  legname <- switch(
    column_text,
    "width" = "Width",
    "height" = "Height",
    "diameter" = "Diameter",
    "thickness" = "Thickness",
    "fresh_weight" = "Fresh Weight",
    'D_div_W' = "Diameter / Width",
    'FW_div_H'='Fresh Weight / Height',
    'FW_div_W'='Fresh Weight / Width',
    'FW_div_D'='Fresh Weight / Diameter',
    'FW_div_T'='Fresh Weight / T'
  )
  #    "add_HW" = "Height + Width",
  #    "add_HWT" = "Height + Width + Thickness",
  #    "add_DT" = "Diameter + Thickness",
  #    "mul_HWDT" = "Height * Weight * Diameter * Thickness",
  #    "fresh_weight" = "Fresh Weight"
  #  )
  return(legname)
}

#Functions for error removal
#tagBySTDDEV--Called by ErrorCleaning before (in?) removeByID, this tags observations
# if they are more than n = 3 standard deviations away from the mean in a given parameter.
# However, may have to remove as fresh weight is not normally distributed.

ACCESSORY.tagBySTDDEV <-
  function(dataset, an_accession) {
    blacksite <- vector()
    #print(blacksite)
    dataset_accessionfilter <-
      filter(dataset, accession == an_accession)
    #for (columnname in c("H_div_W", "FW_div_W", "FW_div_D")) {
    for (columnname in c("D_div_W", "H_div_W")) {
      tempcol <-
        dataset_accessionfilter[[columnname]]
      meanCol <- mean(tempcol)
      STDCol <- sd(tempcol)
      bound_upper <- (meanCol + (STDCol * SWITCHBOARD.STDDEV))
      bound_lower <- (meanCol - (STDCol * SWITCHBOARD.STDDEV))
      detention <-
        filter(
          dataset_accessionfilter, !!as.symbol(columnname) < bound_lower |
            !!as.symbol(columnname) > bound_upper
        )
      blacklist <- as.vector(detention$completeID)
      blacksite <- c(blacksite, blacklist)
      print(c("List:", blacklist))
      print(c("Site:", blacksite))
      print(c("Column:", columnname))
    }
    return(blacksite)
  }

#ErrorCleaning--
ACCESSORY.ErrorCleaning <- function(dataset, doClean) {
  rawblock <- dataset
  #print(nrow(rawblock))
  if (isTRUE(doClean)) {
    for (accession in SWITCHBOARD.strACCESSIONLIST) {
      blacksite <-
        ACCESSORY.tagBySTDDEV(dataset, accession)
      rawblock <-
        filter(rawblock, !(rawblock$completeID %in% blacksite))
      #print(nrow(rawblock))
    }
    finalDataset <- rawblock
    return(finalDataset)
  } else {
    return(dataset)
  }
}

#General Data Subset Function

ACCESSORY.DataSubset <-
  function(filter_data,
           filter_accession,
           dataset_in) {
    if (filter_data == SWITCHBOARD.strALLDATA) {
      F_data <- dataset_in
    } else {
      print(dataset_in)
      F_data <- filter(dataset_in, dataset == filter_data)
    }
    if (filter_accession == SWITCHBOARD.strALLACCESSIONS) {
      F_data_accession <- F_data
    } else {
      F_data_accession <-
        filter(F_data, accession == filter_accession)
    }
    if (nrow(F_data_accession) == 0) {
      return(NULL)
    } else {
      return(F_data_accession)
    }
  }

#Report Generation
###Plotter function which creates labeled figure, regressions

PIPELINE.PlotRegress = function(x,
                       y,
                       year = SWITCHBOARD.strALLDATA,
                       accession = SWITCHBOARD.strALLACCESSIONS,
                       dataset_in) {
  subset <-
    ACCESSORY.DataSubset(year, accession, dataset_in)
  if (is.null(subset) == TRUE) {
    return(NULL)
  }
  xlab = ACCESSORY.colnameToLegend(x)
  ylab = ACCESSORY.colnameToLegend(y)
  subslice = select(subset, x, y)
  x_data = subslice[[x]] #colnameToCol(x, subslice)
  y_data = subslice[[y]] #colnameToCol(y, subslice)
  
  #Double logistic plot currently the only one implemented.
  doublelogline = lm(log(y_data) ~ log(x_data), data = subslice)
  #print(c(species, xlab, " vs. ", ylab))
  #print(AIC(plotline, k = (2+paraNum(x))))
  #print("-----------------------")
  rsquare = rsq(doublelogline, type = "sse")
  intercept <- coefficients(doublelogline)[1]
  slope <- coefficients(doublelogline)[2]
  details <-
    paste0("y_0 = ",
           round(intercept, SWITCHBOARD.roundto),
           ", m = ",
           round(slope, SWITCHBOARD.roundto),
           collapse = " ")
  
  csvR2Frame <<-
    rbind(csvR2Frame,
          data.frame(
            "Dataset" = paste0(year, " Dataset ", accession, collapse = " "),
            "Graph" = paste0(xlab, " vs. ", ylab,
                             collapse = " "),
            "R2" = round(rsquare, SWITCHBOARD.roundto)
          ))
  plot(
    x_data,
    y_data,
    type = "p",
    log = "xy",
    cex = 0.8,
    sub = paste0(
      details,
      ", R^2 =",
      round(rsquare, SWITCHBOARD.roundto),
      ", N = ",
      nrow(subslice)
    ),
    main = paste0(xlab, " vs. ", ylab, "\n", year, " Dataset, ", accession, collapse = " "),
    xlab = xlab,
    ylab = ylab
  )
  lines(x_data, exp(predict(
    doublelogline, newdata = list(x_data = x_data)
  )) , col = "grey")
  #print(cat(paste(species, year, x, " vs. ", y, round(rsquare, SWITCHBOARD.roundto), collapse = "\t", '\n')))
}

###Generator function which generates plots using Plotter over all of the axes combinations

PIPELINE.CorrPlotsGenerator <-
  function(figure_guide,
           year = SWITCHBOARD.strALLDATA,
           accession = SWITCHBOARD.strALLACCESSIONS,
           dataset_in) {
    for (i in 1:nrow(figure_guide)) {
      row <- figure_guide[i,]
      PIPELINE.PlotRegress(row$xname, row$yname, year, accession, dataset_in)
    }
  }

ACCESSORY.qqGen <- function(accession, column, dataset_in, doClean) {
  subsetA <-
    ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA, filter_accession = accession, dataset_in = dataset_in)
  title <-
    paste0(accession, ", ", column)
  png(
    paste0(
      "C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\QQPlot_Images",
      ifelse(doClean, "_Corrected", ""),
      "\\",
      column,
      "--",
      accession,
      ifelse(doClean, "_Corrected", ""),
      ".png"
    )
  )
  #qqnorm(subsetA[[column]], pch = 1, frame = TRUE, main = title)
  #qqline(subsetA[[column]], col = "steelblue", lwd = 2)
  qqPlot(subsetA[[column]], main = title, envelope = 0.95)
  #ggsave(paste0("C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\", accession, "--", column))
  dev.off()
}

ACCESSORY.qqGenLog <- function(accession, column, dataset_in, doClean) {
  subsetA <-
    ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA, filter_accession = accession, dataset_in = dataset_in)
  title <-
    paste0(accession, ", ", column)
  png(
    paste0(
      "C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\QQPlot_Images",
      ifelse(doClean, "_Corrected", ""),
      "\\",
      column,
      "--",
      accession,
      "_LOG",
      ifelse(doClean, "_Corrected", ""),
      ".png"
    )
  )
  #qqnorm(subsetA[[column]], pch = 1, frame = TRUE, main = title)
  #qqline(subsetA[[column]], col = "steelblue", lwd = 2)
  qqPlot(log(subsetA[[column]]), main = title, envelope = 0.95)
  #ggsave(paste0("C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\", accession, "--", column))
  dev.off()
}

ACCESSORY.qqGenSqrt <- function(accession, column, dataset_in = dataset_in, doClean) {
  subsetA <-
    ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA, filter_accession = accession, dataset_in)
  title <-
    paste0(accession, ", ", column)
  png(
    paste0(
      "C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\QQPlot_Images",
      ifelse(doClean, "_Corrected", ""),
      "\\",
      column,
      "--",
      accession,
      "_SQRT", 
      ifelse(doClean, "_Corrected", ""),
      ".png"
    )
  )
  #qqnorm(subsetA[[column]], pch = 1, frame = TRUE, main = title)
  #qqline(subsetA[[column]], col = "steelblue", lwd = 2)
  qqPlot(sqrt(subsetA[[column]]), main = title, envelope = 0.95)
  #ggsave(paste0("C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\", accession, "--", column))
  dev.off()
}

PIPELINE.PlotCompiler <- function(dataset_in, doClean) {
  for (dataset in c(SWITCHBOARD.strALLDATA)) {
    #..., "Year 1", "Year 2"
    for (accession in c(
      SWITCHBOARD.strALLACCESSIONS,
      "242",
      "246",
      "319",
      "325",
      "326",
      "390",
      "572",
      "580",
      "582",
      "584",
      "585",
      "839",
      "845",
      "854"
    )) {
      PIPELINE.CorrPlotsGenerator(SWITCHBOARD.GRAPHS_LIST, dataset, accession, dataset_in)
      
      #print(c("Generating Plots for ", i, j))
    }
  }
  
  for (accession in SWITCHBOARD.strACCESSIONLIST) {
    for (measure in SWITCHBOARD.strQQMEASURESLIST) {
      ACCESSORY.qqGen(accession, measure, dataset_in, doClean)
      #ACCESSORY.qqGenLog(accession, measure, dataset_in, doClean)
      #ACCESSORY.qqGenSqrt(accession, measure, dataset_in, doClean)
    }
  }
}

main <- function() {
  
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  #-----------------------------------------------
  #_______________________________________________
  #DATASET MANIPULATION, FUNCTION EXECUTION SEGMENT
  
  #Load Parlier csv files into memory
  PARLIER <- read.csv(file = 'PARL0_04262021.csv')
  
  #Data compatibility modifications
  PARLIER <- mutate(PARLIER, pad_id = 1)
  PARLIER$dataset <- "Year 1"
  
  #Add Error-correction derived measures
  
  PARLIER <- mutate(PARLIER, H_div_W = height / width)
  PARLIER <- mutate(PARLIER, FW_div_W = fresh_weight / width)
  PARLIER <- mutate(PARLIER, FW_div_D = fresh_weight / diameter)
  PARLIER <- mutate(PARLIER, FW_div_H = fresh_weight / height)
  PARLIER <- mutate(PARLIER, FW_div_T = fresh_weight / thickness)
  PARLIER <- mutate(PARLIER, D_div_W = diameter / width)

  for(doClean in c(FALSE, TRUE)) {
    
    finishedPARLIER <- ACCESSORY.ErrorCleaning(PARLIER, doClean)
    write.csv(
      finishedPARLIER,
      paste0(
        "C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\SAS_Datasets\\PARL0", 
        ifelse(doClean, "C", ""),
        ".csv"
      ),
      row.names = FALSE
      )
    write.csv(
      finishedPARLIER,
      paste0(
        "C:\\Users\\gjang\\Documents\\SASUniversityEdition\\myfolders\\PARL0", 
        ifelse(doClean, "C", ""),
        ".csv"
      ),
      row.names = FALSE
    )
  
    title <- paste(
      "Double-Log Cladode Parameter Cross-Relations",
      ifelse(doClean, " Corrected", ""),
      ".pdf",
      sep = ""
    )
    pdf(title)
    PIPELINE.PlotCompiler(finishedPARLIER, doClean)
    dev.off()
    write.csv(
      csvR2Frame,
      paste(
        "C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\R^2_Records\\Double-Log_Regression_R^2_Values",
        ifelse(doClean, "_Corrected", ""),
        ".csv",
        sep = ""
      ),
      row.names = FALSE
  )
  }
}



#Run Error-cleaning Script

#plot(curve(dweibull(x, shape=2, scale = 1), from=0, to=max(finishedPARLIER$fresh_weight)))

#Functions for creating new columns

###Create Function Measures

main()
