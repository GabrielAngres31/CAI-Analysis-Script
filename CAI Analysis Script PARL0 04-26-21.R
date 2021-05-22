#Establish working directory with target files.
setwd("C:\\Users\\gjang\\Desktop\\Cushman Materials\\")

#Load relevant packages into library
library("xlsx")
library("dplyr")
library("tibble")
library("rsq")
library("qwraps2")
library("ggplot2")
library("car")
library("dunn.test")

#SWITCHBOARD--Hacky way of getting consistent legend generation, also controls value rounding.
CLEAN_DATA <- TRUE
strALLACCESSIONS <-
  "All Accessions" #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
strALLDATA <-
  "Entire"#<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
roundto <- 3
STDDEV <- 3
strACCESSIONLIST <-
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

#List of all plots to be generated.
GRAPHS_LIST <- tribble(
  ~ xname,  ~ yname,
  "width",  "fresh_weight",
  'height',  'fresh_weight',
  'diameter',  'fresh_weight',
  'thickness',  'fresh_weight',
  'DdivW', 'fresh_weight',
  'width',  'height',
  'width',  'diameter',
  'width',  'thickness',
  'height',  'diameter',
  'height',  'thickness',
  'diameter',  'thickness',
  
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

#Just for easy parameterization in the Generator function
colnameToLegend <- function(column_text) {
  legname <- switch(
    column_text,
    "width" = "Width",
    "height" = "Height",
    "diameter" = "Diameter",
    "thickness" = "Thickness",
    "fresh_weight" = "Fresh Weight",
    'DdivW' = "Diameter / Width"
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

tagBySTDDEV <-
  function(dataset, an_accession) {
    #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
    blacksite <- vector()
    #print(blacksite)
    dataset_accessionfilter <-
      filter(dataset, accession == an_accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
    #for (columnname in c("H_div_W", "FW_div_W", "FW_div_D")) {
    for (columnname in c("DdivW", "H_div_W")) {
      tempcol <-
        dataset_accessionfilter[[columnname]] #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
      meanCol <- mean(tempcol)
      STDCol <- sd(tempcol)
      bound_upper <- (meanCol + (STDCol * STDDEV))
      bound_lower <- (meanCol - (STDCol * STDDEV))
      detention <-
        filter(
          dataset_accessionfilter,
          !!as.symbol(columnname) < bound_lower |
            !!as.symbol(columnname) > bound_upper
        ) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
      blacklist <- as.vector(detention$completeID)
      blacksite <- c(blacksite, blacklist)
      print(c("List:", blacklist))
      print(c("Site:", blacksite))
      print(c("Column:", columnname))
    }
    return(blacksite)
  }

#ErrorCleaning--
ErrorCleaning <- function(dataset) {
  rawblock <- dataset
  #print(nrow(rawblock))
  if (isTRUE(CLEAN_DATA)) {
    for (accession in strACCESSIONLIST) {
      #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
      blacksite <-
        tagBySTDDEV(PARLIER, accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
      rawblock <-
        filter(rawblock,!(rawblock$completeID %in% blacksite))
      #print(nrow(rawblock))
    }
    finalPARLIER <- rawblock
    return(finalPARLIER)
  } else {
    return(PARLIER)
  }
}

#General Data Subset Function

DataSubset <-
  function(filter_data = strALLDATA,
           filter_accession = strALLACCESSIONS) {
    #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
    if (filter_data == strALLDATA) {
      F_data <- finishedPARLIER
    } else {
      F_data <- filter(finishedPARLIER, dataset == filter_data)
    }
    if (filter_accession == strALLACCESSIONS) {
      #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
      F_data_accession <- F_data
    } else {
      F_data_accession <-
        filter(F_data, accession == filter_accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
    }
    if (nrow(F_data_accession) == 0) {
      #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
      return(NULL)
    } else {
      return(F_data_accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
    }
  }

#Report Generation
###Plotter function which creates labeled figure, regressions

PlotRegress = function(x,
                       y,
                       year = strALLDATA,
                       accession = strALLACCESSIONS) {
  #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  subset <-
    DataSubset(year, accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  if (is.null(subset) == TRUE) {
    return(NULL)
  }
  xlab = colnameToLegend(x)
  ylab = colnameToLegend(y)
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
           round(intercept, roundto),
           ", m = ",
           round(slope, roundto),
           collapse = " ")
  
  csvR2Frame <<-
    rbind(csvR2Frame,
          data.frame(
            "Dataset" = paste0(year, " Dataset ", accession, collapse = " "),
            "Graph" = paste0(xlab, " vs. ", ylab,
                             collapse = " "),
            "R2" = round(rsquare, roundto)
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
      round(rsquare, roundto),
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
  #print(cat(paste(species, year, x, " vs. ", y, round(rsquare, roundto), collapse = "\t", '\n')))
}

###Generator function which generates plots using Plotter over all of the axes combinations

CorrPlotsGenerator <-
  function(figure_guide,
           year = strALLDATA,
           accession = strALLACCESSIONS) {
    #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
    for (i in 1:nrow(figure_guide)) {
      row <- figure_guide[i, ]
      PlotRegress(row$xname, row$yname, year, accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
    }
  }

qqGen <- function(accession, column) {
  subsetA <-
    DataSubset(filter_accession = accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  title <-
    paste0(accession, ", ", column) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  png(
    paste0(
      "C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\",
      column,
      "--",
      accession,
      ".png"
    )
  ) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  #qqnorm(subsetA[[column]], pch = 1, frame = TRUE, main = title)
  #qqline(subsetA[[column]], col = "steelblue", lwd = 2)
  qqPlot(subsetA[[column]], main = title, envelope = 0.95)
  #ggsave(paste0("C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\", accession, "--", column))
  dev.off()
}

qqGenLog <- function(accession, column) {
  subsetA <-
    DataSubset(filter_accession = accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  title <-
    paste0(accession, ", ", column) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  png(
    paste0(
      "C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\",
      column,
      "--",
      accession,
      "_LOG.png"
    )
  ) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  #qqnorm(subsetA[[column]], pch = 1, frame = TRUE, main = title)
  #qqline(subsetA[[column]], col = "steelblue", lwd = 2)
  qqPlot(log(subsetA[[column]]), main = title, envelope = 0.95)
  #ggsave(paste0("C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\", accession, "--", column))
  dev.off()
}

qqGenSqrt <- function(accession, column) {
  subsetA <-
    DataSubset(filter_accession = accession) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  title <-
    paste0(accession, ", ", column) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  png(
    paste0(
      "C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\",
      column,
      "--",
      accession,
      "_SQRT.png"
    )
  ) #<@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%#@#%>
  #qqnorm(subsetA[[column]], pch = 1, frame = TRUE, main = title)
  #qqline(subsetA[[column]], col = "steelblue", lwd = 2)
  qqPlot(sqrt(subsetA[[column]]), main = title, envelope = 0.95)
  #ggsave(paste0("C:\\Users\\gjang\\Desktop\\Cushman Materials\\QQPlot images\\", accession, "--", column))
  dev.off()
}

PlotCompiler <- function() {
  for (i in c(strALLDATA)) {
    #..., "Year 1", "Year 2"
    for (j in c(
      strALLACCESSIONS,
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
      CorrPlotsGenerator(GRAPHS_LIST, i, j)
      
      #print(c("Generating Plots for ", i, j))
    }
  }
  for (k in strACCESSIONLIST) {
    for (l in c("fresh_weight", "width", "height", "diameter", "thickness")) {
      qqGen(k, l)
      qqGenLog(k, l)
      qqGenSqrt(k, l)
    }
  }
}

main <- function() {
  title <- paste(
    "Double-Log Cladode Parameter Cross-Relations",
    ifelse(CLEAN_DATA, " Corrected", ""),
    ".pdf",
    sep = ""
  )
  pdf(title)
  PlotCompiler()
  dev.off()
  write.csv(
    csvR2Frame,
    paste(
      "C:\\Users\\gjang\\Desktop\\Cushman Materials\\R^2 Records",
      ifelse(CLEAN_DATA, " Corrected", ""),
      ".csv",
      sep = ""
    ),
    row.names = FALSE
  )
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#-----------------------------------------------
#_______________________________________________
#DATASET MANIPULATION, FUNCTION EXECUTION SEGMENT

#Load Parlier csv files into memory
parlierOld <- read.csv(file = 'PARL0_04262021.csv')

#Data compatibility modifications
parlierOld <- mutate(parlierOld, pad_id = 1)
parlierOld$dataset <- "Year 1"

#Merge Y1 and Y2 Parlier Datasets
PARLIER <- rbind(parlierOld) #, parlierNew)

#Add Error-correction ratios

PARLIER <- mutate(PARLIER, H_div_W = height / width)
PARLIER <- mutate(PARLIER, FW_div_W = fresh_weight / width)
PARLIER <- mutate(PARLIER, FW_div_D = fresh_weight / diameter)
PARLIER <- mutate(PARLIER, DdivW = diameter / width)

#Run Error-cleaning Script
finishedPARLIER <- ErrorCleaning(PARLIER)
#plot(curve(dweibull(x, shape=2, scale = 1), from=0, to=max(finishedPARLIER$fresh_weight)))

#Functions for creating new columns

###Create Function Measures
finishedPARLIER <-
  mutate(finishedPARLIER, add_HW = (height + width))

finishedPARLIER <-
  mutate(finishedPARLIER, add_HWT = height + width + thickness)

finishedPARLIER <-
  mutate(finishedPARLIER, add_DT = diameter + thickness)

finishedPARLIER <-
  mutate(finishedPARLIER, mul_HWDT = height * width * diameter * thickness)




main()
