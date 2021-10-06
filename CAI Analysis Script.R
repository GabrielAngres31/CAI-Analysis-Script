#Box Version 1.0

#Establish working directory with target files.
setwd("C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script")

#Load relevant packages into library
#library("xlsx")
library("dplyr")
library("tibble")
library("rsq")
library("car")
library("Rcmdr")
library("hash")

SWITCHBOARD.csvMAINFILE <- read.csv(file = 'PARL0_09092021.csv')

#SWITCHBOARD
#   Hacky way of getting consistent legend generation, also controls value rounding over the program
#     and some constant definitions

SWITCHBOARD.strALLACCESSIONS <- "All Accessions"
SWITCHBOARD.DIRECTORY <-
  "C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\"
SWITCHBOARD.strALLDATA <- "Entire"
SWITCHBOARD.roundto <- 3
SWITCHBOARD.strCLEANONLIST <-
  c("D_div_W")
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

SWITCHBOARD.strMODELLIST <-
  c(
    "Accession",
    "W",
    "H",
    "D",
    "Th",
    "WH",
    "WD",
    "WT",
    "HD",
    "HT",
    "DT",
    "WHD",
    "WHT",
    "WDT",
    "HDT",
    "HWDT"
  )

SWITCHBOARD.FRESH_DRY_RAWS <- list(
  c(9.17, 0.60),
  c(8.48, 0.49),
  c(8.73, 0.60),
  c(10.23, 0.68),
  c(5.76, 0.49),
  c(10.65, 0.89),
  c(10.37, 0.65),
  c(5.47, 0.40),
  c(8.11, 0.57),
  c(9.71, 0.66),
  c(4.56, 0.32),
  c(3.97, 0.33),
  c(7.59, 0.71),
  c(2.31, 0.21)
)

SWITCHBOARD.AVG_FRESH_DRY_WEIGHTS <-
  hash(SWITCHBOARD.strACCESSIONLIST, SWITCHBOARD.FRESH_DRY_RAWS)


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
    "FW_div_T",
    "dry_weight"
  )

#List of all plots  to be generated.
SWITCHBOARD.GRAPHS_LIST <- tribble(
  ~ xname, ~ yname,
  'width',     'fresh_weight',
  'height',    'fresh_weight',
  'diameter',  'fresh_weight',
  'thickness', 'fresh_weight',
  'D_div_W',   'fresh_weight',
  'width',     'height',
  'width',     'diameter',
  'width',     'thickness',
  'height',    'diameter',
  'height',    'thickness',
  'diameter',  'thickness',
  'D_div_W',   'dry_weight'
  #  'add_HW', 'thickness',
  #  'add_HW', 'fresh_weight',
  #  'add_HWT', 'fresh_weight',
  #  'add_DT', 'fresh_weight',
  #  'mul_HWDT', 'fresh_weight'
)

#R^2 Output csv dataframes
csvR2Frame <-
  data.frame(
    Dataset = character(),
    Graph = character(),
    R2 = double(),
    Mode = character()
  )

#List of anomalous datapads--TODO: Establish statistical, heuristic basis for removal
errorPads <-
  vector() #c("PARL 585-1-8", "PARL 319-1-1", "PARL 854-2-3", "PARL 580-1-2", "PARL 584-1-1", "PARL 242-4-3")

#FUNCTION DEFINITIONS

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
    'FW_div_H' = 'Fresh Weight / Height',
    'FW_div_W' = 'Fresh Weight / Width',
    'FW_div_D' = 'Fresh Weight / Diameter',
    'FW_div_T' = 'Fresh Weight / Thickness',
    "dry_weight" = "Dry Weight"
  )
  return(legname)
}

#Functions for error removal
#tagBySTDDEV--Called by ErrorCleaning before (in?) removeByID, this tags observations
# if they are more than n = 3 standard deviations away from the mean in a given parameter.
# However, may have to remove as fresh weight is not normally distributed.

ACCESSORY.tagBySTDDEV <-
  function(dataset, an_accession, threshold) {
    blacksite <- vector()
    #print(blacksite)
    dataset_accessionfilter <-
      filter(dataset, accession == an_accession)
    #for (columnname in c("H_div_W", "FW_div_W", "FW_div_D")) {
    for (columnname in SWITCHBOARD.strCLEANONLIST) {
      tempcol <-
        dataset_accessionfilter[[columnname]]
      meanCol <- mean(tempcol)
      STDCol <- sd(tempcol)
      bound_upper <- (meanCol + (STDCol * as.numeric(threshold)))
      bound_lower <- (meanCol - (STDCol * as.numeric(threshold)))
      detention <-
        filter(
          dataset_accessionfilter,
          !!as.symbol(columnname) < bound_lower |
            !!as.symbol(columnname) > bound_upper
        )
      blacklist <- as.vector(detention$completeID)
      blacksite <- c(blacksite, blacklist)
      #print(c("List:", blacklist))
      #print(c("Site:", blacksite))
      #print(c("Column:", columnname))
    }
    return(blacksite)
  }

#ErrorCleaning--
ACCESSORY.ErrorCleaning <- function(dataset, perc_threshold) {
  rawblock <- dataset
  dev_threshold <- qnorm((100 - ((
    100 - perc_threshold
  ) / 2)) / 100)
  for (accession in SWITCHBOARD.strACCESSIONLIST) {
    blacksite <-
      ACCESSORY.tagBySTDDEV(dataset, accession, dev_threshold)
    rawblock <-
      filter(rawblock,!(rawblock$completeID %in% blacksite))
  }
  finalDataset <- rawblock
  return(finalDataset)
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
                                dataset_in,
                                graphmode)
{
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
  if (graphmode == "Double Log") {
    targetline = lm(log(y_data) ~ log(x_data), data = subslice)
  } else if (graphmode == "Double Linear") {
    targetline = lm(y_data ~ x_data, data = subslice)
  }
  
  rsquare = rsq(targetline, type = "sse")
  intercept <- coefficients(targetline)[1]
  slope <- coefficients(targetline)[2]
  details <-
    paste0(
      "y_0 = ",
      round(intercept, SWITCHBOARD.roundto),
      ", m = ",
      round(slope, SWITCHBOARD.roundto),
      collapse = " "
    )
  
  csvR2Frame <<- rbind(
    csvR2Frame,
    data.frame(
      "Dataset" = paste0(year, " Dataset ", accession, collapse = " "),
      "Graph" = paste0(xlab, " vs. ", ylab,
                       collapse = " "),
      "R2" = round(rsquare, SWITCHBOARD.roundto),
      "Mode" = graphmode
    )
  )
  plot(
    x_data,
    y_data,
    type = "p",
    log = ifelse(graphmode == "Double Log", "xy", ""),
    cex = 0.8,
    sub = paste0(
      details,
      ", R^2 =",
      round(rsquare, SWITCHBOARD.roundto),
      ", N = ",
      nrow(subslice)
    ),
    main = paste0(
      xlab,
      " vs. ",
      ylab,
      "\n",
      year,
      " Dataset, ",
      accession,
      "Mode - ",
      graphmode,
      collapse = " "
    ),
    xlab = xlab,
    ylab = ylab
  )
  if (graphmode == "Double Log") {
    lines(x_data, exp(predict(targetline, newdata = list(x_data = x_data))) , col = "grey")
  } else if (graphmode == "Double Linear") {
    lines(x_data, predict(targetline, newdata = list(x_data = x_data)) , col = "grey")
  }
}

###Generator function which generates plots using Plotter over all of the axes combinations

PIPELINE.CorrPlotsGenerator <-
  function(figure_guide,
           year = SWITCHBOARD.strALLDATA,
           accession = SWITCHBOARD.strALLACCESSIONS,
           dataset_in) {
    for (i in 1:nrow(figure_guide)) {
      for (graphmode in c("Double Log", "Double Linear")) {
        row <- figure_guide[i, ]
        PIPELINE.PlotRegress(row$xname,
                             row$yname,
                             year,
                             accession,
                             dataset_in,
                             graphmode)
      }
    }
  }

ACCESSORY.qqGen <-
  function(accession, column, dataset_in, threshold) {
    subsetA <-
      ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA,
                           filter_accession = accession,
                           dataset_in = dataset_in)
    title <-
      paste0(accession, ", ", column)
    png(
      paste0(
        SWITCHBOARD.DIRECTORY,
        "QQPlot_Images",
        "\\",
        column,
        "--",
        accession,
        "_LIN_",
        threshold,
        "th.png"
      )
    )
    qqPlot(subsetA[[column]], main = title, envelope = 0.95)
    dev.off()
  }

ACCESSORY.qqGenLog <-
  function(accession, column, dataset_in, threshold) {
    subsetA <-
      ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA,
                           filter_accession = accession,
                           dataset_in = dataset_in)
    title <-
      paste0(accession, ", ", column)
    png(
      paste0(
        SWITCHBOARD.DIRECTORY,
        "QQPlot_Images",
        "\\",
        column,
        "--",
        accession,
        "_LOG_",
        threshold,
        "th.png"
      )
    )
    qqPlot(log(subsetA[[column]]), main = title, envelope = 0.95)
    dev.off()
  }

ACCESSORY.qqGenSqrt <-
  function(accession, column, dataset_in = dataset_in, threshold) {
    subsetA <-
      ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA, filter_accession = accession, dataset_in)
    title <-
      paste0(accession, ", ", column)
    png(
      paste0(
        SWITCHBOARD.DIRECTORY,
        "QQPlot_Images",
        "\\",
        column,
        "--",
        accession,
        "_SQRT_",
        threshold,
        "th.png"
      )
    )
    qqPlot(sqrt(subsetA[[column]]), main = title, envelope = 0.95)
    dev.off()
  }

PIPELINE.PlotCompiler <- function(dataset_in, threshold) {
  for (dataset in c(SWITCHBOARD.strALLDATA)) {
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
    }
  }
  
  for (accession in SWITCHBOARD.strACCESSIONLIST) {
    for (measure in SWITCHBOARD.strQQMEASURESLIST) {
      ACCESSORY.qqGen(accession, measure, dataset_in, threshold)
    }
  }
}

ACCESSORY.allSubsetsTable <- function (dataset) {
  models_df <-
    data.frame(Name = character(),
               AdjRsq = numeric(),
               SBC = numeric())
  iden_num <- 1
  for (w in c(FALSE, TRUE)) {
    for (h in c(FALSE, TRUE)) {
      for (d in c(FALSE, TRUE)) {
        for (t in c(FALSE, TRUE)) {
          iter_string <-
            paste0(
              ifelse(w, "width*", ""),
              ifelse(h, "height*", ""),
              ifelse(d, "diameter*", ""),
              ifelse(t, "thickness*", ""),
              sep = ""
            )
          iter_string <- substr(iter_string, 1, nchar(iter_string) - 1)
          if (iter_string != "") {
            iter_formula <- as.formula(paste0("fresh_weight ~ ", iter_string))
            iter_model <- lm(iter_formula, dataset)
            iter_summ <- summary(iter_model)
            models_df <- rbind(models_df,
                               data.frame(
                                 Name = iter_string,
                                 AdjRsq = iter_summ[["adj.r.squared"]],
                                 SBC = BIC(iter_model)
                               ))
          }
        }
      }
    }
  }
  return(models_df)
}

statsArray <- array(1:1350,
                    dim = c(3, 15, 15, 2),
                    dimnames = list(
                      c("100", "95", "90"),
                      c(
                        "All",
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
                      ),
                      c(
                        "width",
                        "height",
                        "diameter",
                        "thickness",
                        "width*height",
                        "width*diameter",
                        "width*thickness",
                        "height*diameter",
                        "height*thickness",
                        "diameter*thickness",
                        "width*height*diameter",
                        "width*height*thickness",
                        "width*diameter*thickness",
                        "height*diameter*thickness",
                        "width*height*diameter*thickness"
                      ),
                      c("AdjRsq", "SBC")
                    ))

main <- function() {
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  #-----------------------------------------------
  #_______________________________________________
  #DATASET MANIPULATION, FUNCTION EXECUTION SEGMENT
  #Load Parlier csv files into memory
  workingFilenames <- list()
  
  PARLIER <- SWITCHBOARD.csvMAINFILE
  
  #Data compatibility modifications
  PARLIER <- mutate(PARLIER, pad_id = 1)
  PARLIER$dataset <- "Year 1"
  
  #Add Error-correction derived measures
  
  PARLIER <- PARLIER %>%
    mutate(H_div_W = height / width) %>%
    mutate(FW_div_W = fresh_weight / width) %>%
    mutate(FW_div_D = fresh_weight / diameter) %>%
    mutate(FW_div_H = fresh_weight / height) %>%
    mutate(PARLIER, FW_div_T = fresh_weight / thickness) %>%
    mutate(PARLIER, D_div_W = diameter / width)
  
  areaFrame <-
    read.csv('Pad_Area_Estimations.csv', fileEncoding = 'UTF-8-BOM')
  PARLIER <- merge(PARLIER, areaFrame, by = "completeID")
  
  #Add Dry Weight measure for further analysis
  
  PARLIER$dry_weight <- 0
  for (accession in keys(SWITCHBOARD.AVG_FRESH_DRY_WEIGHTS)) {
    accession_set <-
      values(SWITCHBOARD.AVG_FRESH_DRY_WEIGHTS[accession])
    PARLIER$dry_weight[PARLIER$accession == accession] <-
      PARLIER$fresh_weight * (accession_set[2] / accession_set[1])
  }
  
  PARLIER <- PARLIER %>%
    mutate(Eccentricity = ((height^2 - width^2)/height)) %>%
    mutate(PartRatio = ((height-width)/(height+width))^2) %>%
    mutate(Pade_Peri = pi*(height+width)*(64-3*PartRatio)/(64-16*PartRatio)) %>%
    mutate(Pade_Derived_Diam = Pade_Peri/pi) %>%
    mutate(diamRatio = Pade_Derived_Diam/diameter)
  
  
  #print(PARLIER$fresh_weight*as.numeric(accession_set[3])/as.numeric(accession[1]))
  #print(PARLIER$dry_weight1)
  
  #print(PARLIER)
  
  for (threshold in c(100, 95, 90)) {
    parlVersion <- ACCESSORY.ErrorCleaning(PARLIER, threshold)
    print(head(parlVersion))    
    fileTitle <- paste0(SWITCHBOARD.DIRECTORY,
                        "ref_PARL0",
                        ifelse(threshold != 100, paste0("C_", threshold), ""),
                        ".csv")
    write.csv(parlVersion,
              fileTitle,
              row.names = FALSE)
    workingFilenames <- append(workingFilenames, fileTitle)
    #---------------------------------------------------------------
    for (accession in SWITCHBOARD.strACCESSIONLIST) {
      loadup_df <- ACCESSORY.allSubsetsTable(ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA, accession, parlVersion))
      for (model in c(
        "width",
        "height",
        "diameter",
        "thickness",
        "width*height",
        "width*diameter",
        "width*thickness",
        "height*diameter",
        "height*thickness",
        "diameter*thickness",
        "width*height*diameter",
        "width*height*thickness",
        "width*diameter*thickness",
        "height*diameter*thickness",
        "width*height*diameter*thickness"
      )) {
        for (stat in c("AdjRsq", "SBC")) {
          statsArray[as.character(threshold),
                     accession,
                     model,
                     stat] <- loadup_df[loadup_df$Name == model, stat]
        }
      }
    }
    #---------------------------------------------------------------
    
    pdftitle <- paste(
      "Double-Log Cladode Parameter Cross-Relations",
      ifelse(threshold != 100, paste0(" Corrected - ", threshold), ""),
      ".pdf",
      sep = ""
    )
    pdf(pdftitle)
    PIPELINE.PlotCompiler(parlVersion, threshold)
    dev.off()
    write.csv(
      csvR2Frame,
      paste(
        SWITCHBOARD.DIRECTORY,
        "R^2_Records\\Double-Log_Regression_R^2_Values",
        ifelse(threshold != 100, paste0("_Corrected_", threshold), ""),
        ".csv",
        sep = ""
      ),
      row.names = FALSE
    )
    csvR2Frame <<-
      data.frame(
        Dataset = character(),
        Graph = character(),
        R2 = double(),
        Mode = character()
      )
    
    statFrame <- data.frame(matrix(ncol = 16, nrow = 0))
    
    RsqFrame <- statFrame
    SBCFrame <- statFrame
    
    dimKeyPair <- list(
      "width" = "W",
      "height" = "H",
      "diameter" = "D",
      "thickness" = "Th",
      "width*height" = "WH",
      "width*diameter" = "WD",
      "width*thickness" = "WT",
      "height*diameter" = "HD",
      "height*thickness" = "HT",
      "diameter*thickness" = "DT",
      "width*height*diameter" = "WHD",
      "width*height*thickness" = "WHT",
      "width*diameter*thickness" = "WDT",
      "height*diameter*thickness" = "HDT",
      "width*height*diameter*thickness" = "WHDT"
    )
    
    modelaliases <- c(
      "width",
      "height",
      "diameter",
      "thickness",
      "width*height",
      "width*diameter",
      "width*thickness",
      "height*diameter",
      "height*thickness",
      "diameter*thickness",
      "width*height*diameter",
      "width*height*thickness",
      "width*diameter*thickness",
      "height*diameter*thickness",
      "width*height*diameter*thickness"
    )
    
    for (accession in SWITCHBOARD.strACCESSIONLIST) {
      RSQinsert <- c(accession)
      SBCinsert <- c(accession)
      for (model in c(1:length(modelaliases))) {
        RSQinsert[model + 1] <-
          statsArray[toString(threshold), toString(accession), modelaliases[model], "AdjRsq"]
        SBCinsert[model + 1] <-
          statsArray[toString(threshold), toString(accession), modelaliases[model], "SBC"]
      }
      
      RsqFrame <- rbind(RsqFrame, RSQinsert)
      SBCFrame <- rbind(SBCFrame, SBCinsert)
    }
    
    colnames(RsqFrame) <- SWITCHBOARD.strMODELLIST
    write.csv(
      RsqFrame,
      paste0(
        SWITCHBOARD.DIRECTORY,
        "Correlation_Analyses\\Rsq_T_",
        threshold,
        ".csv"
      ),
      row.names = FALSE
    )
    
    colnames(SBCFrame) <- SWITCHBOARD.strMODELLIST
    write.csv(
      SBCFrame,
      paste0(
        SWITCHBOARD.DIRECTORY,
        "Correlation_Analyses\\SBC_T_",
        threshold,
        ".csv"
      ),
      row.names = FALSE
    )
  }
}

main()
