#CAI ANALYSIS 2.0 SNIPPETS
#Scrap and test code used to diagnose program faults and debug.
#Every other file (CALC, UTIL, COMPILE, GRAPHER, and the main script) should be run before this if possible.

test_mult_model <- CALC.ModelFetch(DATA_OBJECTS.EllipticalRegressions_DF_UNFILT, 242, "fresh_weight~Theo_Area*thickness*Pade_Derived_Diam")
test_model <- CALC.ModelFetch(DATA_OBJECTS.LinearRegressions_DF, 242, "fresh_weight~height")

#]]DEBUG
#GRAPHER.tukeyHypergraph(CALC.tukeyTester(SWITCHBOARD.csvAUGMFILE, "height"), "height")

#]]DEBUG
#GRAPHER.intermeasure(test_aug_filt, DATA_OBJECTS.LinearRegressions_DF, 242, "height~width")

#]]DEBUG
#GRAPHER.R2heatmap(test_heatmap_in, "Test_Heatmap_R2", "Heatmaps")

#]]DEBUG
#GRAPHER.SBCheatmap(test_heatmap_in_SBC, "Test_Heatmap_SBC", "Heatmaps")

#]]DEBUG
#GRAPHER.tukeyBoxplot(SWITCHBOARD.csvAUGMFILE, "height")

# MAKER.BoxplotMaker <- function(dataset, measure, colors) {
#   TUKEY <- CALC.tukeyTester(dataset, measure)
#   LABELS <- CALC.tukeyLabel_DF(TUKEY)
#   target_boxplot <- boxplot(dataset[[measure]] ~ dataset$accession, ylim=c(min(dataset[[measure]]) , 1.1*max(dataset[[measure]])) , col=colors[as.numeric(LABELS[,1])] , ylab="value" , main="")
#   target_boxplot
#   return(target_boxplot)
# }

#MAKER.BoxplotMaker(SWITCHBOA RD.csvAUGMFILE, "height", EXAMPLE_COLORS)




# # MAINFILE | QQ Plots
# cat("::| Compiling QQ Plots (MAINFILE)\n")
# COMPILE.qqGen(SWITCHBOARD.csvMAINFILE)
# 
# # INTERMEASURE | Intermeasure
# cat("::| Compiling Intermeasure Plots (MAINFILE)\n")
# COMPILE.intermeasure(DATA_OBJECTS.LinearRegressions_DF)
# 
# # MAINFILE | Tukey Intercomparison Plots
# cat("::| Compiling Tukey Intercomparison Plots (MAINFILE)\n")
# COMPILE.tukeygraphs(SWITCHBOARD.csvMAINFILE)




# Rejected Code -----------------------------------------------------------

DATA_OBJECTS.blacklist_2PN_RAWHAND <- c("PARL 319-4-8",
                                        "PARL 325-4-4",
                                        "PARL 325-2-4",
                                        "PARL 390-3-7",
                                        "PARL 390-4-7",
                                        "PARL 572-1-8",
                                        "PARL 572-3-3",
                                        "PARL 580-3-5",
                                        "PARL 584-3-1",
                                        "PARL 839-4-5")

DATA_OBJECTS.blacklist_2PN_CHECK <- UTIL.Collect_ID_NPN(DATA_OBJECTS.LinearRegressions_DF[,-which(names(DATA_OBJECTS.LinearRegressions_DF) %in% c("thickness~width", "thickness~height", "thickness~diameter"))], ">2p/n Threshold")
DATA_OBJECTS.csvAUGM_2PN <- UTIL.DropOnNPN(DATA_OBJECTS.csvAUGM_3PN, DATA_OBJECTS.blacklist_2PN_CHECK)



#' @name UTIL.DropOnHybrid_3PN_and2PNVisual
#' @description
#' @param data_df
#' @param model_df
#' @param droplist_3PN
#' @param droplist_2PN
#' @example

# UTIL.DropOnHybrid_3PN_and2PNVisual <- function(data_df, model_df, droplist_3PN, droplist_2PN) {
#   check_2PN <- UTIL.Collect_ID_NPN(model_df, ">2p/n Threshold")
#   for (ID in droplist_2PN) {
#     if(!(ID %in% check_2PN)) {
#       print(paste0("Warning: ", ID, " is not in 2p/n list and will not be removed."))
#       droplist_2PN <- droplist_2PN[-which(droplist_2PN == ID)]
#     }
#   }
#   final_droplist <- unique(c(droplist_3PN, droplist_2PN))
#   data_out <- data_df[-which(final_droplist %in% data_df$completeID),]
#   print(nrow(data_out))
#   return(data_out)
# }

#' @name UTIL.CollectFilterModelIDs
#' @description A test function designed to check all datapoints which pass the > 3p/n leverage threshold in any model, for each accession. A diagnostic tool.
#' @param data_df The source dataframe from which the models are derived.
#' @param fits_df The dataframe of models by accession.
#' @example UTIL.CollectFilterModelIDs(SWITCHBOARD.csvMAINFILE, DATA_OBJECTS.LinearRegressions_DF)

# UTIL.CollectFilterModelIDs <- function(data_df, fits_df) {
#   
#   ID_message_collection <- c()
#   ID_collection <- c()
#   
#   df_errors <- data.frame(ID = NA, Model = NA)
#   df_errors <- df_errors[-c(1),]
#   
#   all_models <- names(fits_df[2:length(names(fits_df))])
#   for (accession in SWITCHBOARD.ACCESSIONLIST) {
#     for (model_name in all_models) {
#       model_in <- CALC.ModelFetch(fits_df, accession, model_name)
#       ID_vec <- model_in %>%
#         CALC.ModelIndexHighLeverage() %>%
#         CALC.ModelValuePairsFromIndex(model_in) %>%
#         CALC.PadIDfromModelValuePairs(data_df, accession)
#       for(ID in ID_vec) {
#         df_errors[nrow(df_errors) + 1,] <- c(ID, model_name)
#         ID_collection <- c(ID_collection, ID)
#         string_out <- paste0("[", accession, "]{", model_name, "}: ", ID)
#         #print(string_out)
#         ID_message_collection <- c(ID_message_collection, string_out)
#       }
#     }
#   }
#   unique_IDs <- unique(ID_collection)
#   unsorted_unique_ID_occurrences <- c()
#   for(ID in unique_IDs) {
#     occurrences <- length(which(ID_collection == ID))
#     unsorted_unique_ID_occurrences <- c(unsorted_unique_ID_occurrences, paste0("[", occurrences, "]:", ID))
#   }
#   sorted_unique_ID_occurrences <- sort(unsorted_unique_ID_occurrences)
#   df_errors <- df_errors %>%
#     mutate(Accession = substring(ID,6,8),
#            ID = as.factor(ID),
#            Model = as.factor(Model)) %>% 
#     add_count(ID)
#   print(df_errors)
#   return(df_errors)
# }

#' @name UTIL.TagFilterByID
#' @description
#' @param df_in
#' @param df_errors
#' @example

# UTIL.TagFilterByID <- function(df_in, df_errors) {
#   df_out <- df_in %>%
#     mutate(tagged_hiLeverage = F)
#   df_out$tagged_hiLeverage[df_out$completeID %in% df_errors$ID] <- T
#   head(df_out)
#   print(nrow(df_out))
#   return(df_out)
# }


 
#' @name UTIL.Collect_ID_NPN
#' @description Collects 
#' @param model_df
#' @param tag_text
#' @example

# UTIL.Collect_ID_NPN <- function(model_df, tag_text) {
#   ID_vec_out <- c()
#   for (accession in rownames(model_df)) {
#     for (modelname in colnames(model_df)[-1]) {
#       model <- CALC.ModelFetch(model_df, accession, modelname)
#       model_aug <- model %>%
#         UTIL.AugmentModelWithLeverageTagging() %>%
#         UTIL.AugmentModelWithID(accession)
#       ID_index <- which(model_aug$model$tagged == tag_text)
#       ID_vec <- model_aug$model$ID[ID_index]
#       ID_vec_out <- unique(c(ID_vec_out, ID_vec))
#     }
#   }
#   #print(str(ID_vec_out))
#   return(ID_vec_out) #droplist
# }

#' @name CALC.ModelIndexMidLeverage
#' @description Given a linear model, returns a vector of indices of data points whose leverage value exceeds a 2p/n threshold,
#'   where p is the number of parameters in the model, 
#'   and n is the number of data points in the model.
#' @param model_in A linear model object.

CALC.ModelIndexMidLeverage <- function(model_in) {
  threshold <- 2*length(model_in$coefficients)/length(model_in$model[[1]])
  result_index <- as.vector(which(hatvalues(model_in) > threshold))
  return(result_index)
}

#' @name CALC.ModelIndexMidLeverage_Exclusive
#' @description Given two index lists of 2p/n and 3p/n threshold values, returns solely the indexes of datapoints which
#'   pass the 2p/n threshold without passing the 3p/n threshold. See CALC.ModelNumHighLeverage and CALC.ModelNumMidLeverage
#' @param result_index_MID A vector list of indices, typically generated by CALC.ModelIndexMidLeverage
#' @param result_index_HIGH A vector list of indices, typically generated by CALC.ModelIndexHighLeverage

CALC.ModelIndexMidLeverage_Exclusive <- function(result_index_MID, result_index_HIGH) {
  return(result_index_MID[-which(result_index_HIGH %in% result_index_MID)])
}


#' @name UTIL.AugmentModelWithLeverageTagging
#' @description Adds data on whether or not each datapoint exceeds the 3P/N leverage threshold for the particular model.
#' @param model_in The linear model object to be augmented.
#' @returns An augmented linear model object with a tagged column within the model_obj$model property.

UTIL.AugmentModelWithLeverageTagging <- function(model_in) {
  model_out <- model_in
  
  indices_HI <- model_out %>%
    CALC.ModelIndexHighLeverage()
  
  # indices_MID <- model_out %>%
  #   CALC.ModelIndexMidLeverage() %>%
  #   CALC.ModelIndexMidLeverage_Exclusive(indices_HI)
  
  model_out$model$tagged <- "N/A"
  #model_out$model$tagged[indices_MID] <- ">2p/n Threshold"
  model_out$model$tagged[indices_HI] <- ">3p/n Threshold"
  
  return(model_out)
}

