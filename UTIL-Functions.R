#' @name UTIL.CollectFilterModelIDs
#' @description A test function designed to check all datapoints which pass the > 3p/n leverage threshold in any model, for each accession. A diagnostic tool.
#' @param data_df The source dataframe from which the models are derived.
#' @param fits_df The dataframe of models by accession.
#' @example UTIL.CollectFilterModelIDs(SWITCHBOARD.csvMAINFILE, DATA_OBJECTS.LinearRegressions_DF)

UTIL.CollectFilterModelIDs <- function(data_df, fits_df) {
  
  ID_message_collection <- c()
  ID_collection <- c()
  
  df_errors <- data.frame(ID = NA, Model = NA)
  df_errors <- df_errors[-c(1),]
  
  all_models <- names(fits_df[2:length(names(fits_df))])
  for (accession in SWITCHBOARD.ACCESSIONLIST) {
    for (model_name in all_models) {
      model_in <- CALC.ModelFetch(fits_df, accession, model_name)
      ID_vec <- model_in %>%
        CALC.ModelIndexHighLeverage() %>%
        CALC.ModelValuePairsFromIndex(model_in) %>%
        CALC.PadIDfromModelValuePairs(data_df, accession)
      for(ID in ID_vec) {
        df_errors[nrow(df_errors) + 1,] <- c(ID, model_name)
        ID_collection <- c(ID_collection, ID)
        string_out <- paste0("[", accession, "]{", model_name, "}: ", ID)
        #print(string_out)
        ID_message_collection <- c(ID_message_collection, string_out)
      }
    }
  }
  unique_IDs <- unique(ID_collection)
  unsorted_unique_ID_occurrences <- c()
  for(ID in unique_IDs) {
    occurrences <- length(which(ID_collection == ID))
    unsorted_unique_ID_occurrences <- c(unsorted_unique_ID_occurrences, paste0("[", occurrences, "]:", ID))
  }
  sorted_unique_ID_occurrences <- sort(unsorted_unique_ID_occurrences)
  df_errors <- df_errors %>%
    mutate(Accession = substring(ID,6,8),
           ID = as.factor(ID),
           Model = as.factor(Model)) %>% 
    add_count(ID)
  print(df_errors)
  return(df_errors)
}

#' @name UTIL.TagFilterByID 
#' @description 
#' @param df_in
#' @param df_errors
#' @example

UTIL.TagFilterByID <- function(df_in, df_errors) {
  df_out <- df_in %>%
    mutate(tagged_hiLeverage = F)
  df_out$tagged_hiLeverage[df_out$completeID %in% df_errors$ID] <- T
  head(df_out)
  print(nrow(df_out))
  return(df_out)
}

#' @name UTIL.AugmentModelWithID
#' @description
#' @param model_in
#' @param accession
#' @example

UTIL.AugmentModelWithID <- function(model_in, accession) {
  model_out <- model_in
  ID_vector <- c()
  ID_array <- array(c(1:32), dim = c(8, 4))
  ID_rows <- row(ID_array)
  ID_cols <- col(ID_array)
  for (n in c(1:32)) {
    ID_1 <- ID_cols[n]
    ID_2 <- ID_rows[n]
    ID_vector <- c(ID_vector, paste0("PARL ", accession, "-", ID_1, "-", ID_2))
  }
  model_out$model$ID <- ID_vector
  return(model_out)
}

#' @name UTIL.toModelString
#' @description Takes two strings and assembles them into a model name parseable by an lm() call.
#' @param xname A string denoting the independent variable.
#' @param yname A string denoting the denpendent variable.
#' @example UTIL.toModelString("width", "height") >> "height~width"

UTIL.toModelString <- function(xname, yname) {
  return(paste0(yname, "~", xname))
}

#' @name UTIL.quickPNG
#' @description
#' @param filenamestring
#' @param subfolder
#' @example

UTIL.quickPNG <- function(filenamestring, subfolder) {
  png(paste0(SWITCHBOARD.DIRECTORY, "\\", subfolder, "\\", filenamestring, ".png"))
}

#' @name UTIL.quickCSV
#' @description
#' @param df_csv
#' @param filenamstring
#' @param subfolder
#' @example

UTIL.quickCSV <- function(df_csv, filenamestring, subfolder) {
  write.csv(
    df_csv,
    paste0(SWITCHBOARD.DIRECTORY, "\\", subfolder, "\\", filenamestring, ".csv"),
    row.names = FALSE
  )
}

#'
#'
#'
#'
# UTIL.DropData <- function(data_in, data_drop_IDs) {
#   data_out <- data_in[-which(data_in$completeID %in% data_drop_IDs),]
#   return(data_out)
# }

#' @name UTIL.Collect_ID_NPN
#' @description
#' @param model_df
#' @param tag_text
#' @example

UTIL.Collect_ID_NPN <- function(model_df, tag_text) {
  ID_vec_out <- c()
  for (accession in rownames(model_df)) {
    for (modelname in colnames(model_df)[-1]) {
      model <- CALC.ModelFetch(model_df, accession, modelname)
      model_aug <- model %>%
        UTIL.AugmentModelWithLeverageTagging() %>%
        UTIL.AugmentModelWithID(accession)
      ID_index <- which(model_aug$model$tagged == tag_text)
      ID_vec <- model_aug$model$ID[ID_index]
      ID_vec_out <- unique(c(ID_vec_out, ID_vec))
    }
  }
  #print(str(ID_vec_out))
  return(ID_vec_out) #droplist
}

#' @name UTIL.DropOnNPN
#' @description
#' @param data_in
#' @param droplist_NPN
#' @example

UTIL.DropOnNPN <- function(data_in, droplist_NPN) {
  data_out <- data_in[-which(data_in$completeID %in% droplist_NPN),]
  #print(nrow(data_out))
  return(data_out)
}

#' @name UTIL.DropOnHybrid_3PN_and2PNVisual
#' @description
#' @param data_in
#' @param model_df
#' @param droplist_3PN
#' @param droplist_2PN
#' @example




# UTIL.DropOnHybrid_3PN_and2PNVisual <- function(data_in, model_df, droplist_3PN, droplist_2PN) {
#   check_2PN <- UTIL.Collect_ID_NPN(model_df, ">2p/n Threshold")
#   for (ID in droplist_2PN) {
#     if(!(ID %in% check_2PN)) {
#       print(paste0("Warning: ", ID, " is not in 2p/n list and will not be removed."))
#       droplist_2PN <- droplist_2PN[-which(droplist_2PN == ID)]
#     }
#   }
#   final_droplist <- unique(c(droplist_3PN, droplist_2PN))
#   data_out <- data_in[-which(final_droplist %in% data_in$completeID),]
#   print(nrow(data_out))
#   return(data_out)
# }

#' @name UTIL.AugmentModelWithLeverageTagging
#' @description
#' @param model_in
#' @example

UTIL.AugmentModelWithLeverageTagging <- function(model_in) {
  model_out <- model_in
  
  indices_HI <- model_out %>%
    CALC.ModelIndexHighLeverage()
  
  indices_MID <- model_out %>%
    CALC.ModelIndexMidLeverage() %>%
    CALC.ModelIndexMidLeverage_Exclusive(indices_HI)
  
  model_out$model$tagged <- "N/A"
  model_out$model$tagged[indices_MID] <- ">2p/n Threshold"
  model_out$model$tagged[indices_HI] <- ">3p/n Threshold"
  
  return(model_out)
}


UTIL.StoreCSV <- function(df_in, filenamestring, subfolder) {
  write.csv(df_in, paste0(SWITCHBOARD.DIRECTORY, "\\", subfolder, "\\", filenamestring, ".csv"), row.names = FALSE)
}

UTIL.AllAccessionToMain <- function(main_df, add_list) {
  rowaddloc <- nrow(main_df) + 1 
  main_df[rowaddloc, 1] <- "All_Accessions"
  for (model in names(add_list)) {
    coladdloc <- which(names(add_list) == model)
    main_df[rowaddloc, coladdloc+1] <- add_list[[model]]
  }
  
  return(main_df)
}


UTIL.StoreHeatmapsAndFile <- function(data_in, fill_models_df, filename, subfolder, abbreviations) {
  
  all_R2 <- CALC.ModelGenerator_AllAccessions_R2(data_in, fill_models_df)
  all_SBC <- CALC.ModelGenerator_AllAccessions_SBC(data_in, fill_models_df)
  
  R2_frame_name <- paste0(filename, "_R2")
  SBC_frame_name <- paste0(filename, "_SBC")
  
  # model_df <- CALC.ModelGenerator(data_in, fill_models_df)
  # R2_frame <- CALC.df_R2(model_df)
  # R2_frame <- R2_frame %>%
  #   UTIL.AllAccessionToMain(all_R2)
  # COMPILE.R2heatmap(R2_frame, filenamestring, subfolder, abbreviations)
  
  data_in %>% 
    CALC.ModelGenerator(fill_models_df) %>%
    CALC.df_R2 %>%
    UTIL.AllAccessionToMain(all_R2) %>%
    `row.names<-`(NULL) %T>%
    COMPILE.R2heatmap(R2_frame_name, subfolder, abbreviations) %>%
    UTIL.StoreCSV(R2_frame_name, subfolder)
  
  data_in %>% 
    CALC.ModelGenerator(fill_models_df) %>%
    CALC.df_SBC %>%
    UTIL.AllAccessionToMain(all_SBC) %>%
    `row.names<-`(NULL) %T>%
    COMPILE.SBCheatmap(SBC_frame_name, subfolder, abbreviations) %>%
    UTIL.StoreCSV(SBC_frame_name, subfolder)
}