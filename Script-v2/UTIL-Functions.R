#' @name UTIL.toModelString
#' @description Takes two strings and assembles them into a model name parseable by an lm() call.
#' @param xname A string denoting the independent variable.
#' @param yname A string denoting the denpendent variable.
#' @example UTIL.toModelString("width", "height") >> "height~width"
#' @returns A string containing a modelname.

UTIL.toModelString <- function(xname, yname) {
  return(paste0(yname, "~", xname))
}

# Filtration Utilities ----------------------------------------------------

#' @name UTIL.AugmentModelWithID
#' @description Adds IDs for each datapoint with a given accession number.
#' @param model_in The linear model object to be augmented.
#' @param accession A string containing the accession ID for the augmentation.
#' @returns An augmented linear model object with a ID column within the model_obj$model property.
 
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

#' @name Assemble3PNDroplist
#' @description Collects all datapoints which exceed the 3P/N leverage threshold in the models selected, by accession.
#' @param fits_df A dataframe containing pregenerated models from CALC.ModelGenerator.
#' @param data_df A dataframe containing all measured data points with IDs and their measures.
#' @returns A vector with string IDs.

UTIL.Assemble3PNDroplist <- function(fits_df, data_df) {
  blacklist_3pn <- c()

  for (modelname in colnames(fits_df)[-1]) {

    for (accession in SWITCHBOARD.ACCESSIONLIST) {

      model <- CALC.ModelFetch(fits_df, accession, modelname)
      target_loc <- CALC.ModelIndexHighLeverage(model)

      for (loc in target_loc) {
        pairframe <- CALC.ModelValuePairsFromIndex(loc, model)
        modelID <- CALC.PadIDfromModelValuePairs(pairframe, data_df, accession)
        blacklist_3pn <- c(blacklist_3pn, modelID)
      }
    }
  }
  blacklist_3pn <- sort(unique(blacklist_3pn))
  return(blacklist_3pn)
}

#' @name UTIL.DropOnNPN
#' @description Removes IDs in a droplist from a larger dataset.
#' @param data_df A dataframe containing all measured data points with IDs and their measures.
#' @param droplist_NPN A vector of strings containing IDs of the datapoints to drop.
#' @returns A dataframe excluding the removed datapoints.

UTIL.DropOnNPN <- function(data_df, droplist_NPN) {
  data_out <- data_df[-which(data_df$completeID %in% droplist_NPN),]
  #print(nrow(data_out))
  return(data_out)
}

#' @name UTIL.AugmentModelWithLeverageTagging
#' @description Adds data on whether or not each datapoint exceeds the 3P/N leverage threshold for the particular model.
#' @param model_in The linear model object to be augmented.
#' @returns An augmented linear model object with a tagged column within the model_obj$model property.

UTIL.AugmentModelWithLeverageTagging <- function(model_in) {
  model_out <- model_in
  
  indices_HI <- model_out %>%
    CALC.ModelIndexHighLeverage()
  
  model_out$model$tagged <- "N/A"
  model_out$model$tagged[indices_HI] <- ">3p/n Threshold"
  
  return(model_out)
}

#' @name UTIL.AllAccessionToMain
#' @description Appends a row of summary statistic values (for models applied over the entire dataset) to a dataframe containing the by-accession table of values.
#' @param model_df The summary statistic dataframe to be augmented.
#' @param add_list A list of values to be appended to the dataframe.
#' @returns A dataframe with an additional row and the same amount of columns.

UTIL.AllAccessionToMain <- function(model_df, add_list) {
  rowaddloc <- nrow(model_df) + 1 
  model_df[rowaddloc, 1] <- "All"
  for (model in names(add_list)) {
    coladdloc <- which(names(add_list) == model)
    model_df[rowaddloc, coladdloc+1] <- add_list[[model]]
  }
  
  return(model_df)
}


# File Creation Utilities -------------------------------------------------

#' @name UTIL.FolderPathCheck
#' @description Checks to see if a folder exists, and makes it if it doesn't.
#' @param subfolder_string A string or path containing a valid subfolder.

UTIL.FolderPathCheck <- function(subfolder_string) {
  path_components <- str_split(subfolder_string, "\\\\")[[1]]
  for (position in c(1:length(path_components))) {
    checkpath <- path_components[1:position] %>%
      paste0(collapse = "\\\\")
    if(!(file.exists(checkpath))) {
      dir.create(file.path(checkpath))
    } 
  }
}

#' @name UTIL.quickPNG
#' @description Creates a .png object with a given name to store an image in.
#' @param filenamestring The name given for the .png file.
#' @param subfolder The folder that the file should be stored in.

UTIL.quickPNG <- function(filenamestring, subfolder, dimensions = c(600, 4, 4)) {
  
  png_ppi = dimensions[1]
  png_height = dimensions[2]
  png_width = dimensions[3]
  
  UTIL.FolderPathCheck(subfolder)
  filepath <- file.path(SWITCHBOARD.DIRECTORY, subfolder, filenamestring, fsep = "\\")
  
  png(paste0(filepath, ".png"), 
      width = png_width * png_ppi, 
      height = png_height * png_ppi, 
      res = png_ppi)
}

#' @name UTIL.quickCSV
#' @description Creates a .csv file from a dataframe and stores it in a particular path.
#' @param data_df A dataframe containing all measured data points with IDs and their measures.
#' @param filenamestring The name given for the .png file.
#' @param subfolder The folder that the file should be stored in.

UTIL.quickCSV <- function(data_df, filenamestring, subfolder) {
  
  UTIL.FolderPathCheck(subfolder)
  filepath <- file.path(SWITCHBOARD.DIRECTORY, subfolder, filenamestring, fsep = "\\")
  
  write.csv(
    data_df,
    paste0(filepath, ".csv"), 
    row.names = FALSE
  )
}

#' @name UTIL.StoreHeatmapsAndFile
#' @description Generates and stores a Heatmap image as a .png and its corresponding table as a .csv file.
#' @param data_df A dataframe containing all measured data points with IDs and their measures.
#' @param fill_models_df A vector of model strings.
#' @param filenamestring The name given for the .png file.
#' @param subfolder The folder that the file should be stored in. 
#' @param abbreviations A vector of string abbreviations for the models for intelligible graphing.

UTIL.StoreHeatmapsAndFile <- function(data_df, fill_models_df, filenamestring, subfolder, heatmap_title, abbreviations) {
  
  if (!(dir.exists(subfolder))) {
    dir.create(subfolder)
  }
  
  all_R2 <- CALC.ModelGenerator_AllAccessions_R2(data_df, fill_models_df)
  all_SBC <- CALC.ModelGenerator_AllAccessions_SBC(data_df, fill_models_df)
  
  R2_frame_name <- paste0(filenamestring, "_R2")
  SBC_frame_name <- paste0(filenamestring, "_SBC")
  
  data_df %>% 
    CALC.ModelGenerator(fill_models_df) %>%
    CALC.df_R2 %>%
    UTIL.AllAccessionToMain(all_R2) %>%
    `row.names<-`(NULL) %T>%
    COMPILE.R2heatmap(R2_frame_name, paste0(subfolder, "\\Graphs"), heatmap_title, abbreviations) %>%
    UTIL.quickCSV(R2_frame_name, paste0(subfolder, "\\Tables"))
  
  data_df %>% 
    CALC.ModelGenerator(fill_models_df) %>%
    CALC.df_SBC %>%
    UTIL.AllAccessionToMain(all_SBC) %>%
    `row.names<-`(NULL) %T>%
    COMPILE.SBCheatmap(SBC_frame_name, paste0(subfolder, "\\Graphs"), heatmap_title, abbreviations) %>%
    UTIL.quickCSV(SBC_frame_name, paste0(subfolder, "\\Tables"))
}

#' @name UTIL.SplitModelListOnInteractions
#' @description Splits a master model list into two lists with separate + and * models for more visually accessible heatmap graphing.
#' @param model_list A vector of model strings.
#' @param num_keep THe number of initial models to put in both output lists.
#' @returns A list with two string vectors.

UTIL.SplitModelListOnInteractions <- function(model_list, num_keep) {
  #print(model_list)
  keep_names <- model_list$models[1:num_keep]
  
  #no_interactions_piece <- model_list[which("+" %in% model_list$models),]
  #interactions_piece <- model_list[which("*" %in% model_list$models),]
  
  no_interactions_piece <- model_list[which(model_list$models %like% "[+]"),]
  interactions_piece <- model_list[which(model_list$models %like% "[*]"),]
  
  no_inter_list <- as.data.frame(c(keep_names, no_interactions_piece)) %>%
    `colnames<-`(c("models"))
  inter_list <- as.data.frame(c(keep_names, interactions_piece)) %>%
    `colnames<-`(c("models"))
  
  return(list(no_inter_list, inter_list))
}