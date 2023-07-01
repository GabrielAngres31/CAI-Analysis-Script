# Set repo
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

library(devtools)
# Package list
#install_version("vctrs", version = "0.5.2", repos = r["CRAN"])

REQUIRED_LIBRARIES <- c("dplyr", 
                        "tibble", 
                        "rsq", 
                        "car", 
                        "Rcmdr",
                        "hash",
                        "multcompView",
                        "HyperG",
                        "htmlwidgets", 
                        "magrittr",
                        "multiApply",
                        "stringr",
                        "viridis",
                        "RColorBrewer",
                        "data.table",
                        "rgl",
                        "knitr",
                        "rglwidget") 

REQUIRED_LIBRARIES <- c("dplyr", "tibble", "rsq", "car", "Rcmdr", "hash", "multcompView", "HyperG", "plotly", "htmlwidgets", "magrittr", "tidyr", "multiApply", "stringr", "viridis", "RColorBrewer", "data.table")

# library("dplyr")
# library("tibble")
# library("rsq")
# library("car")
# library("Rcmdr")
# library("hash")
# library("multcompView")
# library("HyperG")
# library("plotly")
# library("htmlwidgets")
# library("magrittr")
# library("tidyr")
# library("multiApply")
# library("stringr")
# library("viridis")
# library("pheatmap")
# library("RColorBrewer")
# library("data.table")

require(devtools)
#install_version("ggplot2", version = "0.9.1", repos = r["CRAN"])

# Install packages not yet installed
installed_packages <- REQUIRED_LIBRARIES %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(REQUIRED_LIBRARIES[!installed_packages])
}

# Packages loading
invisible(lapply(REQUIRED_LIBRARIES, library, character.only = TRUE))

# Working directory check
if (!(basename(getwd()) == "CAI-Analysis-Script")) {
  cat("Please make sure you have set your working directory to the \"CAI-Analysis-Script\" folder.\n")
} else {
  cat("Working directory verified as \"CAI-Analysis-Script\".\n")
}


# Data Loading -----------------

# Loading Dataset
# NOTE: All dataset augmentation must take place before the program is run.
#       The user must check that all data necessary for the analysis has been generated.


# Checks that a file path exists, and creates it if it does not.
UTILITY.folderpathcheck <- function(subfolder_string, silenced = FALSE) {
  # Get all parts of the file path to check
  path_components <- str_split(subfolder_string, "\\\\")[[1]]
  
  # Go through the entire file path, subdirectory by subdirectory
  for (position in c(1:length(path_components))) {
    
    # Iteratively build up the file path
    checkpath <- path_components[1:position] %>%
      paste0(collapse = "\\\\")
    
    # If the path doesn't exist, make it
    if(!(file.exists(checkpath))) {
      dir.create(file.path(checkpath))
      if(!silenced) {cat(paste0("\"", checkpath, "\" has been CREATED...\n"))}
    } 
  }
  
  # Affirm that the request file path exists
  if(!silenced) {cat(paste0("\"", subfolder_string, "\" exists\n"))}
}

# Creates a PNG at a given file path with default dimensions 600ppi, 4 inches height and width.
# Relies on there being an existing plot object to store.
UTILITY.quickPNG <- function(filenamestring, subfolder, dimensions = c(600, 4, 4)) {
  
  # Set the dimensions of the image to store.
  # This should be configured according to the requirements of the submission journal
  png_ppi = dimensions[1]
  png_height = dimensions[2]
  png_width = dimensions[3]
  
  # Check that the filepath for the PNG exists, and make it if it doesn't
  UTILITY.folderpathcheck(subfolder, silenced = TRUE)
  filepath <- file.path(subfolder, filenamestring, fsep = "\\")
  
  # Create the PNG
  png(paste0(filepath, ".png"), 
      width = png_width * png_ppi, 
      height = png_height * png_ppi, 
      res = png_ppi)
}

# Check if the dataframe contains the columns needed for the analysis
UTILITY.columncheck <- function(columns_for_analysis, source_data = DATA.unfiltered, display = FALSE) {
  
  # Check whether the required columns (passed as a parameter) are named columns in the dataframe
  if (all(columns_for_analysis %in% colnames(source_data))) {
    
    # Output a diagnostic message if "display" is TRUE
    if (display) {
      cat("All required columns are present in the dataframe: \"")
      cat(columns_for_analysis)
      cat("\"\n")
    }
  } else { # Break out of the program with a warning message
    cat("Some or all of the required columns are missing from the dataframe\n")
    stop(paste0("Please check that all required columns [", columns_for_analysis, "] are present in your input file..."))
  }
}

# Checks to see if the given data contains valid accessions.
# Data must have an accession column for this function to work.
UTILITY.accessioncheck <- function(target_accessions, source_data = DATA.unfiltered) {
  
  # Get the name of the function that invoked this one
  calling_function_name <- deparse(sys.calls()[[sys.nframe()-1]])
  
  # Check that accession information is available in the data
  UTILITY.columncheck(c("accession"), source_data)
  
  # Check for a valid accession parameter
  accessions_possible <- source_data$accession %>% unique %>% sort
  
  # If an invalid accession is passed, stop the program with a custom error message with the specific function call.
  if (!(all(target_accessions %in% c("All", accessions_possible)))) {
    stop(paste0("MESSAGE: Please check that you have a valid input for \"accession\" for ", calling_function_name, ". You passed [", target_accessions, "]\n"))
  }
} 


# Data Filtration -----------------


# This function takes an list of hatvalues, given the size of the data and the number of parameters,
#   and returns a vector of Booleans which state whether each datapoints exceeds a threshold of 3P/N.
GENERATOR.leverageFlagger <- function(model_object, leverage_constant = 3) {
  
  # Get the hatvalues of the model object
  model_hatvalues <- hatvalues(model_object)
  
  # Get the number of parameters used in the model (P)
  num_params <- model_object %>%
    coef %>%
    length
  
  # Get the number of datapoints that the model is fit over (N)
  data_size <- model_object %>%
    .$model %>%
    nrow
  
  # Produce a boolean vector showing whether particular datapoints in the model exceed the leverage value
  leverageDecision <- (model_hatvalues > leverage_constant*num_params/data_size)
  return(leverageDecision)
}

# This function takes a data file and certain parameters for filtering:
#   Taking a set of models to evaluate individual leverage on for each accession, 
#   this function returns a list of all the datapoints that produced at least [max_flags] flags (1 by default)
#   using a filter on leverage (3P/N by default) using GENERATOR.leverageFlagger()
UTILITY.datapointFlagger <- function(max_flags = 1, leverage_constant = 3, source_data = DATA.unfiltered, file_models_for_filtering = "leverage_filter_models.txt") {
  # Get possible accessions
  accessions_possible <- source_data$accession %>% unique %>% sort
  
  # Generate Filtered Dataset
  
  #   Get the models to filter on for leverage criteria
  models_for_filtering <- scan(paste0("parameters\\", file_models_for_filtering), character())
  #   Prepare a dataframe to store specimen IDs and the model that they had a too-high leverage value for
  table_of_errors_by_model <- data.frame(matrix(nrow = 0, ncol = 2), row.names = NULL) %>%
    setnames(c("completeID", "model_for_filtering"))
  
  #   Fill the dataframe
  for (model in models_for_filtering) {
    for (acc in accessions_possible) {
      leverage_decision <- GENERATOR.modelFinder(source_data, model, acc, return_lm_object = TRUE) %>%
        GENERATOR.leverageFlagger(leverage_constant)
      flagged_IDs <- source_data %>%
        filter(accession == acc) %>%
        .$completeID %>%
        .[leverage_decision]
      
      for (ID in flagged_IDs) {
        table_of_errors_by_model[nrow(table_of_errors_by_model)+1,] <- list(completeID = ID, model_for_filtering = model)
      }
      
    } 
  }
  
  #   Get list of unique specimen IDs
  unique_problem_IDs <- table_of_errors_by_model$completeID %>%
    unique()
  
  #   Prepare matrix-style dataframe showing whether a particular unique point was flagged in a given model
  table_of_errors_expanded <- data.frame(matrix(0, nrow = length(unique_problem_IDs), ncol = length(models_for_filtering)), row.names = unique_problem_IDs) %>%
    setnames(c(models_for_filtering))
  
  #   Fill the dataframe
  for(row in 1:nrow(table_of_errors_by_model)) {
    problem_ID <- table_of_errors_by_model[row, "completeID"]
    model_problem <- table_of_errors_by_model[row, "model_for_filtering"]
    
    table_of_errors_expanded[problem_ID, model_problem] <- 1
  }
  
  #   Calculate the number of models that each ID was flagged in...
  table_of_errors_expanded <- table_of_errors_expanded %>% 
    mutate(sum_flags = rowSums(.)) %>%
    #   ...and if the max flags is below the threshold max_flags, it's kept on the list of outliers to disclude.
    filter(sum_flags >= max_flags)
  
  # Get the points to be excluded
  data_blacklist <- table_of_errors_expanded %>% 
    .$sum_flags %>%
    names()
  
  # Return the list
  return(data_blacklist)
}

# Model Generation ----------------

GENERATOR.modelFinder <- function(source_data, model_string, accession_choice = "All", return_params = c("All"), return_lm_object = FALSE) {
  
  # Check that accession information is available in the data
  UTILITY.columncheck(c("accession"), source_data)
  
  # Check for a valid accession parameter
  accessions_possible <- source_data$accession %>% unique %>% sort
  
  # Check if the accessions passed are valid and stop the execution if they aren't
  UTILITY.accessioncheck(accession_choice)
  
  #TODO: Check that model is correct, and stop it if it isn't!
  
  # Ensure that "accession" is a vector with one element
  accession_choice <- c(accession_choice)
  
  # Select data to model on, from source_data
  if (identical(accession_choice, "All")) {
    selected_data <- source_data
  } else {
    selected_data <- source_data %>%
      filter(accession == accession_choice)
  }
  
  # Fit a Linear Model to the selected data, according to the given model
  model_object <- lm(model_string, selected_data)
  
  # If the function call specified returning the linear model object,
  #   we can just return the lm() object and break out of the function.
  if (return_lm_object) {
    return(model_object)
  }
  
  # Extract various characteristics of the model
  #   Establish valid parameters
  valid_parameters = c("R^2", "Adj. R^2", "SBC", "Leverage", "P-value", "Coefficients", "Model")
  
  #   Stop the function execution if the user passes an invalid parameter
  if (identical(return_params, c("All"))) {
    return_params <- valid_parameters
  }
  if (!(all(return_params %in% valid_parameters))) {
    stop(paste0("MESSAGE: One or more parameters in your call of modelFinder were invalid. You passed: ", return_params))
  }
  
  #   Generate all possible valid parameters
  #     Generate the model summary statistics
  model_summary  <- summary(model_object)
  
  #     R^2
  model_R2       <- model_summary$r.squared
  
  #     Adjusted R^2
  model_adj_R2   <- model_summary$adj.r.squared
  
  #     SBC
  model_SBC      <- BIC(model_object)
  
  #     List of leverage values with IDs
  model_leverage <- hatvalues(model_object)
  
  #     p-value
  model_pvalue   <- coefficients(model_summary)[,4]
  
  #     list of fitted model parameters with names
  model_coeff    <- coefficients(model_summary)[,1]
  
  #     Full model data
  model_data     <- model_object$model
  
  #     X and Y data for the regression
  X_data <- model_object$model[[2]]
  Y_data <- model_object$model[[1]]
  
  #     Names of dependent and independent variables
  X_name <- model_object$model %>% names() %>% .[[2]]
  Y_name <- model_object$model %>% names() %>% .[[1]]
  
  # Put them in a named list
  return_values_full <- list(model_R2, model_adj_R2, model_SBC, model_leverage, model_pvalue, model_coeff, model_data)
  names(return_values_full) <- valid_parameters
  
  # Subset the list according to what the user asked for
  #   Returns all parameters by default
  subset_request <- return_values_full[return_params]
  
  # Output the desired parameters
  return(subset_request)
  
}


# Plot Makers --------------

# This function makes a QQ plot of a single measure by accession.
# NOTE: It doesn't make sense to make a qqplot of all the data given that they have different distributions, 
#       so passing "All" into this function for the accession value is not possible.
GENERATOR.qqGen <- function(measure, target_accession, source_data = DATA.unfiltered) {
  # Select a particular accession from the set
  target_data <- source_data %>% filter(accession == target_accession) #[[measure]][source_data$accession == accession]
  # Generate title text for the qqPlot
  qq_plot_name <- paste0(measure, "--", target_accession, "_QQ-Plot")
  # Create a PNG at the appropriate folder using UTILITY.quickPNG
  UTILITY.quickPNG(qq_plot_name, "image_output\\qq_plots")
  # Generate the qqPlot
  qqPlot(target_data[[measure]], ylab = measure, main = paste0("QQPlot of ", measure, " for accession ", target_accession), envelope = 0.95)
  # Turn off the plotting device
  dev.off()
}

# This function takes an input dataframe and a choice of measures, and stores a violin plot of the measure across each accession.
GENERATOR.violin_plot <- function(plot_data, measure, tag) {
  
  # Prepare file to write plot image to
  plot_filename  <- paste0(measure, "_distribution_", tag)
  plot_subfolder <- "image_output\\violin_plots"
  UTILITY.quickPNG(plot_filename, plot_subfolder)
  
  # Generate the violin plot file in ggplot
  violin_plot <- ggplot(plot_data, aes(get(measure), as.factor(accession), fill = as.factor(accession))) +
    geom_violin(trim=F, show.legend = F) +
    geom_boxplot(width = 0.2) +
    xlab(measure) +
    ylab("Accession") + 
    ggtitle(paste0("Distribution of ", measure, " by Accession")) + 
    theme(legend.position = "none")
  
  # Plot the plot
  plot(violin_plot)
  
  # Close graphing device
  dev.off()
}

# This function take a source data file and a dataframe containing models and their abbreviations for display, 
#   and returns a heatmap with accessions for the rows and abbreviated model names for the columns

#TODO: Find out why passing default parameter results in unexpected behavior
GENERATOR.heatmap_table <- function(models_and_abbreviations_df, model_statistic, source_data) {
  
  # Verify that the input models dataframe only contains the model and abbreviation columns that we need
  models_and_abbreviations_df %<>% select(c("model", "abbreviation"))
  
  # This vector must be adjusted by the user if they wish to add features that can be used in a heatmap.
  allowed_parameters <- c("R^2", "Adj. R^2", "SBC")
  
  # Verify that model_statistic is valid (i.e. a single-value measure)
  if (!(model_statistic %in% allowed_parameters)) {
    stop(paste0("You have put in an invalid measure. Currently, only [", allowed_parameters, "] are permitted. You passed [", model_statistic, "]"))
  }
  
  # Get valid accession
  accessions_possible <- source_data$accession %>% unique %>% sort %>% c("All", .)
  
  # Get models and abbreviations as separate vectors
  models <- models_and_abbreviations_df$model
  abbreviations <- models_and_abbreviations_df$abbreviation
  
  # Create a named list (abbreviations named by their full model text) for easy iteration
  models_abbr_list <- as.list(abbreviations) %>%
    `setNames`(models)
  
  # Initialize the return dataframe
  heatmap_table <- matrix(0, nrow = length(accessions_possible), ncol = length(models), dimnames = list(accessions_possible, abbreviations)) %>%
    data.frame()
  
  #Fill the dataframe with the given value
  for (model in models) {
    for (accession in accessions_possible) {
      # Get the appropriate abbreviation for the model
      abbr <- models_abbr_list[[model]]
      
      # Get the required value for the particular accession and model
      value <- GENERATOR.modelFinder(source_data, model, accession_choice = accession, return_params = model_statistic)
      
      # Set the cell value
      heatmap_table[accession, abbr] <- value
    }
  }
  
  # Account for strange bug where empty columns are added
  heatmap_table %<>% select(all_of(abbreviations))
  
  return(heatmap_table)
}

# This function takes a source data file, the desired models with abbreviations required for the heatmap (as a named list), 
#     and the target statistic, and returns a heatmap plot.
#     NOTE: This function is designed to work with the output of GENERATOR.heatmap_table, so be sure to have a table of model abbreviations by accesions as input!
GENERATOR.heatmap_plot <- function(heatmap_table, model_statistic, heatmap_title, plot_type = "Absolute") {
  
  # Get the model abbreviations as a factor ordering for ease of display
  model_abbr <- heatmap_table %>% names()
  factor_order_model_abbr <- model_abbr %>% factor(levels = model_abbr)
  
  # What happens next depends on the type of plot we use.
  # If the measured statistic is an absolute in a range, e.g. -1 to 0 to 1, 
  #     then we should use the "Absolute" display option.
  
  if (plot_type == "Absolute") {
    
    return_plot <- heatmap_table %>% 
      rownames_to_column("accession") %>%
      pivot_longer(cols = -1, names_to = "abbr_model", values_to = "statistic") %>%
      mutate(abbr_model = factor(abbr_model, levels = factor_order_model_abbr)) %>%
      ggplot(aes(abbr_model, as.factor(accession), fill = statistic)) +
      geom_tile() +
      xlab("Models") +
      ylab("Accession") +
      ggtitle(heatmap_title) +
      theme(axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
      scale_fill_viridis(direction = -1, begin = max(heatmap_table[-1]), end = ifelse(min(heatmap_table[-1]) < 0, 0, min(heatmap_table[-1]))) +
      labs(fill = model_statistic)
    plot(return_plot)
    
    # However, if the measure statistic is only meaningful relative to other 
    #     values derived in other models, then the "Relative to Model" option
    #     is more appropriate.
  } else if (plot_type == "Relative to Model") {
    
    # Take the heatmap table,
    return_plot <- heatmap_table %>%
      # Normalize over all rows (NOTE: This transposes the dataframe)
      apply(1, function(x)(x-min(x))/(max(x)-min(x))) %>%
      as.data.frame() %>%
      rownames_to_column("abbr_model") %>%
      pivot_longer(cols = -1, names_to = "accession", values_to = "statistic") %>%
      mutate(abbr_model = factor(abbr_model, levels = factor_order_model_abbr)) %>%
      ggplot(aes(abbr_model, as.factor(accession), fill = statistic)) +
      geom_tile() +
      xlab("Models") +
      ylab("Accession") +
      ggtitle(heatmap_title) +
      theme(axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
      scale_fill_viridis(256, direction = -1, begin = 0, end = 1, option = "A") +
      labs(fill = model_statistic)
    plot(return_plot)
  }
  
}

# This function generates a difference-in-means plot for a given measure over all accessions
GENERATOR.tukeydiff_plot <- function(measure, source_data = DATA.unfiltered) {
  # Start with a linear model of the association of the accession with the measure under scrutiny
  lm(source_data[[measure]] ~ factor(source_data$accession)) %>%
    # Use ANOVA to get the residuals
    aov %>%
    # Get the differences in means and the statistical significance of those differences
    TukeyHSD(., 'factor(source_data$accession)', conf.level=0.95) %>%
    # Plot the differences
    plot(las = 1, col = "black", cex.axis = 0.2)
}

# This function generates a hypergraph of differences in means for a given measure over all accessions
GENERATOR.tukeygroup_plot <- function(measure, source_data = DATA.unfiltered) {
  # Start with a linear model of the association of the accession with the measure under scrutiny
  lm(source_data[[measure]] ~ factor(source_data$accession)) %>% 
    # Use ANOVA to get the residuals
    aov %>%
    # Get the differences in means and the statistical significance of those differences
    TukeyHSD(., 'factor(source_data$accession)', conf.level=0.95) %>%
    # Unlist the data
    .[["factor(source_data$accession)"]] %>%
    # Get the "p adj" column
    .[,4] %>%
    # Create list of group membership by accession
    multcompLetters %>% 
    # Split into rectangular matrix of boolean values
    .$LetterMatrix %>%
    # Convert to dataframe
    data.frame %>% 
    # Pull out the rownames into an "accession" column
    rownames_to_column("accession") %>%
    # Sort accessions from lowest to highest
    arrange(accession) %>%
    # Create long table of boolean values
    pivot_longer(c(-1), names_to = "group") %>%
    # Filter only for rows that state group membership of an accession
    filter(value) %>%
    # Remove the boolean values
    select(-value) %>%
    # Split dataframe by group into a list of dataframes
    split(., .$group) %>%
    # Convert list of dataframes into a list of vectors, thereby creating an edgelist
    lapply(pull, accession) %>%
    # Create hypergraph from edgelist
    hypergraph_from_edgelist() %>%
    # Plot hypergraph
    plot()
  text(0, 1.05, pos = 4, cex = 1.2, paste0("Hypergraph of ", measure))
}

# This function generates a scatterplot with a fitted line from a desired model and over a given accession
GENERATOR.fitline_plot <- function(model_string, accession_choice, source_data = DATA.unfiltered, source_data_name = "Unfiltered Data") {
  
  # Acquire linear model object
  model_to_plot <- GENERATOR.modelFinder(source_data, model_string, accession_choice, return_params = c("Adj. R^2", "P-value", "Coefficients", "Model"))
  
  # Get model data for ggplot call
  model_data <- model_to_plot$Model 
  
  # Get individual variables for the plot
  y_data <- model_data$Y_data
  x_data <- model_data$X_data
  
  # Get names of dependent and independent variables
  y_name <- model_data %>% names() %>% .[[1]]
  x_name <- model_data %>% names() %>% .[[2]]
  
  # Get regression line information
  intercept <- model_to_plot$Coefficients[1]
  slope     <- model_to_plot$Coefficients[2]
  
  # Get adjusted R^2 value of the plot
  adj_rsq <- model_to_plot$`Adj. R^2`
  
  # Create the plot
  intermeasure_plot_out <- ggplot(model_data, aes(!!as.symbol(x_name), !!as.symbol(y_name))) +
    geom_point() +
    xlab(x_data) +
    ylab(y_data) +
    labs(title = (paste0("Correlation between\n", y_name, " and ", x_name, "\nACCESSION:", accession_choice, "\nADJ.R^2 = ", adj_rsq, "\n", source_data_name))) +
    # geom_text(hjust=1.05, vjust=0) +
    # geom_text(hjust=-.20, vjust=0) +
    geom_abline(intercept = intercept, slope = slope)
  
  plot(intermeasure_plot_out)
}

# This function generates a 3D scatterplot over three measures of choice from the given data.
# Can be run independently of program execution.
# To obtain the plots shown in the paper, use ("height", "width", "thickness") and ("thickness", "diameter", "fresh_weight") as parameters.
GENERATOR.3Dscatter_plot <- function (x_axis, y_axis, z_axis, source_data = DATA.unfiltered) {
  
  accessions <- source_data$accession %>% unique %>% sort
  
  plotpalette_3D <- c("#009a00", # Green
                      "#ef1101", # Bright Red
                      "#5c0300", # Burgundy
                      "#fd9206", # Pale Orange
                      "#a45d00", # Light Brown
                      "#fef636", # Pale Yellow
                      "#8f8a00", # Dirty Yellow
                      "#0006a4", # Dark Blue
                      "#0efd0e", # Bright Green
                      "#010afb", # Bright Blue
                      "#ba20fd", # Lavender
                      "#7600a8", # Purple
                      "#080f0f", # Black
                      "#faffff"  # White
  ) 
  
  plot_ly(x = x_axis, y = y_axis, z = z_axis, color = as.factor(accessions), colors =  plotpalette_3D)
  
}







### DEFINITIONS AND CALLS BOUNDARY -------------------------
# 

# All code after this section and before the "Program End" section is custom code relating specifically to the CAI analysis.
# This code is not guaranteed to work with other datasets that do not share the same initial measures.
# To create a custom workflow with new data, use the functions above as a guide, and the function calls below as examples.
###

DATA_SELECTION <- commandArgs(defaults = list(datafile = "data_files\\full_data.csv"))[1]

# Data Integrity Check and Non-Heatmap Plots ----------------

# Default/Custom Dataset Input



while(TRUE) {
  
  if (DATA_SELECTION == "data_files\\full_data.csv") {
    break
  }
  
  default_file_decision <- readline("MESSAGE: Would you like to use the default dataset? y/N: ")
  
  if (default_file_decision == "y") {
    
    DATA.unfiltered <- read.csv("data_files\\full_data.csv")
    cat("Using default dataset \"full_data.csv\"\n")
    
  } else if (default_file_decision == "N") {
    
    data_file_name <- readline("Please enter the name of your data file as shown in the \"data_files\" folder: ") 
    DATA.unfiltered <- data_file_name %>%
      paste0("data_files\\", .) %>% 
      read.csv(file = .)
    
    cat(paste0("Using dataset ", data_file_name, "\n"))
    
  } else {
    
    print("Invalid Entry - Please Try Again")
    
    next
  }
  break
}

# Begin time checking
STARTTIME <- Sys.time()

# Data Integrity Check

columns_for_analysis <- read.csv("parameters/required_columns.txt") %>%
  unlist() %>%
  as.vector()
cat("MESSAGE: Column check: ")

UTILITY.columncheck(columns_for_analysis, DATA.unfiltered, TRUE)

# Produce Filtered Dataset on 3P/N Criterion using UTILITY.datapointFlagger and GENERATOR.leverageFlagger

#   NOTE: UTILITY.datapointFlagger is run with default parameters 
#   max_flags = 1 
#   leverage_constant = 3
#   source_data = DATA.unfiltered
#   file_models_for_filtering = "leverage_filter_models.txt"

cat("\n")
data_blacklist <- UTILITY.datapointFlagger()

DATA.filtered <- DATA.unfiltered %>%
  filter(!(completeID %in% data_blacklist))

# File and Table Generation

#   Generate QQPlots on Unfiltered Data
DATA.qqplot_measures <- scan("parameters\\initial_measures.txt", character())
DATA.possible_accessions <- DATA.unfiltered$accession %>% unique %>% sort

cat("Generating QQ Plots on Unfiltered Data...\n")
for (accession in DATA.possible_accessions) {
  for (measure in DATA.qqplot_measures) {
    GENERATOR.qqGen(measure, accession)
  }
}

# Generate Violin Charts on Filtered and Unfiltered Data
DATA.initial_measures <- scan("parameters\\initial_measures.txt", character())

cat("Generating Violin Plots on Unfiltered Data...\n")
for (measure in DATA.initial_measures) {
  GENERATOR.violin_plot(DATA.unfiltered, measure, "UNF")
  GENERATOR.violin_plot(DATA.filtered, measure, "FIL")
}


# Generate Tukey plots on Filtered and Unfiltered Data
cat("Generating Tukey Plots and Hypergraphs...\n")

for (measure in DATA.initial_measures) {
  
  tukey_png_dimensions = c(600, 6, 6)
  
  # Unfiltered Data
  #   Difference in Means Plot
  UTILITY.quickPNG((paste0("Means Comparison-", measure, "--UNF")), "image_output\\tukeyplots\\mean_differences")
  GENERATOR.tukeydiff_plot(measure, source_data = DATA.unfiltered)
  dev.off()
  
  #   Hypergraph
  UTILITY.quickPNG((paste0("Hypergraph-", measure, "--UNF")), "image_output\\tukeyplots\\hypergraphs", dimensions = tukey_png_dimensions)
  GENERATOR.tukeygroup_plot(measure, source_data = DATA.unfiltered)
  dev.off()
  
  # Filtered Data
  #   Difference in Means Plot
  UTILITY.quickPNG((paste0("Means Comparison-", measure, "--FIL")), "image_output\\tukeyplots\\mean_differences")
  GENERATOR.tukeydiff_plot(measure, source_data = DATA.filtered)
  dev.off()
  
  #   Hypergraph
  UTILITY.quickPNG((paste0("Hypergraph-", measure, "--FIL")), "image_output\\tukeyplots\\hypergraphs", dimensions = tukey_png_dimensions)
  GENERATOR.tukeygroup_plot(measure, source_data = DATA.filtered)
  dev.off()
  
}

DATA.intermeasure_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\intermeasure_models.csv")

cat("Generating Intermeasure Correlations...\n")
for (accession in DATA.possible_accessions) {
  for (plot_model in DATA.intermeasure_models_and_abbrs$model) {
    
    model_abbreviation <- DATA.intermeasure_models_and_abbrs %>%
      filter(model == plot_model) %>%
      .$abbreviation
    
    UTILITY.quickPNG(paste0("Intermeasure-", accession, "-", model_abbreviation, "--UNF"), "image_output\\intermeasure\\unfiltered")
    GENERATOR.fitline_plot(plot_model, accession, source_data = DATA.unfiltered, source_data_name = "Unfiltered Data")
    dev.off()
    
    UTILITY.quickPNG(paste0("Intermeasure-", accession, "-", model_abbreviation, "--FIL"), "image_output\\intermeasure\\filtered")
    GENERATOR.fitline_plot(plot_model, accession, source_data = DATA.filtered,   source_data_name = "Filtered Data")
    dev.off()
  }
}

# Intermeasure Heatmaps---------------------

#   Generate intermeasure correlations on Unfiltered and Filtered data

#   Generate intermeasure correlations - adjusted R^2 table and heatmap
cat("Generating Intermeasure Heatmaps...\n")

UTILITY.quickPNG((paste0("Intermeasure Adj-R2--UNF")), "image_output\\heatmaps\\intermeasure")
GENERATOR.heatmap_table(DATA.intermeasure_models_and_abbrs, "Adj. R^2", source_data = DATA.unfiltered) %T>%
  write.csv("table_output\\intermeasure\\Intermeasure Adj-R2.csv") %>%
  GENERATOR.heatmap_plot("Adj. R^2", "Intermeasure Adj. R^2", plot_type = "Absolute")
dev.off()

# Main Models Heatmaps - BOX, FITTING, and ELLIPSE - FULL ---------------
#     Load the relevant model sets to be displayed in the heatmap
DATA.BOX_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\box_models.csv")
DATA.FIT_ELL_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\fitted_and_ellipse_models.csv")

#     Generate versions of these model sets containing only the terms themselves (A+B+C)
#       or their interaction terms (A*B*C)
DATA.BOX.no_interactions <- DATA.BOX_models_and_abbrs %>%
  filter(!(grepl("[*]", model)))

DATA.BOX.interactions <- DATA.BOX_models_and_abbrs %>%
  filter(!(grepl("[+]", model)))

DATA.FITTING_and_ELLIPSE.no_interactions <- DATA.FIT_ELL_models_and_abbrs %>%
  filter(!(grepl("[*]", model)))

DATA.FITTING_and_ELLIPSE.interactions <- DATA.FIT_ELL_models_and_abbrs %>%
  filter(!(grepl("[+]", model)))


#     Generate the actual heatmaps
UTILITY.heatmap_helper <- function(png_title, png_folder, models_abbr_object, statistic, source_data, csv_filepath, heatmap_title, heatmap_plottype) {
  UTILITY.quickPNG(png_title, png_folder)
  GENERATOR.heatmap_table(models_abbr_object, statistic, source_data) %T>%
    write.csv(csv_filepath) %>%
    GENERATOR.heatmap_plot(statistic, heatmap_title, plot_type = heatmap_plottype)
  dev.off()
}

#   Generate Box model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Box Model Heatmaps...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--UNF", "image_output\\heatmaps\\BOX", DATA.BOX_models_and_abbrs, "Adj. R^2", DATA.unfiltered, "table_output\\BOX\\Box Model Adj-R2 - UNF.csv", "Box Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--UNF", "image_output\\heatmaps\\BOX", DATA.BOX_models_and_abbrs, "SBC", DATA.unfiltered, "table_output\\BOX\\Box Model SBC - UNF.csv", "Box Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--FIL", "image_output\\heatmaps\\BOX", DATA.BOX_models_and_abbrs, "Adj. R^2", DATA.filtered, "table_output\\BOX\\Box Model Adj-R2 - FIL.csv", "Box Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--FIL", "image_output\\heatmaps\\BOX", DATA.BOX_models_and_abbrs, "SBC", DATA.filtered, "table_output\\BOX\\Box Model SBC - FIL.csv", "Box Model SBC - FIL", "Relative to Model")

#   Generate Fitting and Ellipse model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Fitting and Ellipse Heatmaps...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--UNF", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FIT_ELL_models_and_abbrs, "Adj. R^2", DATA.unfiltered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model Adj-R2 - UNF.csv", "Fit/Ell Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--UNF", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FIT_ELL_models_and_abbrs, "SBC", DATA.unfiltered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model SBC - UNF.csv", "Fit/Ell Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--FIL", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FIT_ELL_models_and_abbrs, "Adj. R^2", DATA.filtered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model Adj-R2 - FIL.csv", "Fit/Ell Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--FIL", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FIT_ELL_models_and_abbrs, "SBC", DATA.filtered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model SBC - FIL.csv", "Fit/Ell Model SBC - FIL", "Relative to Model")

# Main Models Heatmaps - BOX, FITTING, and ELLIPSE - NO INTERACTIONS ---------------------

#   Generate Box model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Box Model Heatmaps with No Interactions...\n")

DATA.BOX.no_interactions <- DATA.BOX_models_and_abbrs %>%
  filter(!(grepl("[*]", model)))

DATA.FITTING_and_ELLIPSE.no_interactions <- DATA.FIT_ELL_models_and_abbrs %>%
  filter(!(grepl("[*]", model)))


#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--UNF [NO Interactions]", "image_output\\heatmaps\\BOX", DATA.BOX.no_interactions, "Adj. R^2", DATA.unfiltered, "table_output\\BOX\\Box Model Adj-R2 - UNF [NO Interactions].csv", "Box Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--UNF [NO Interactions]", "image_output\\heatmaps\\BOX", DATA.BOX.no_interactions, "SBC", DATA.unfiltered, "table_output\\BOX\\Box Model SBC - UNF [NO Interactions].csv", "Box Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj.R^2--FIL [NO Interactions]", "image_output\\heatmaps\\BOX", DATA.BOX.no_interactions, "Adj. R^2", DATA.filtered, "table_output\\BOX\\Box Model Adj-R2 - FIL [NO Interactions].csv", "Box Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--FIL [NO Interactions]", "image_output\\heatmaps\\BOX", DATA.BOX.no_interactions, "SBC", DATA.filtered, "table_output\\BOX\\Box Model SBC - FIL [NO Interactions].csv", "Box Model SBC - FIL", "Relative to Model")

#   Generate Fitting and Ellipse model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Fitting and Ellipse Heatmaps with No Interactions...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--UNF [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.no_interactions, "Adj. R^2", DATA.unfiltered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model Adj-R2 - UNF [NO Interactions].csv", "Fit/Ell Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--UNF [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.no_interactions, "SBC", DATA.unfiltered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model SBC - UNF [NO Interactions].csv", "Fit/Ell Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--FIL [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.no_interactions, "Adj. R^2", DATA.filtered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model Adj-R2 - FIL [NO Interactions].csv", "Fit/Ell Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--FIL [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.no_interactions, "SBC", DATA.filtered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model SBC - FIL [NO Interactions].csv", "Fit/Ell Model SBC - FIL", "Relative to Model")


# Main Models Heatmaps - BOX, FITTING, and ELLIPSE - INTERACTIONS ---------------------
#   Generate Box model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Box Model Heatmaps WIth Interactions...\n")

DATA.BOX.interactions <- DATA.BOX_models_and_abbrs %>%
  filter(!(grepl("[+]", model)))

DATA.FITTING_and_ELLIPSE.interactions <- DATA.FIT_ELL_models_and_abbrs %>%
  filter(!(grepl("[+]", model)))

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--UNF [INTERACTIONS]", "image_output\\heatmaps\\BOX", DATA.BOX.interactions, "Adj. R^2", DATA.unfiltered, "table_output\\BOX\\Box Model Adj-R2 - UNF [INTERACTIONS].csv", "Box Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--UNF [INTERACTIONS]", "image_output\\heatmaps\\BOX", DATA.BOX.interactions, "SBC", DATA.unfiltered, "table_output\\BOX\\Box Model SBC - UNF [INTERACTIONS].csv", "Box Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--FIL [INTERACTIONS]", "image_output\\heatmaps\\BOX", DATA.BOX.interactions, "Adj. R^2", DATA.filtered, "table_output\\BOX\\Box Model Adj-R2 - FIL [INTERACTIONS].csv", "Box Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--FIL [INTERACTIONS]", "image_output\\heatmaps\\BOX", DATA.BOX.interactions, "SBC", DATA.filtered, "table_output\\BOX\\Box Model SBC - FIL [INTERACTIONS].csv", "Box Model SBC - FIL", "Relative to Model")

cat("Generating Fitting and Ellipse Heatmaps With Interactions...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--UNF [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.interactions, "Adj. R^2", DATA.unfiltered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model Adj-R2 - UNF [INTERACTIONS].csv", "Fit/Ell Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--UNF [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.interactions, "SBC", DATA.unfiltered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model SBC - UNF [INTERACTIONS].csv", "Fit/Ell Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--FIL [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.interactions, "Adj. R^2", DATA.filtered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model Adj-R2 - FIL [INTERACTIONS].csv", "Fit/Ell Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--FIL [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE", DATA.FITTING_and_ELLIPSE.interactions, "SBC", DATA.filtered, "table_output\\FITTING and ELLIPSE\\Fit-Ell Model SBC - FIL [INTERACTIONS].csv", "Fit/Ell Model SBC - FIL", "Relative to Model")

# Program end ------------------------

ENDTIME <- Sys.time()
TOTALTIME <- ENDTIME - STARTTIME

cat("Total Time Elapsed: ")
cat(TOTALTIME)
cat(" Minutes")