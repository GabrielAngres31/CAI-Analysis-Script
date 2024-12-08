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
                        "plotly", 
                        "htmlwidgets", 
                        "magrittr", 
                        "tidyr", 
                        "multiApply", 
                        "stringr", 
                        "viridis", 
                        "RColorBrewer", 
                        "data.table",
                        "here",
                        "minpack.lm")

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
UTILITY.folderpathcheck <- function(subfolder_string, silenced = TRUE) {
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
  # If a file isn't saving properly, set ` silenced = FALSE ` in the UTILITY.folderpathcheck call
  UTILITY.folderpathcheck(subfolder, silenced = TRUE)
  filepath <- file.path(subfolder, filenamestring, fsep = "\\")
  
  # Create the PNG
  png(paste0(filepath, ".png"), 
      width = png_width * png_ppi, 
      height = png_height * png_ppi, 
      res = png_ppi)
}

# Creates a CSV at a given file path.
# If a file isn't saving properly, set ` silenced = FALSE ` in the UTILITY.folderpathcheck call
UTILITY.quickCSV <- function(csv_object, filenamestring, subfolder) {
  UTILITY.folderpathcheck(subfolder, silenced = TRUE)
  write.csv(csv_object, paste0(subfolder, "\\", filenamestring))
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
    
    if(min(heatmap_table[-1]) < 0) {
      cat("[")
      cat(heatmap_title)
      cat("]'s minimum value is lower than 0\n")
      
    }
    
    #color_vector <- c( mako(length(heatmap_table[-1])), rep("grey80", times = length(heatmap_table[-1])) )
    return_plot <- heatmap_table %>% 
      rownames_to_column("accession") %>%
      pivot_longer(cols = -1, names_to = "abbr_model", values_to = "statistic") %>%
      mutate(abbr_model = factor(abbr_model, levels = factor_order_model_abbr)) %>%
      ggplot(aes(abbr_model, as.factor(accession), fill = statistic)) +
      geom_tile() +
      xlab("Models") +
      ylab("Accession") +
      theme(plot.title = element_text(size = 10)) +
      ggtitle(heatmap_title) +
      theme(axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
      # scale_color_gradientn(colors = color_vector,
      #                       values = c(0, 1), na.value = "pink",
      #                       breaks = seq(0, 1, 0.005),limits = c(0, 1)) +
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
      theme(plot.title = element_text(size = 10)) +
      ggtitle(heatmap_title) +
      theme(axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
      scale_fill_viridis(256, direction = -1, begin = 0, end = 1, option = "A") +
      labs(fill = model_statistic)
    plot(return_plot)
  }
  
}

# This function calls the GENERATOR functions heatmap_table and heatmap_plot in sequence, allowing a table and its corresponding heatmap to be generated in the same call.
# A reference parameter can be added in the form c(variable_name, value) to generate a single-value column in the heatmap as a reference.
# If confirmation_dialogue is TRUE, the console will ask whether you want to generate the heatmap with the given title (with readline). Any input other than "Y" will abort the call.

UTILITY.heatmap_helper <- function(png_title, png_folder, 
                                   models_abbr_object, statistic, 
                                   source_data, csv_title, csv_filepath, heatmap_title, 
                                   heatmap_plottype, 
                                   reference = NULL, confirmation_dialogue = FALSE) {
  
  if(confirmation_dialogue) {
    decision <- readline(paste0("Generate: ", heatmap_title, " ? | "))
    if(!(decision == "Y")) {
      return(NULL)
    }
  }
 
  
  UTILITY.quickPNG(png_title, png_folder)
  
  GENERATOR.heatmap_table(models_abbr_object, statistic, source_data) %>%
    {if(!(is.null(reference))) mutate(., !! reference[1] := as.numeric(reference[2])) else .} %T>%
    UTILITY.quickCSV(csv_title, csv_filepath) %>%
    GENERATOR.heatmap_plot(statistic, heatmap_title, plot_type = heatmap_plottype)
  
  dev.off()
}

UTILITY.heatmap_delta <- function(table_1_path, 
                                  table_2_path,
                                  png_title, png_folder, 
                                  models_abbr_object, statistic, 
                                  csv_title, csv_filepath, 
                                  heatmap_title, heatmap_plottype, 
                                  confirmation_dialogue = TRUE) {
  
  # Confirm that the delta table should be generated
  
  if(confirmation_dialogue) {
    decision <- readline(paste0("Generate delta table: ", heatmap_title, " ? | "))
    if(!(decision == "Y")) {
      return(NULL)
    }
  }
  
  # Read the tables
  
  table_1 <- read.csv(table_1_path) %>% setnames(c("accession", models_abbr_object$abbreviation))
  table_2 <- read.csv(table_2_path) %>% setnames(c("accession", models_abbr_object$abbreviation))
  
  # Assert that the tables are sufficiently identically formatted to compare
  
  stopifnot(rownames(table_1) == rownames(table_2))
  stopifnot(    nrow(table_1) ==     nrow(table_2))
  stopifnot(    ncol(table_1) ==     ncol(table_2))
  
  # Generate the deltas
  
  table_out <- table_1 
  table_out[-1] <- table_1[-1] - table_2[-1]
  
  # Generate a PNG for the heatmap to go into
  
  UTILITY.quickPNG(png_title, png_folder)
  
  # Write the .csv and .png files
  
  model_abbr <- table_out %>% names()

  factor_order_model_abbr <- model_abbr %>% factor(levels = model_abbr)
  
  if(heatmap_plottype == "Absolute") {
    return_plot <- table_out %>%
      column_to_rownames(var="accession") %T>%
      UTILITY.quickCSV(csv_title, csv_filepath) %>%
      rownames_to_column("accession") %>%
      pivot_longer(cols = -1, names_to = "abbr_model", values_to = "statistic") %>%
      mutate(abbr_model = factor(abbr_model, levels = factor_order_model_abbr)) %>%
      ggplot(aes(abbr_model, as.factor(accession), fill = statistic)) +
      geom_tile() +
      xlab("Models") +
      ylab("Accession") +
      theme(plot.title = element_text(size = 10)) +
      ggtitle(heatmap_title) +
      theme(axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
      # scale_color_gradientn(colors = color_vector,
      #                       values = c(0, 1), na.value = "pink",
      #                       breaks = seq(0, 1, 0.005),limits = c(0, 1)) +
      
      #scale_fill_viridis(direction = -1, begin = max(table_out[-1]), end = ifelse(min(table_out[-1]) < 0, 0, min(table_out[-1]))) +
      
      scale_fill_gradient2(low = "#FF0000",
                          mid = "#FFFFFF",
                          high = "#075AFF") +
      
      labs(fill = statistic)
    plot(return_plot)
    
    # However, if the measure statistic is only meaningful relative to other 
    #     values derived in other models, then the "Relative to Model" option
    #     is more appropriate.
  } else if (plot_type == "Relative to Model") {
  
  # Take the heatmap table,
    return_plot <- table_out %>%
      UTILITY.quickCSV(csv_title, csv_filepath) %>%
      ggplot(aes(abbr_model, as.factor(accession), fill = statistic)) +
      geom_tile() +
      xlab("Models") +
      ylab("Accession") +
      theme(plot.title = element_text(size = 10)) +
      ggtitle(heatmap_title) +
      theme(axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
      scale_fill_viridis(256, direction = -1, begin = 0, end = 1, option = "A") +
      labs(fill = statistic)
    plot(return_plot)
  }
    
  dev.off()
  
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
# TODO: THIS ISN'T DISPLAYING A GRAPH
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
  source_data %>%
  plot_ly(x = x_axis, y = y_axis, z = z_axis, color = as.factor(accessions), colors =  plotpalette_3D) %>%
    add_markers(marker = list(size = 10, opacity = 1)) %>%
    layout(scene = list(xaxis = list(title = "X-axis"),
                        yaxis = list(title = "Y-axis"),
                        zaxis = list(title = "Z-axis")),
           title = "3D Scatter Plot")
}







### DEFINITIONS AND CALLS BOUNDARY -------------------------
# 

# All code after this section and before the "Program End" section is custom code relating specifically to the CAI analysis.
# This code is not guaranteed to work with other datasets that do not share the same initial measures.
# To create a custom workflow with new data, use the functions above as a guide, and the function calls below as examples.
###

DATA_SELECTION <- commandArgs()[1]

# Data Integrity Check and Non-Heatmap Plots ----------------

# Default/Custom Dataset Input



while(TRUE) {
  
  if (DATA_SELECTION == "data_files\\full_data.csv") {
    break
  }
  
  default_file_decision <- readline("MESSAGE: Would you like to use the default dataset? y/N: ")
  
  if (default_file_decision == "y") {
    
    # Check to see if the full dataset for use has been generated from a prior run of the program.
    # In this case, this is configured for the CAI dataset.
    # If not, it is generated from existing, complementary datasets.
    # For your own analyses, place your own code in here.
    # If one or more datasets are missing, or if there are no datasets, the program will stop.
    
    if (!file.exists("full_data.csv")) {
      
      columns_for_analysis <- read.csv("parameters\\required_columns_raw_in.txt") %>%
        unlist() %>%
        as.vector()
      
      read.csv("data_files\\raw_data.csv") %T>% 
        
        UTILITY.columncheck(columns_for_analysis, ., TRUE) %>%
        
        mutate(H_div_W = height/width,
               FW_div_W = fresh_weight/width,
               FW_div_D = fresh_weight/diameter,
               FW_div_H = fresh_weight/height,
               FW_div_T = fresh_weight/thickness,
               D_div_W = diameter/width,
               Theo_Area = pi*height*width*0.25,
               PartRatio = ((height-width)/(height+width))^2,
               Pade_Peri = pi*(height+width)*(64-3*PartRatio^2)/(64-16*PartRatio),
               Pade_Derived_Diam = Pade_Peri/pi) %>%
        
        merge(read.csv("data_files\\pad_area_estimations.csv"), by = "completeID") %>%
        
        merge(read.csv("data_files\\water_content_data.csv"  ), by = "completeID") %>%
        
        mutate(dry_weight = fresh_weight*dry_proportion,
               water_weight = fresh_weight*water_proportion) %>%
        
        write.csv("data_files\\full_data.csv")
    }
    
    # Load initial program dataset ("Full" dataset) with column check for data integrity
    
    DATA.unfiltered <- read.csv("data_files\\full_data.csv")
    
    read.csv("parameters/required_columns_full_data.txt") %>%
      unlist() %>%
      as.vector() %>%
      UTILITY.columncheck(DATA.unfiltered, TRUE)
    
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

cat("Generating Violin Plots on Unfiltered/Filtered Data...\n")
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
DATA.intermeasure_with_dry_weight_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\intermeasure_models_with_dry_weight.csv")

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
  
  for (plot_model in DATA.intermeasure_with_dry_weight_models_and_abbrs$model) {
    
    model_abbreviation <- DATA.intermeasure_with_dry_weight_models_and_abbrs %>%
      filter(model == plot_model) %>%
      .$abbreviation
    
    UTILITY.quickPNG(paste0("Intermeasure-", accession, "-", model_abbreviation, "--UNF"), "image_output\\intermeasure\\unfiltered_DW")
    GENERATOR.fitline_plot(plot_model, accession, source_data = DATA.unfiltered, source_data_name = "Unfiltered Data")
    dev.off()
    
    UTILITY.quickPNG(paste0("Intermeasure-", accession, "-", model_abbreviation, "--FIL"), "image_output\\intermeasure\\filtered_DW")
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
  UTILITY.quickCSV("Intermeasure Adj-R2.csv", "table_output\\intermeasure") %>%
  GENERATOR.heatmap_plot("Adj. R^2", "Intermeasure Adj. R^2", plot_type = "Absolute")
dev.off()

UTILITY.quickPNG((paste0("Intermeasure Adj-R2 DW--UNF")), "image_output\\heatmaps\\intermeasure")
GENERATOR.heatmap_table(DATA.intermeasure_with_dry_weight_models_and_abbrs, "Adj. R^2", source_data = DATA.unfiltered) %T>%
  UTILITY.quickCSV("Intermeasure Adj-R2 DW.csv", "table_output\\intermeasure") %>%
  GENERATOR.heatmap_plot("Adj. R^2", "Intermeasure Adj. R^2 DW", plot_type = "Absolute")
dev.off()

# Main Models Heatmaps - BOX, FITTING, and ELLIPSE - FULL ---------------
#     Load the relevant model sets to be displayed in the heatmap

DATA.BOX_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\box_models.csv")
DATA.BOX_DRY_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\box_models_using_dry_weight.csv")

DATA.FIT_ELL_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\fitted_and_ellipse_models.csv")
DATA.FIT_ELL_DRY_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\fitted_and_ellipse_models_using_dry_weight.csv")

#     Generate versions of these model sets containing only the terms themselves (A+B+C)
#       or their interaction terms (A*B*C)

DATA.BOX_models_and_abbrs %T>%
  {assign("DATA.BOX.interactions",    value = filter(., !(grepl("[+]", model))), envir = .GlobalEnv)} %T>%
  {assign("DATA.BOX.no_interactions", value = filter(., !(grepl("[*]", model))), envir = .GlobalEnv)}

DATA.BOX_DRY_models_and_abbrs %T>%
  {assign("DATA.BOX_DRY.interactions",    value = filter(., !(grepl("[+]", model))), envir = .GlobalEnv)} %T>%
  {assign("DATA.BOX_DRY.no_interactions", value = filter(., !(grepl("[*]", model))), envir = .GlobalEnv)}

DATA.FIT_ELL_models_and_abbrs %T>%
  {assign("DATA.FITTING_and_ELLIPSE.interactions",    value = filter(., !(grepl("[+]", model))), envir = .GlobalEnv)} %T>%
  {assign("DATA.FITTING_and_ELLIPSE.no_interactions", value = filter(., !(grepl("[*]", model))), envir = .GlobalEnv)}

DATA.FIT_ELL_DRY_models_and_abbrs %T>%
  {assign("DATA.FITTING_and_ELLIPSE_DRY.interactions",    value = filter(., !(grepl("[+]", model))), envir = .GlobalEnv)} %T>%
  {assign("DATA.FITTING_and_ELLIPSE_DRY.no_interactions", value = filter(., !(grepl("[*]", model))), envir = .GlobalEnv)}

#     Generate the actual heatmaps

#   Generate Box model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Box Model Heatmaps...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--UNF", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX_models_and_abbrs, "Adj. R^2", DATA.unfiltered, "Box Model Adj-R2 - UNF.csv", "table_output\\BOX\\Fresh Weight", "Box Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--UNF", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX_models_and_abbrs, "SBC", DATA.unfiltered, "Box Model SBC - UNF.csv", "table_output\\BOX\\Fresh Weight", "Box Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--FIL", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX_models_and_abbrs, "Adj. R^2", DATA.filtered, "Box Model Adj-R2 - FIL.csv", "table_output\\BOX\\Fresh Weight", "Box Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--FIL", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX_models_and_abbrs, "SBC", DATA.filtered, "Box Model SBC - FIL.csv", "table_output\\BOX\\Fresh Weight", "Box Model SBC - FIL", "Relative to Model")

#--

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2 DRY--UNF", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY_models_and_abbrs, "Adj. R^2", DATA.unfiltered, "Box Model Adj-R2 DRY - UNF.csv", "table_output\\BOX\\Dry Weight", "Box Model Adj. R^2 DRY - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC DRY--UNF", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY_models_and_abbrs, "SBC", DATA.unfiltered, "Box Model SBC DRY - UNF.csv", "table_output\\BOX\\Dry Weight", "Box Model SBC DRY - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2 DRY--FIL", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY_models_and_abbrs, "Adj. R^2", DATA.filtered, "Box Model Adj-R2 DRY - FIL.csv", "table_output\\BOX\\Dry Weight", "Box Model Adj. R^2 DRY - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC DRY--FIL", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY_models_and_abbrs, "SBC", DATA.filtered, "Box Model SBC DRY - FIL.csv", "table_output\\BOX\\Dry Weight", "Box Model SBC DRY - FIL", "Relative to Model")

#--

#     Unfiltered Adj. R^2 Box Model Deltas
UTILITY.heatmap_delta("table_output\\BOX\\Dry Weight\\Box Model Adj-R2 DRY - UNF.csv",
                      "table_output\\BOX\\Fresh Weight\\Box Model Adj-R2 - UNF.csv",
                      "Box Model Adj-R2 UNF Fresh-Dry Deltas", "image_output\\heatmaps\\BOX\\Deltas", 
                      DATA.BOX_models_and_abbrs, "Adj. R^2", 
                      "Box Model Adj-R2 UNF Deltas.csv", "table_output\\BOX\\Deltas",  
                      "Box Model Adj. R^2 UNF Deltas", "Absolute", 
                      confirmation_dialogue = FALSE)

#     Filtered Adj. R^2 Box Model Deltas
UTILITY.heatmap_delta("table_output\\BOX\\Dry Weight\\Box Model Adj-R2 DRY - FIL.csv",
                      "table_output\\BOX\\Fresh Weight\\Box Model Adj-R2 - FIL.csv",
                      "Box Model Adj-R2 FIL Fresh-Dry Deltas", "image_output\\heatmaps\\BOX\\Deltas", 
                      DATA.BOX_models_and_abbrs, "Adj. R^2", 
                      "Box Model Adj-R2 FIL Deltas.csv", "table_output\\BOX\\Deltas",  
                      "Box Model Adj. R^2 FIL Deltas", "Absolute", 
                      confirmation_dialogue = FALSE)


###### ---
###### ---
###### ---

#   Generate Fitting and Ellipse model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Fitting and Ellipse Heatmaps...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--UNF", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FIT_ELL_models_and_abbrs, "Adj. R^2", DATA.unfiltered, "Fit-Ell Model Adj-R2 - UNF.csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--UNF", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FIT_ELL_models_and_abbrs, "SBC", DATA.unfiltered, "Fit-Ell Model SBC - UNF.csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--FIL", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FIT_ELL_models_and_abbrs, "Adj. R^2", DATA.filtered, "Fit-Ell Model Adj-R2 - FIL.csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--FIL", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FIT_ELL_models_and_abbrs, "SBC", DATA.filtered, "Fit-Ell Model SBC - FIL.csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model SBC - FIL", "Relative to Model")

#---

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2 DRY--UNF", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FIT_ELL_DRY_models_and_abbrs, "Adj. R^2", DATA.unfiltered, "Fit-Ell Model Adj-R2 DRY - UNF.csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model Adj. R^2 DRY - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC DRY--UNF", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FIT_ELL_DRY_models_and_abbrs, "SBC", DATA.unfiltered, "Fit-Ell Model SBC DRY - UNF.csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model SBC DRY - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2 DRY--FIL", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FIT_ELL_DRY_models_and_abbrs, "Adj. R^2", DATA.filtered, "Fit-Ell Model Adj-R2 DRY - FIL.csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model Adj. R^2 DRY - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC DRY--FIL", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FIT_ELL_DRY_models_and_abbrs, "SBC", DATA.filtered, "Fit-Ell Model SBC DRY - FIL.csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model SBC DRY - FIL", "Relative to Model")

#---

#   Unfiltered Adj. R^2 Deltas
UTILITY.heatmap_delta("table_output\\FITTING AND ELLIPSE\\Dry Weight\\Fit-Ell Model Adj-R2 DRY - UNF.csv",
                      "table_output\\FITTING AND ELLIPSE\\Fresh Weight\\Fit-Ell Model Adj-R2 - UNF.csv",
                      "Fit-Ell Model Adj-R2 UNF Fresh-Dry Deltas", "image_output\\heatmaps\\FITTING AND ELLIPSE\\Deltas",
                      DATA.FIT_ELL_models_and_abbrs, "Adj. R^2",
                      "Fit-Ell Model Adj-R2 UNF Deltas.csv", "table_output\\FITTING AND ELLIPSE\\Deltas",
                      "Fit-Ell Model Adj. R^2 UNF Deltas", "Absolute",
                      confirmation_dialogue = FALSE)

#   Filtered Adj. R^2 Deltas
UTILITY.heatmap_delta("table_output\\FITTING AND ELLIPSE\\Dry Weight\\Fit-Ell Model Adj-R2 DRY - FIL.csv",
                      "table_output\\FITTING AND ELLIPSE\\Fresh Weight\\Fit-Ell Model Adj-R2 - FIL.csv",
                      "Fit-Ell Model Adj-R2 FIL Fresh-Dry Deltas", "image_output\\heatmaps\\FITTING AND ELLIPSE\\Deltas",
                      DATA.FIT_ELL_models_and_abbrs, "Adj. R^2",
                      "Fit-Ell Model Adj-R2 FIL Deltas.csv", "table_output\\FITTING AND ELLIPSE\\Deltas",
                      "Fit-Ell Model Adj. R^2 FIL Deltas", "Absolute",
                      confirmation_dialogue = FALSE)



###### ---
###### ---
###### ---

# Main Models Heatmaps - BOX, FITTING, and ELLIPSE - NO INTERACTIONS ---------------------

#   Generate Box model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Box Model Heatmaps with No Interactions...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--UNF [NO Interactions]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.no_interactions, "Adj. R^2", DATA.unfiltered, "Box Model Adj-R2 - UNF [NO Interactions].csv", "table_output\\BOX\\Fresh Weight", "Box Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--UNF [NO Interactions]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.no_interactions, "SBC", DATA.unfiltered, "Box Model SBC - UNF [NO Interactions].csv", "table_output\\BOX\\Fresh Weight", "Box Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj.R^2--FIL [NO Interactions]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.no_interactions, "Adj. R^2", DATA.filtered, "Box Model Adj-R2 - FIL [NO Interactions].csv", "table_output\\BOX\\Fresh Weight", "Box Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--FIL [NO Interactions]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.no_interactions, "SBC", DATA.filtered, "Box Model SBC - FIL [NO Interactions].csv", "table_output\\BOX\\Fresh Weight", "Box Model SBC - FIL", "Relative to Model")

#--
#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2 DRY--UNF [NO Interactions]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.no_interactions, "Adj. R^2", DATA.unfiltered, "Box Model Adj-R2 DRY - UNF [NO Interactions].csv", "table_output\\BOX\\Dry Weight", "Box Model Adj. R^2 DRY - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC DRY--UNF [NO Interactions]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.no_interactions, "SBC", DATA.unfiltered, "Box Model SBC DRY - UNF [NO Interactions].csv", "table_output\\BOX\\Dry Weight", "Box Model SBC DRY - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj.R^2 DRY--FIL [NO Interactions]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.no_interactions, "Adj. R^2", DATA.filtered, "Box Model Adj-R2 DRY - FIL [NO Interactions].csv", "table_output\\BOX\\Dry Weight", "Box Model Adj. R^2 DRY - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC DRY--FIL [NO Interactions]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.no_interactions, "SBC", DATA.filtered, "Box Model SBC DRY - FIL [NO Interactions].csv", "table_output\\BOX\\Dry Weight", "Box Model SBC DRY - FIL", "Relative to Model")

#---

# Unfiltered Adj. R^2 Box Model Deltas - NO INTERACTIONS
UTILITY.heatmap_delta("table_output\\BOX\\Dry Weight\\Box Model Adj-R2 DRY - UNF [NO Interactions].csv",
                      "table_output\\BOX\\Fresh Weight\\Box Model Adj-R2 - UNF [NO Interactions].csv",
                      "Box Model Adj-R2 UNF Fresh-Dry Deltas [NO Interactions]", "image_output\\heatmaps\\BOX\\Deltas", 
                      DATA.BOX.no_interactions, "Adj. R^2", 
                      "Box Model Adj-R2 UNF Deltas [NO Interactions].csv", "table_output\\BOX\\Deltas",  
                      "Box Model Adj. R^2 UNF Deltas [NO Interactions]", "Absolute", 
                      confirmation_dialogue = FALSE)

# Filtered Adj. R^2 Box Model Deltas - NO INTERACTIONS
UTILITY.heatmap_delta("table_output\\BOX\\Dry Weight\\Box Model Adj-R2 DRY - FIL [NO Interactions].csv",
                      "table_output\\BOX\\Fresh Weight\\Box Model Adj-R2 - FIL [NO Interactions].csv",
                      "Box Model Adj-R2 FIL Fresh-Dry Deltas [NO Interactions]", "image_output\\heatmaps\\BOX\\Deltas", 
                      DATA.BOX.no_interactions, "Adj. R^2", 
                      "Box Model Adj-R2 FIL Deltas [NO Interactions].csv", "table_output\\BOX\\Deltas",  
                      "Box Model Adj. R^2 FIL Deltas [NO Interactions]", "Absolute", 
                      confirmation_dialogue = FALSE)


###### ---
###### ---
###### ---

#   Generate Fitting and Ellipse model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Fitting and Ellipse Heatmaps with No Interactions...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--UNF [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.no_interactions, "Adj. R^2", DATA.unfiltered, "Fit-Ell Model Adj-R2 - UNF [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--UNF [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.no_interactions, "SBC", DATA.unfiltered, "Fit-Ell Model SBC - UNF [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--FIL [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.no_interactions, "Adj. R^2", DATA.filtered, "Fit-Ell Model Adj-R2 - FIL [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--FIL [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.no_interactions, "SBC", DATA.filtered, "Fit-Ell Model SBC - FIL [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model SBC - FIL", "Relative to Model")

#---

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2 DRY--UNF [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.no_interactions, "Adj. R^2", DATA.unfiltered, "Fit-Ell Model Adj-R2 DRY - UNF [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model Adj. R^2 DRY - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC DRY--UNF [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.no_interactions, "SBC", DATA.unfiltered, "Fit-Ell Model SBC DRY - UNF [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model SBC DRY - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2 DRY--FIL [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.no_interactions, "Adj. R^2", DATA.filtered, "Fit-Ell Model Adj-R2 DRY - FIL [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model Adj. R^2 DRY - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC DRY--FIL [NO Interactions]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.no_interactions, "SBC", DATA.filtered, "Fit-Ell Model SBC DRY - FIL [NO Interactions].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model SBC DRY - FIL", "Relative to Model")

#---

#   Unfiltered Adj. R^2 Deltas - NO INTERACTIONS
UTILITY.heatmap_delta("table_output\\FITTING AND ELLIPSE\\Dry Weight\\Fit-Ell Model Adj-R2 DRY - UNF [NO Interactions].csv",
                      "table_output\\FITTING AND ELLIPSE\\Fresh Weight\\Fit-Ell Model Adj-R2 - UNF [NO Interactions].csv",
                      "Fit-Ell Model Adj-R2 UNF Fresh-Dry Deltas [NO Interactions]", "image_output\\heatmaps\\FITTING AND ELLIPSE\\Deltas",
                      DATA.FITTING_and_ELLIPSE.no_interactions, "Adj. R^2",
                      "Fit-Ell Model Adj-R2 UNF Deltas [NO Interactions].csv", "table_output\\FITTING AND ELLIPSE\\Deltas",
                      "Fit-Ell Model Adj. R^2 UNF Deltas [NO Interactions]", "Absolute",
                      confirmation_dialogue = FALSE)

#   Filtered Adj. R^2 Deltas - NO INTERACTIONS
UTILITY.heatmap_delta("table_output\\FITTING AND ELLIPSE\\Dry Weight\\Fit-Ell Model Adj-R2 DRY - FIL [NO Interactions].csv",
                      "table_output\\FITTING AND ELLIPSE\\Fresh Weight\\Fit-Ell Model Adj-R2 - FIL [NO Interactions].csv",
                      "Fit-Ell Model Adj-R2 FIL Fresh-Dry Deltas [NO Interactions]", "image_output\\heatmaps\\FITTING AND ELLIPSE\\Deltas",
                      DATA.FITTING_and_ELLIPSE.no_interactions, "Adj. R^2",
                      "Fit-Ell Model Adj-R2 FIL Deltas [NO Interactions].csv", "table_output\\FITTING AND ELLIPSE\\Deltas",
                      "Fit-Ell Model Adj. R^2 FIL Deltas [NO Interactions]", "Absolute",
                      confirmation_dialogue = FALSE)

###### ---
###### ---
###### ---

# Main Models Heatmaps - BOX, FITTING, and ELLIPSE - INTERACTIONS ---------------------
#   Generate Box model Heatmaps, Filtered and Unfiltered, Adj. R^2 and SBC
cat("Generating Box Model Heatmaps WIth Interactions...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--UNF [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.interactions, "Adj. R^2", DATA.unfiltered, "Box Model Adj-R2 - UNF [INTERACTIONS].csv", "table_output\\BOX\\Fresh Weight", "Box Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--UNF [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.interactions, "SBC", DATA.unfiltered, "Box Model SBC - UNF [INTERACTIONS].csv", "table_output\\BOX\\Fresh Weight", "Box Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2--FIL [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.interactions, "Adj. R^2", DATA.filtered, "Box Model Adj-R2 - FIL [INTERACTIONS].csv", "table_output\\BOX\\Fresh Weight", "Box Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC--FIL [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Fresh Weight", DATA.BOX.interactions, "SBC", DATA.filtered, "Box Model SBC - FIL [INTERACTIONS].csv", "table_output\\BOX\\Fresh Weight", "Box Model SBC - FIL", "Relative to Model")
#---

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2 DRY--UNF [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.interactions, "Adj. R^2", DATA.unfiltered, "Box Model Adj-R2 DRY - UNF [INTERACTIONS].csv", "table_output\\BOX\\Dry Weight", "Box Model Adj. R^2 DRY - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Box Model SBC DRY--UNF [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.interactions, "SBC", DATA.unfiltered, "Box Model SBC DRY - UNF [INTERACTIONS].csv", "table_output\\BOX\\Dry Weight", "Box Model SBC DRY - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Box Model Adj-R2 DRY--FIL [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.interactions, "Adj. R^2", DATA.filtered, "Box Model Adj-R2 DRY - FIL [INTERACTIONS].csv", "table_output\\BOX\\Dry Weight", "Box Model Adj. R^2 DRY - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Box Model SBC DRY--FIL [INTERACTIONS]", "image_output\\heatmaps\\BOX\\Dry Weight", DATA.BOX_DRY.interactions, "SBC", DATA.filtered, "Box Model SBC DRY - FIL [INTERACTIONS].csv", "table_output\\BOX\\Dry Weight", "Box Model SBC DRY - FIL", "Relative to Model")

#---

# Unfiltered Adj. R^2 Box Model Deltas - INTERACTIONS
UTILITY.heatmap_delta("table_output\\BOX\\Dry Weight\\Box Model Adj-R2 DRY - UNF [Interactions].csv",
                      "table_output\\BOX\\Fresh Weight\\Box Model Adj-R2 - UNF [Interactions].csv",
                      "Box Model Adj-R2 UNF Fresh-Dry Deltas [Interactions]", "image_output\\heatmaps\\BOX\\Deltas", 
                      DATA.BOX.interactions, "Adj. R^2", 
                      "Box Model Adj-R2 UNF Deltas [Interactions].csv", "table_output\\BOX\\Deltas",  
                      "Box Model Adj. R^2 UNF Deltas [Interactions]", "Absolute", 
                      confirmation_dialogue = FALSE)

# Filtered Adj. R^2 Box Model Deltas - INTERACTIONS
UTILITY.heatmap_delta("table_output\\BOX\\Dry Weight\\Box Model Adj-R2 DRY - FIL [Interactions].csv",
                      "table_output\\BOX\\Fresh Weight\\Box Model Adj-R2 - FIL [Interactions].csv",
                      "Box Model Adj-R2 FIL Fresh-Dry Deltas [Interactions]", "image_output\\heatmaps\\BOX\\Deltas", 
                      DATA.BOX.interactions, "Adj. R^2", 
                      "Box Model Adj-R2 FIL Deltas [Interactions].csv", "table_output\\BOX\\Deltas",  
                      "Box Model Adj. R^2 FIL Deltas [Interactions]", "Absolute", 
                      confirmation_dialogue = FALSE)

###### ---
###### ---
###### ---

cat("Generating Fitting and Ellipse Heatmaps With Interactions...\n")

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--UNF [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.interactions, "Adj. R^2", DATA.unfiltered, "Fit-Ell Model Adj-R2 - UNF [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model Adj. R^2 - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--UNF [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.interactions, "SBC", DATA.unfiltered, "Fit-Ell Model SBC - UNF [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model SBC - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2--FIL [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.interactions, "Adj. R^2", DATA.filtered, "Fit-Ell Model Adj-R2 - FIL [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model Adj. R^2 - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC--FIL [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Fresh Weight", DATA.FITTING_and_ELLIPSE.interactions, "SBC", DATA.filtered, "Fit-Ell Model SBC - FIL [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Fresh Weight", "Fit/Ell Model SBC - FIL", "Relative to Model")

#---

#     Unfiltered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2 DRY--UNF [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.interactions, "Adj. R^2", DATA.unfiltered, "Fit-Ell Model Adj-R2 DRY - UNF [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model Adj. R^2 DRY - UNF", "Absolute")

#     Unfiltered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC DRY--UNF [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.interactions, "SBC", DATA.unfiltered, "Fit-Ell Model SBC DRY - UNF [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model SBC DRY - UNF", "Relative to Model")

#     Filtered Adj. R^2
UTILITY.heatmap_helper(
  "Fit-Ell Model Adj-R2 DRY--FIL [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.interactions, "Adj. R^2", DATA.filtered, "Fit-Ell Model Adj-R2 DRY - FIL [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model Adj. R^2 DRY - FIL", "Absolute")

#     Filtered SBC
UTILITY.heatmap_helper(
  "Fit-Ell Model SBC DRY--FIL [INTERACTIONS]", "image_output\\heatmaps\\FITTING and ELLIPSE\\Dry Weight", DATA.FITTING_and_ELLIPSE_DRY.interactions, "SBC", DATA.filtered, "Fit-Ell Model SBC DRY - FIL [INTERACTIONS].csv", "table_output\\FITTING and ELLIPSE\\Dry Weight", "Fit/Ell Model SBC DRY - FIL", "Relative to Model")

# ---

#   Unfiltered Adj. R^2 Deltas - INTERACTIONS
UTILITY.heatmap_delta("table_output\\FITTING AND ELLIPSE\\Dry Weight\\Fit-Ell Model Adj-R2 DRY - UNF [INTERACTIONS].csv",
                      "table_output\\FITTING AND ELLIPSE\\Fresh Weight\\Fit-Ell Model Adj-R2 - UNF [INTERACTIONS].csv",
                      "Fit-Ell Model Adj-R2 UNF Fresh-Dry Deltas [INTERACTIONS]", "image_output\\heatmaps\\FITTING AND ELLIPSE\\Deltas",
                      DATA.FITTING_and_ELLIPSE.interactions, "Adj. R^2",
                      "Fit-Ell Model Adj-R2 UNF Deltas [INTERACTIONS].csv", "table_output\\FITTING AND ELLIPSE\\Deltas",
                      "Fit-Ell Model Adj. R^2 UNF Deltas [INTERACTIONS]", "Absolute",
                      confirmation_dialogue = FALSE)

#   Filtered Adj. R^2 Deltas - INTERACTIONS
UTILITY.heatmap_delta("table_output\\FITTING AND ELLIPSE\\Dry Weight\\Fit-Ell Model Adj-R2 DRY - FIL [INTERACTIONS].csv",
                      "table_output\\FITTING AND ELLIPSE\\Fresh Weight\\Fit-Ell Model Adj-R2 - FIL [INTERACTIONS].csv",
                      "Fit-Ell Model Adj-R2 FIL Fresh-Dry Deltas [INTERACTIONS]", "image_output\\heatmaps\\FITTING AND ELLIPSE\\Deltas",
                      DATA.FITTING_and_ELLIPSE.interactions, "Adj. R^2",
                      "Fit-Ell Model Adj-R2 FIL Deltas [INTERACTIONS].csv", "table_output\\FITTING AND ELLIPSE\\Deltas",
                      "Fit-Ell Model Adj. R^2 FIL Deltas [INTERACTIONS]", "Absolute",
                      confirmation_dialogue = FALSE)

###### ---
###### ---
###### ---

# Prior Work - REIS et. al. -----------------

cat("Generating REIS et. al. Comparison Heatmaps...\n")

DATA.REIS_models_and_abbrs.AREA <- read.csv("parameters\\heatmap_sources\\REIS_area_models.csv")

#     Unfiltered AREA Adj. R^2
UTILITY.heatmap_helper(
  "Reis et al AREA Adj-R2--UNF", "image_output\\heatmaps\\REIS", DATA.REIS_models_and_abbrs.AREA, "Adj. R^2", DATA.unfiltered, "REIS et al Comparisons - AREA - UNF.csv", "table_output\\REIS", "REIS et al AREA Comparison - UNF", "Absolute", c("Reis Area", 0.91))

#     Filtered AREA Adj. R^2
UTILITY.heatmap_helper(
  "Reis et al AREA Adj-R2--FIL", "image_output\\heatmaps\\REIS", DATA.REIS_models_and_abbrs.AREA, "Adj. R^2", DATA.filtered, "REIS et al Comparisons - AREA - FIL.csv", "table_output\\REIS", "REIS et al AREA Comparison - FIL", "Absolute", c("Reis Area", 0.91))


DATA.REIS_models_and_abbrs.DRYWEIGHT <- read.csv("parameters\\heatmap_sources\\REIS_dryweight_models.csv")

#     Unfiltered DRY WEIGHT Adj. R^2
UTILITY.heatmap_helper(
  "Reis et al DW Adj-R2--UNF", "image_output\\heatmaps\\REIS", DATA.REIS_models_and_abbrs.DRYWEIGHT, "Adj. R^2", DATA.unfiltered, "REIS et al Comparisons - DW - UNF.csv", "table_output\\REIS", "REIS et al DW Comparison - UNF", "Absolute", c("Reis Dry Weight", 0.72))

#     Filtered DRY WEIGHT Adj. R^2
UTILITY.heatmap_helper(
  "Reis et al DW Adj-R2--FIL", "image_output\\heatmaps\\REIS", DATA.REIS_models_and_abbrs.DRYWEIGHT, "Adj. R^2", DATA.filtered, "REIS et al Comparisons - DW - FIL.csv", "table_output\\REIS", "REIS et al DW Comparison - FIL", "Absolute", c("Reis Dry Weight", 0.72))

# Prior Work - LUCENA et. al. ---------------

# POWER Law Comparisons

cat("Generating LUCENA et. al. Comparison Heatmaps...\n")

DATA.LUCENA_models_and_abbrs.POWER <- read.csv("parameters\\heatmap_sources\\LUCENA_power_law.csv")
# Adj. R^2
UTILITY.heatmap_helper("Lucena et. al. POWER Adj-R2--UNF", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.POWER, "Adj. R^2", DATA.unfiltered,   "LUCENA et al Comparisons - POWER Adj R^2 - UNF.csv", "table_output\\LUCENA", "LUCENA et al POWER Adj R^2 Comparison - UNF", "Absolute")
UTILITY.heatmap_helper("Lucena et. al. POWER Adj-R2--FIL", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.POWER, "Adj. R^2", DATA.filtered,   "LUCENA et al Comparisons - POWER Adj R^2 - FIL.csv", "table_output\\LUCENA", "LUCENA et al POWER Adj R^2 Comparison - FIL", "Absolute")
# SBC
UTILITY.heatmap_helper("Lucena et. al. POWER SBC--UNF", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.POWER, "SBC", DATA.unfiltered, "LUCENA et al Comparisons - POWER SBC - UNF.csv", "table_output\\LUCENA", "LUCENA et al POWER SBC Comparison - UNF", "Relative to Model")
UTILITY.heatmap_helper("Lucena et. al. POWER SBC--FIL", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.POWER, "SBC", DATA.filtered,   "LUCENA et al Comparisons - POWER SBC - FIL.csv", "table_output\\LUCENA", "LUCENA et al POWER SBC Comparison - FIL", "Relative to Model")


# GAMMA Model Comparisons
DATA.LUCENA_models_and_abbrs.GAMMA <- read.csv("parameters\\heatmap_sources\\LUCENA_gamma_law.csv")
# Adj. R^2
UTILITY.heatmap_helper("Lucena et. al. GAMMA Adj-R2--UNF", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.GAMMA, "Adj. R^2", DATA.unfiltered,   "LUCENA et al Comparisons - GAMMA Adj R^2 - UNF.csv", "table_output\\LUCENA", "LUCENA et al GAMMA Adj R^2 Comparison - UNF", "Absolute")
UTILITY.heatmap_helper("Lucena et. al. GAMMA Adj-R2--FIL", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.GAMMA, "Adj. R^2", DATA.filtered,   "LUCENA et al Comparisons - GAMMA Adj R^2 - FIL.csv", "table_output\\LUCENA", "LUCENA et al GAMMA Adj R^2 Comparison - FIL", "Absolute")
# SBC
UTILITY.heatmap_helper("Lucena et. al. GAMMA SBC--UNF", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.GAMMA, "SBC", DATA.unfiltered, "LUCENA et al Comparisons - GAMMA SBC - UNF.csv", "table_output\\LUCENA", "LUCENA et al GAMMA SBC Comparison - UNF", "Relative to Model")
UTILITY.heatmap_helper("Lucena et. al. GAMMA SBC--FIL", "image_output\\heatmaps\\LUCENA", DATA.LUCENA_models_and_abbrs.GAMMA, "SBC", DATA.filtered,   "LUCENA et al Comparisons - GAMMA SBC - FIL.csv", "table_output\\LUCENA", "LUCENA et al GAMMA SBC Comparison - FIL", "Relative to Model")

# Prior Work - COLOGGERO & PARRERA --------------- TODO

cat("Generating Cologgero & Parrera Comparison Heatmaps...\n")

DATA.COLOGGERO_PARRERA_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\caloggero_parrera_drybyrect.csv")

# Adj. R^2
UTILITY.heatmap_helper("Cologgero & Parrera Adj-R2--UNF", "image_output\\heatmaps\\COLOGGERO_PARRERA", DATA.COLOGGERO_PARRERA_models_and_abbrs, "Adj. R^2", DATA.unfiltered,   "Cologgero & Parrera Comparisons - Adj R^2 - UNF.csv", "table_output\\COLOGGERO_PARRERA", "COLOGGERO & PARRERA Adj R^2 Comparison - UNF", "Absolute")
UTILITY.heatmap_helper("Cologgero & Parrera Adj-R2--FIL", "image_output\\heatmaps\\COLOGGERO_PARRERA", DATA.COLOGGERO_PARRERA_models_and_abbrs, "Adj. R^2", DATA.filtered,   "Cologgero & Parrera Comparisons - Adj R^2 - FIL.csv", "table_output\\COLOGGERO_PARRERA", "COLOGGERO & PARRERA Adj R^2 Comparison - FIL", "Absolute")
# SBC
UTILITY.heatmap_helper("Cologgero & Parrera SBC--UNF", "image_output\\heatmaps\\COLOGGERO_PARRERA", DATA.COLOGGERO_PARRERA_models_and_abbrs, "SBC", DATA.unfiltered,   "Cologgero & Parrera Comparisons - SBC - UNF.csv", "table_output\\COLOGGERO_PARRERA", "COLOGGERO & PARRERA SBC Comparison - UNF", "Relative to Model")
UTILITY.heatmap_helper("Cologgero & Parrera SBC--FIL", "image_output\\heatmaps\\COLOGGERO_PARRERA", DATA.COLOGGERO_PARRERA_models_and_abbrs, "SBC", DATA.filtered,   "Cologgero & Parrera Comparisons - SBC - FIL.csv", "table_output\\COLOGGERO_PARRERA", "COLOGGERO & PARRERA SBC Comparison - FIL", "Relative to Model")

# Prior Work - DE CORTAZAR & NOBEL --------------- TODO

cat("Generating de Cortazar & Nobel Comparison Heatmaps...\n")

DATA.DE_CORTAZAR_NOBEL_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\de_cortazar_nobel_ofarea.csv")

# Adj. R^2
UTILITY.heatmap_helper("de Cortazar & Nobel Adj-R2--UNF", "image_output\\heatmaps\\DE_CORTAZAR_NOBEL", DATA.DE_CORTAZAR_NOBEL_models_and_abbrs, "Adj. R^2", DATA.unfiltered,   "de Cortazar and Nobel Comparisons - Adj R^2 - UNF.csv", "table_output\\DE_CORTAZAR_NOBEL", "DE CORTAZAR & NOBEL Adj R^2 Comparison - UNF", "Absolute")
UTILITY.heatmap_helper("de Cortazar & Nobel Adj-R2--FIL", "image_output\\heatmaps\\DE_CORTAZAR_NOBEL", DATA.DE_CORTAZAR_NOBEL_models_and_abbrs, "Adj. R^2", DATA.filtered,   "de Cortazar and Nobel Comparisons - Adj R^2 - FIL.csv", "table_output\\DE_CORTAZAR_NOBEL", "DE CORTAZAR & NOBEL Adj R^2 Comparison - FIL", "Absolute")
# SBC
UTILITY.heatmap_helper("de Cortazar & Nobel SBC--UNF", "image_output\\heatmaps\\DE_CORTAZAR_NOBEL", DATA.DE_CORTAZAR_NOBEL_models_and_abbrs, "SBC", DATA.unfiltered,   "de Cortazar and Nobel Comparisons - SBC - UNF.csv", "table_output\\DE_CORTAZAR_NOBEL", "DE CORTAZAR & NOBEL SBC Comparison - UNF", "Relative to Model")
UTILITY.heatmap_helper("de Cortazar & Nobel SBC--FIL", "image_output\\heatmaps\\DE_CORTAZAR_NOBEL", DATA.DE_CORTAZAR_NOBEL_models_and_abbrs, "SBC", DATA.filtered,   "de Cortazar and Nobel Comparisons - SBC - FIL.csv", "table_output\\DE_CORTAZAR_NOBEL", "DE CORTAZAR & NOBEL SBC Comparison - FIL", "Relative to Model")


# Prior Work - LUO & NOBEL ---------------
cat("Generating Luo & Nobel Comparison Heatmaps...\n")

DATA.LUO_NOBEL_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\luo_nobel_drybyvol.csv")

# Adj. R^2
UTILITY.heatmap_helper("Luo & Nobel Adj-R2--UNF", "image_output\\heatmaps\\LUO_NOBEL", DATA.LUO_NOBEL_models_and_abbrs, "Adj. R^2", DATA.unfiltered,   "Luo and Nobel Comparisons - Adj R^2 - UNF.csv", "table_output\\LUO_NOBEL", "LUO & NOBEL Adj R^2 Comparison - UNF", "Absolute")
UTILITY.heatmap_helper("Luo & Nobel Adj-R2--FIL", "image_output\\heatmaps\\LUO_NOBEL", DATA.LUO_NOBEL_models_and_abbrs, "Adj. R^2", DATA.filtered,   "Luo and Nobel Comparisons - Adj R^2 - FIL.csv", "table_output\\LUO_NOBEL", "LUO & NOBEL Adj R^2 Comparison - FIL", "Absolute")
# SBC
UTILITY.heatmap_helper("Luo & Nobel SBC--UNF", "image_output\\heatmaps\\LUO_NOBEL", DATA.LUO_NOBEL_models_and_abbrs, "SBC", DATA.unfiltered,   "Luo and Nobel Comparisons - SBC - UNF.csv", "table_output\\LUO_NOBEL", "LUO & NOBEL SBC Comparison - UNF", "Relative to Model")
UTILITY.heatmap_helper("Luo & Nobel SBC--FIL", "image_output\\heatmaps\\LUO_NOBEL", DATA.LUO_NOBEL_models_and_abbrs, "SBC", DATA.filtered,   "Luo and Nobel Comparisons - SBC - FIL.csv", "table_output\\LUO_NOBEL", "LUO & NOBEL SBC Comparison - FIL", "Relative to Model")


# Prior Work - INGLESE et. al. ---------------
cat("Generating Inglese et. al. Comparison Heatmaps...\n")

DATA.INGLESE_models_and_abbrs <- read.csv("parameters\\heatmap_sources\\inglese_AbyD2L2.csv")

# Adj. R^2
UTILITY.heatmap_helper("Inglese et. al. Adj-R2--UNF", "image_output\\heatmaps\\INGLESE", DATA.INGLESE_models_and_abbrs, "Adj. R^2", DATA.unfiltered,   "Inglese et al Comparisons - Adj R^2 - UNF.csv", "table_output\\INGLESE", "INGLESE et al Adj R^2 Comparison - UNF", "Absolute")
UTILITY.heatmap_helper("Inglese et. al. Adj-R2--FIL", "image_output\\heatmaps\\INGLESE", DATA.INGLESE_models_and_abbrs, "Adj. R^2", DATA.filtered,   "Inglese et al Comparisons - Adj R^2 - FIL.csv", "table_output\\INGLESE", "INGLESE et al Adj R^2 Comparison - FIL", "Absolute")
# SBC
UTILITY.heatmap_helper("Inglese et. al. SBC--UNF", "image_output\\heatmaps\\INGLESE", DATA.INGLESE_models_and_abbrs, "SBC", DATA.unfiltered,   "Inglese et al Comparisons - SBC - UNF.csv", "table_output\\INGLESE", "INGLESE et al SBC Comparison - UNF", "Relative to Model")
UTILITY.heatmap_helper("Inglese et. al. SBC--FIL", "image_output\\heatmaps\\INGLESE", DATA.INGLESE_models_and_abbrs, "SBC", DATA.filtered,   "Inglese et al Comparisons - SBC - FIL.csv", "table_output\\INGLESE", "INGLESE et al SBC Comparison - FIL", "Relative to Model")

# Prior Work - Silva et. al. ---------------




UTILITY.nonlinear_model_SILVA_bolton <- function(source_data) {
  
  filtered_data=source_data%>%filter(Area != "NA")
  # print(nrow(filtered_data))
  
  SILVA_hw <- function(data) {
    return(nls(data=data, "Area~a*(1-exp(b*height*width))/-b", start=list(a=100,b=0.001)))
  }
  SILVA_peri <- function(data) {
    return(nls(data=data, "Area~a*(1-exp(b*Pade_Peri))/-b", start=list(a=100,b=0.001)))
  }
  
  SILVA_hw_model_ALL = SILVA_hw(filtered_data)
  SILVA_peri_model_ALL = SILVA_peri(filtered_data)
  
  approx_r2_nonlinear <- function(data, model) {
    pred_values = predict(model)
    depvar = as.character(formula(model))[2]
    actu_values = data[[depvar]]
    r2 = 1-(sum((actu_values - pred_values)^2)/sum((actu_values - mean(actu_values))^2))
    print(r2)
    return(r2)
  }
  
  SILVA_hw_df <- data.frame(
    accession = "All",
    R2 = approx_r2_nonlinear(filtered_data, SILVA_hw_model_ALL),
    SBC = BIC(SILVA_hw_model_ALL)
  )
  SILVA_peri_df <- data.frame(
    accession = "All",
    R2 = approx_r2_nonlinear(filtered_data, SILVA_peri_model_ALL),
    SBC = BIC(SILVA_peri_model_ALL)
  )
  
  for (acc in DATA.possible_accessions) {
    print(acc)
    #print(approx_r2_nonlinear(filtered_data, SILVA_hw_model))
    data_acc = filtered_data %>% filter(accession == acc)
    SILVA_hw_df   %<>% add_row(accession = as.character(acc), R2 = approx_r2_nonlinear(data_acc, SILVA_hw(data_acc)),   SBC = BIC(SILVA_hw(data_acc)))
    SILVA_peri_df %<>% add_row(accession = as.character(acc), R2 = approx_r2_nonlinear(data_acc, SILVA_peri(data_acc)), SBC = BIC(SILVA_peri(data_acc)))
  }
  print(SILVA_hw_df)
  print(SILVA_peri_df)
}

UTILITY.nonlinear_model_SILVA_bolton(DATA.unfiltered)

# Program end ------------------------

ENDTIME <- Sys.time()
TOTALTIME <- ENDTIME - STARTTIME

cat("Total Time Elapsed: ")
cat(TOTALTIME)
cat(" Minutes")