#!/usr/bin/env Rscript

# These scripts were written to automate on a large scale the statistical analyses
# done in R, in order to save time - so that minor changes to analytical processes
# may be quickly reflected across the entire analysis. This allows hours of work
# to be repeated, with or without changes, in less than two minutes while analyzing 
# thousands of datapoints. 
#
# This also permits the data analysis process to be inspected by reviewers at each 
# step of the analysis process.
#
# The program is interspersed with diagnostic output messages ("cat(...)") that 
# output to the console and indicate the progression of the program.
# 
# Contact gangres@nevada.unr.edu for questions regarding the code's function and 
# behavior.


#This script runs four custom function libraries, a file-loading script, and a file- and figure-generating in order.
#UTIL, GRAPHER, COMPILE, and CALC must be run first, in no particular order.
#Then, "File Loading and Prep" and "File Generation" must be run, in that specific order.


# |} Establish working directory with target files.
# Workstation Setup -------------------------------------------------------

STARTTIME <- Sys.time()

cat("Setting Working Directory...\n")

#REMEMBER TO OMIT THIS PRIOR TO PROGRAM SUBMISSION
## Argument from the USER
### https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
### don't use hard coding. If you want to do hard coding, then use github location
setwd("C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\Script-v2")

cat("Locating Dataset File...\n")

## I think you might seperate the location of input. such as data/ or othernames
DATA_FILE <- 'PARL0_06202022.csv'

cat("--------------------------------\n")
cat("--------------------------------\n")
cat("LOADING UTILITY FUNCTIONS")
cat("\n")
## I think you might seperate the location of input such as script/ or bin/ 
source("UTIL-Functions.R")
cat("\n")

cat("--------------------------------\n")
cat("--------------------------------\n")
cat("LOADING CALCULATOR FUNCTIONS\n")
cat("\n")
## I think you might seperate the location of input such as script/ or bin/ 
source("CALC-Functions.R")
cat("\n")

cat("--------------------------------\n")
cat("--------------------------------\n")
cat("LOADING GRAPHER FUNCTIONS\n")
cat("\n")
## I think you might seperate the location of input such as script/ or bin/ 
source("GRAPHER-Functions.R")
cat("\n")

cat("--------------------------------\n")
cat("--------------------------------\n")
cat("LOADING COMPILATION FUNTIONS\n")
cat("\n")
## I think you might seperate the location of input such as script/ or bin/ 
source("COMPILE-Functions.R")
cat("\n")

cat("--------------------------------\n")
cat("--------------------------------\n")
cat("CREATING FILES\n")
cat("\n")
## I think you might seperate the location of input such as script/ or bin/ 
source("File Loading and Prep.R")
cat("\n")

cat("--------------------------------\n")
cat("--------------------------------\n")
cat("CREATING OUTPUT GRAPHS\n")
cat("\n")
source("File Generation.R")
cat("\n")

ENDTIME <- Sys.time()
TOTALTIME <- ENDTIME - STARTTIME

cat("Total Time Elapsed: ")
cat(TOTALTIME)
cat(" Minutes")
