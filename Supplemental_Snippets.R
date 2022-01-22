# This is a collection of code snippets which can be used to derive other 
#   information in the paper after running the program once and creating a copy 
#   of the test dataset.

# WARNING: these snippets CANNOT be run as their own program, and will NOT 
#   return a useful result if executed. These snippets must be used in the
#   command line of your R environment.


# Creating a copy of the test dataset:

# Snippet used to obtain individual Adj. R^2 values for a given model.
#   Requires an ACCESSION value 
#     SWITCHBOARD.strALLACCESSIONS 
#     -or- 
#     Any 3-digit code corresponding to an accession
#   and a THRESHOLD value
#     Ranging from 100 (for everything), 95, and 90. Only these values will 
#     work.

#   Example for all accessions and all data:
summary(lm(Area ~ Theo_Area, data = ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA, SWITCHBOARD.strALLACCESSIONS, ACCESSORY.ErrorCleaning(db_shell, 100)))) 

#   Example for accession 242 on a 95% filter:

summary(lm(Area ~ Theo_Area, data = ACCESSORY.DataSubset(SWITCHBOARD.strALLDATA, 242, ACCESSORY.ErrorCleaning(db_shell, 95)))) 

