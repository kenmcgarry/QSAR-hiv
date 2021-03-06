# qsar_main.R
# building MLP regression models using quantitative structure-activity relationships (QSAR) data
# Feb 2019 - updated.

setwd("C:/common_laptop/R-files/QSAR")  # point to where Rcode and data are located
numberPC <- 15
cat("\nnumber of PCs set to ",numberPC)

# 1. load in the data (187 compounds for train/test), performs PCA, divides 
#    them into train/test sets
source("qsar_data_preprocess.R")  

# 2. load in functions for RMSE etc calcs
source("qsar_functions.R")

# 3. train NN, SVM, RF and RBF on 187 compound data
source("qsar_models.R")
#source("qsar_models_crossvalidated.R")

# 4. load in subsets of ZINC database (2 x 5,000 coumpounds) to search for novel compounds
#   use the best model.
source("qsar_5k.R")


