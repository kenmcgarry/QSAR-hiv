# qsar_5k.R
# Jan 2017, Feb 2019
# Get the 2 x 5,000 compound files into good shape for analysis

setwd("C:/common_laptop/R-files/QSAR")  # point to where file is located
#load("C:/common_laptop/R-files/QSAR/neuralnet187.RData")   # load up trained neural network
sdf5k1 <- read.SDFset("5000compounds.sdf")   # load in the old 5K compounds
sdf5k2 <- read.SDFset("2nd-5000compounds.sdf") # load in the new 5K compounds

valid <- validSDF(sdf5k1); 
sdf5k1 <- sdf5k1[valid] # remove invalid data if we have any.
cat("\n ",length(sdf5k1)," valid compounds in 1st batch")

valid <- validSDF(sdf5k2); 
sdf5k2 <- sdf5k2[valid] # remove invalid data if we have any.
cat("\n ",length(sdf5k2)," valid compounds in 2nd batch")


# get the data out of SDF into matrix form for datamining
blockmatrix5 <- datablock2ma(datablocklist=datablock(sdf5k1)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix5) # Splits to numeric and character matrix 
data5 <- numchar[1]
data5 <-as.data.frame(data5)
names(data5) <- gsub("numMA.", "", names(data5)) # gets rid of "numMA." which is prefixed to column names
dim(data5)
data5 <- na.omit(data5)  # get rid of NA if any
X <- as.matrix(data5)

# ------------- REMOVE UNINFORMATIVE VARIABLES i.e. ANY VARIABLE WITH 90% OR MORE ZEROS OR ONES --------
# ------------- if we dont it causes PCA to crash -------------
limit <- nrow(X)*.90
dim(X)
X <- X[, colSums(X != 1) > limit] 
X <- X[, colSums(X != 0) > limit] 
dim(X)

# -------------  PCA bit, NOTE: will generate a warning --------------------------
pcaSDF <- prcomp(X[,2:ncol(X)], cor = TRUE,scale=TRUE, center=TRUE)
X1 <- pcaSDF$x
n <- nrow(X1)

X1 <- X1[,1:numberPC]

########## try models on the compounds
nn_pred <- compute(nn_model,X1)
rf_pred <- predict(rf_model,X1)
svm_pred <-predict(svm_model,X1)
rbf_pred <- predict(rbf_model,X1)
pls_pred <- predict(pls_model,X1)

nn_temp <- data.frame(Compounds=rownames(nn_pred$net.result),PIC50=nn_pred$net.result)
nn_temp <- nn_temp %>% 
  dplyr::arrange(desc(PIC50))
head(nn_temp,10)

rf_temp <- data.frame(Compounds=names(rf_pred),PIC50=rf_pred)
rf_temp <- rf_temp %>% 
  dplyr::arrange(desc(PIC50))
head(rf_temp,10)

svm_temp <- data.frame(Compounds=names(svm_pred),PIC50=svm_pred)
svm_temp <- svm_temp %>% 
  dplyr::arrange(desc(PIC50))
head(svm_temp,10)

rbf_temp <- data.frame(Compounds=rownames(rbf_pred),PIC50=rbf_pred)
rbf_temp <- rbf_temp %>% 
  dplyr::arrange(desc(PIC50))
head(rbf_temp,10)

pls_temp <- data.frame(Compounds=rownames(pls_pred),PIC50=pls_pred)
rbf_temp <- pls_temp %>% 
  dplyr::arrange(desc(PIC50))
head(pls_temp,10)


