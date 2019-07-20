# amer_preprocess.R
# Preprocess Marks and Amers SDF files, chemical compound data
# Sept-Nov 2016
# Feb 2019 - updated.
# June 2019 - bug fixed

library(ChemmineR)
library(fmcsR)
library(glmnet)
library(plotmo)
library(scatterplot3d) 
library(gplots) 
library(ggplot2)
library(pls)
library(MASS)
library(caret)
library(randomForest)
library(dplyr)
library(neuralnet)
library(data.table)
library(qcc)
library(factoextra)
library(FactoMineR)
library(pcaMethods)
library(xtable)
library(e1071)
library(RSNNS)
library(xlsx)

sdf187 <- read.SDFset("QSAR_187.sdf")    # This is the 187 compounds chosen from literature

# get the data out of SDF into form for datamining
blockmatrix1 <- datablock2ma(datablocklist=datablock(sdf187)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix1) # Splits to numeric and character matrix 
data1 <- numchar[1]
data1 <-as.data.frame(data1)
names(data1) <- gsub("numMA.", "", names(data1)) # gets rid of "numMA." which is prefixed to column names
dim(data1)
data1 <- na.omit(data1)  # get rid of NA if any
X <- as.matrix(data1)

# ------------- REMOVE UNINFORMATIVE VARIABLES i.e. ANY VARIABLE WITH 90% OR MORE ZEROS OR ONES --------
limit <- nrow(X)*.90
dim(X)
X <- X[, colSums(X != 1) > limit] 
X <- X[, colSums(X != 0) > limit] 
dim(X)

# -------------  PCA bit, NOTE: will generate a warning --------------------------
pcaSDF <- prcomp(X[,2:ncol(X)], cor = TRUE,scale=TRUE, center=TRUE)
X1 <- pcaSDF$x
summary(pcaSDF)
n <- nrow(X1)

# ------------- RANDOM SPLIT OF DATA INTO TRAINING AND TEST SETS  --------------
index <- sample(1:n, size = round(0.80*n), replace=FALSE)
train <- X1[index,]
test  <- X1[-index,]
xtest <- test[,1:ncol(X1)]; # do not use 1st column of X value, (the Activity) in the test data
xtrain <- train[,1:ncol(X1)]; # do not use 1st X value, (the Activity) in the training data

trainB <- X1[index,]
testB  <- X1[-index,]
#ytrainB <- trainB[,1] # 1st column is the Activity value (PIC50) i.e. the training label
#ytestB <- testB[,1] # 1st column 
ytrainB <- X[index,1] # 1st column is the Activity value (PIC50)i.e. the training label
ytestB <-  X[-index,1] # 1st column 

xtrain<- as.data.frame((xtrain[,1:numberPC]))
xtest<- as.data.frame((xtest[,1:numberPC]))

ytrainB <- as.data.frame(ytrainB)
ytestB <- as.data.frame(ytestB)
colnames(ytrainB)<- "PIC50"
colnames(ytestB) <- "PIC50"

trainingdata <- cbind(xtrain,ytrainB)
testdata <- cbind(xtest,ytestB)

