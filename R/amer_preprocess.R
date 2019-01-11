# amer_preprocess.R
# Preprocess Marks and Amers SDF files, chemical compound data
# Sept-Nov 2016

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

setwd("C:/R-files/QSAR")  # point to where file is located

sdf187 <- read.SDFset("QSAR_187.sdf")    # load it in a data structure
#sdf187 <- read.SDFset("QSAR_187_NEW.sdf")

#sdfTrain <- read.SDFset("Train150.sdf")  
#sdfTest  <- read.SDFset("Test37.sdf") 
# how many compounds do we have in this data structure?
#length(sdfTest)
#length(sdfTrain)

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


# ------------- RANDOM SPLIT OF DATA INTO TRAINING AND TEST SETS  --------------

set.seed(123)  # set seed so others can duplicate if necessary
n <- nrow(X)
index <- sample(1:n, size = round(0.80*n), replace=FALSE)
train <- X[index,]
test  <- X[-index,]

xtrain <- train[,2:ncol(X)]; # do not use 1st X value, (the Activity) in the training data
ytrain <- train[,1] # 1st column is the Activity value i.e. the training label

xtest <- test[,2:ncol(X)]; # do not use 1st column of X value, (the Activity) in the test data
ytest <- test[,1] # 1st column 

# ------------- FEATURE SELECTION USING CARET --------------
set.seed(12213)
xtrain <- scale(xtrain)
index <- createFolds(xtrain, k = 10, returnTrain = T)

ctrl <- rfeControl(functions = lmFuncs, method = "repeatedcv",repeats = 5,verbose = FALSE,index = index)
subsets <- c(1,2,4,5,10,15,20,25)
lmprofile <- rfe(xtrain,ytrain,sizes = subsets, rfeControl = ctrl)

for(i in 1:length(index)){
  crrltn = cor(xtrain[index[[i]],])     
  findCorrelation(crrltn, cutoff = .90, names = T, verbose = T)
}

# find highly correlated predictor variables and remove
x2 <- model.matrix(~., data = xtrain)[,-1]
nzv <- nearZeroVar(x2)
x3 <- x2[, -nzv]
corr_mat <- cor(x3)
too_high <- findCorrelation(corr_mat, cutoff = .9)
x4 <- x3[, -too_high]
c(ncol(x2), ncol(x3), ncol(x4))

crrltn = findCorrelation(x4, cutoff = .90)
if (length(crrltn) != 0)
  x4 <- x4[,-crrltn]

xtrain<- as.data.frame(xtrain)
xtest<- as.data.frame(xtest)
ytrain <- as.matrix(ytrain)
ytest <- as.matrix(ytest)
x4 <- as.data.frame(x4)
  
plot(t,y,type="l",col = 2)
lines(t[-c(1:2)],ps, col=3)
legend(1.5, 80, c("y", "pred"), cex=1.5, fill=2:3)



