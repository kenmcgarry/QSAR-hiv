# hiv.R
# Analysis of Marks and Amers SDF files, chemical compound data
# Sept 2016

library(ChemmineR)
library(fmcsR)
library(glmnet)

setwd("C:/R-files/QSAR")

sdftest <- read.SDFset("qsar_test.sdf")
sdfset <- read.SDFset("QSAR39.sdf") 
sdf61 <- read.SDFset("QSAR_61.sdf") 

# how many molecules do we have in each file?
length(sdfset)
length(sdftest)
length(sdf61)


# get the TRAINING data into form for datamining
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix) # Splits to numeric and character matrix 
data1 <- numchar[1]
data1<-as.data.frame(data1)
dim(data1)

# get the TEST data into form for datamining
blockmatrix2 <- datablock2ma(datablocklist=datablock(sdftest)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix2) # Splits to numeric and character matrix 
data2 <- numchar[1]
data2<-as.data.frame(data2)
# names(data3) <- gsub("numMA.", "", names(data1)) # gets rid of "numMA." which is appended to column names
dim(data2)

data1[,1] <- log(data1[,1]) # get log of Activity (EC50) for TRAINING
data2[,1] <- log(data2[,1]) # get log of Activity (EC50) for TEST

# -log(data1[,1]) # what to do with initial negative values? These cause mathematical errors!

#----------------------------- glmnet stuff ----------------------------------------
# Y=response, X=predictors, both need to be in matrix format NOT a dataframe.
# compound 38 has NA values (missing data) causing crashes, next few lines excludes it.
# X1[,c('numMA.vsurf_S','numMA.vsurf_W8')]   # this line grabs matrix columns by name if I have to.

# TRAINING data
Y <- as.matrix(data1[c(1:37),1])
X <- as.matrix(data1[c(1:37),2:328])

# TEST DATUM
X1 <- as.matrix(data1[39,2:328]) #Only one test datum (CMP39)
Y1 <- as.matrix(data1[39,1])  # This is the correct answer

# TEST data
X2 <- as.matrix(data2[c(1:length(sdftest)),2:length(data2)]) # 22 test cases
Y2 <- as.matrix(data2[,1])  # These are the correct answers

grid<- 10^seq(10,-2,length=100)

fit = glmnet(X, Y, alpha = 1, nlambda = 150,standardize=TRUE,family="gaussian")

# FIT: the number of nonzero coefficients (Df), the percent (of null) deviance explained (%dev) and the value of Lambda.
fit

# Using a Lambda of 0.03825 (found in fit structure) which has 33 non-zero coeffs giving an explained variance of 0.99350
predict(fit,s=0.03825,newx=X1)

dim(coef(fit))

coef.exact = coef(fit, s = 0.5, exact = TRUE)
coef.apprx = coef(fit, s = 0.5, exact = FALSE)
cbind2(coef.exact, coef.apprx)
plot(fit, xvar = "dev", label = TRUE)
plot(fit, xvar = "lambda", label = TRUE)

# ------------ glmnet on QSAR_61.SDF ------------------------------------
setwd("C:/R-files/QSAR")

sdf61 <- read.SDFset("QSAR_61.sdf") 

# get the data into form for datamining
blockmatrix3 <- datablock2ma(datablocklist=datablock(sdf61)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix3) # Splits to numeric and character matrix 
data3 <- numchar[1]
data3<-as.data.frame(data3)
names(data3) <- gsub("numMA.", "", names(data3)) # gets rid of "numMA." which is prefixed to column names
dim(data3)

data3[,1] <- -log10(data3[,1]) # get -log of Activity (EC50) for TRAINING and TESTING
data3 <- na.omit(data3)  # get rid of NA if any
X <- as.matrix(data3[c(1:nrow(data3)),1:ncol(data3)])

Y<- X[,1] # 1st column is the Activity value i.e. the training label
X<- X[,2:ncol(X)] # do use 1st X value, (the Activity) in the training data

fit61 = glmnet(X[1:50,], Y[1:50], alpha = 1, nlambda = 150,standardize=TRUE,family="gaussian")
fit61

# identify the key independent variables 
vals <- coef(cv1, s = "lambda.min")[which(coef(cv1, s = "lambda.min") != 0)]
cvnames <- colnames(X)[which(coef(cv1, s = "lambda.min") != 0)]
independents <- which(coef(cv1, s = "lambda.min") != 0)
cbind(cvnames,vals)

# now do a linear regression model with variables identified by Lasso
X<-as.data.frame(X)
Y<-as.numeric(as.matrix(Y))
reg <- lm(Y[1:50] ~ X[1:50,1]+X[1:50,21]+X[1:50,33]+X[1:50,35]+X[1:50,42]+
            X[1:50,57]+X[1:50,87]+X[1:50,136]+X[1:50,199]+X[1:50,241]+X[1:50,288]+X[1:50,296])
summary(reg)
lmpred <- predict.lm(reg,new=X[51:61,independents],interval="confidence")
cbind(lmpred,Y[51:61])

# ----------------------------------------------------------------------------

# test the model using compounds 40-61
predicted <- predict(fit61,s="lambda.min",newx=X[51:61,])
cbind2(predicted, Y[51:61])

# use cross-validation this time.
cvfit = cv.glmnet(X[1:39,], Y[1:39], type.measure = "mse", nfolds = 10)
cvpred <- predict(cvfit,type = "response",s = "lambda.min",newx=X[40:61,])
cbind2(cvpred, Y[40:61])

#----------- glmnet with cross-fold validation -----------
foldid=sample(1:10,size=length(Y[1:60]),replace=TRUE)
cv1=cv.glmnet(X[1:60,],Y[1:60],foldid=foldid,alpha=1,type.measure = "mse")
cv.5=cv.glmnet(X[1:60,],Y[1:60],foldid=foldid,alpha=.5,type.measure = "mse")
cv0=cv.glmnet(X[1:60,],Y[1:60],foldid=foldid,alpha=0,type.measure = "mse")

par(mfrow=c(2,2))
plot(cv1,xlab="log(Lambda) for Alpha=1");plot(cv.5,xlab="log(Lambda) for Alpha=0.5");plot(cv0,xlab="log(Lambda) for Alpha=0")
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name,xlim=c(-4,7),ylim=c(0,7))
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
legend("topright",legend=c("alpha= 1","alpha= 0.5","alpha=0"),pch=19,col=c("red","grey","blue"),cex = 0.75,pt.cex = .8,y.intersp=0.3)

coef(cv0, s = "lambda.min")
coef(cv.5, s = "lambda.min")
coef(cv1, s = "lambda.min")
#------------------------------------------------------------------------------------
data2 <- data1[,which(!apply(data1==0,2,all))]  # remove columns with all zeros
data2 <- data2[,which(!apply(data2==1,2,all))]  # remove columns with all ones

results1 <- prcomp(data2[2:37,],scale=TRUE) # data are scaled with mean=0 and variance=1
summary(results1)
names(results1)
  
pcaCharts(results1)
biplot(results1,scale=0, cex=.7)
          
          
# Eigenvalues
eig <- (results1$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.vals <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
head(eig.vals)

# now do regression with PC's
y<- data2[1:36,1]
reg <- lm(y~results1$x[,1]+results1$x[,2]+results1$x[,3]+results1$x[,4]+results1$x[,5]+results1$x[,6])
summary(reg)




#plot(sdfset[1:4], print=FALSE) # Plots structures to R graphics device 
sdf.visualize(sdfset[1:4])

sdfid(sdfset[1:10])
cid(sdfset[1:10])

propma <- atomcountMA(sdfset, addH=FALSE) 
boxplot(propma, col="lightgreen", main="Atom Frequency") 

MW(sdfset[1:10], addH=FALSE)
MF(sdfset[1:10], addH=FALSE) 
groups(sdfset[1:10], groups="fctgroup", type="countMA") 


# EC50 is dependent variable i.e. activity
# get negative log of activity

# I got this function to plot 2 x 2 charts from....
# http://rstudio-pubs-static.s3.amazonaws.com/27823_dbc155ba66444eae9eb0a6bacb36824f.html
pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  screeplot(x)
  screeplot(x,type="l")
  par(mfrow=c(1,1))
}




# divide up training and test data randomly
#set.seed(1)
#indexes = sample(1:nrow(X), size=0.7*nrow(X))
# Split data
#test = X[indexes,]
#dim(test)  
#train = X[-indexes,]
#dim(train) 
#dt = sort(sample(nrow(X), nrow(X)*.7))
#train<-X[dt,]
#test<-X[-dt,]

