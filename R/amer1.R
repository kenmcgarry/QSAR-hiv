# hiv.R
# Analysis of Marks and Amers SDF files, chemical compound data
# Sept 2016

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
sdfTrain <- read.SDFset("Train150.sdf")  
sdfTest  <- read.SDFset("Test37.sdf") 
# how many compounds do we have in this data structure?
length(sdfTest)
length(sdfTrain)


# get the data out of SDF into form for datamining
blockmatrix1 <- datablock2ma(datablocklist=datablock(sdf187)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix1) # Splits to numeric and character matrix 
data1 <- numchar[1]
data1 <-as.data.frame(data1)
names(data1) <- gsub("numMA.", "", names(data1)) # gets rid of "numMA." which is prefixed to column names
dim(data1)
data1 <- na.omit(data1)  # get rid of NA if any
X <- as.matrix(data1)

# Random split of data into train and test sets
set.seed(123)
#data(X)
n <- nrow(X)
index <- sample(1:n, size = round(0.80*n), replace=FALSE)
train <- X[index,]
test  <- X[-index,]

xtrain <- train[,2:ncol(X)]; # do not use 1st X value, (the Activity) in the training data
ytrain <- train[,1] # 1st column is the Activity value i.e. the training label

xtest <- test[,2:ncol(X)]; # do not use 1st column of X value, (the Activity) in the test data
ytest <- test[,1] # 1st column 

# ------------ Process QSAR.SDF ------------------------------------

# get the TRAINING data into form for datamining
blockmatrix3 <- datablock2ma(datablocklist=datablock(sdfTrain)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix3) # Splits to numeric and character matrix 
data3 <- numchar[1]
data3<-as.data.frame(data3)
names(data3) <- gsub("numMA.", "", names(data3)) # gets rid of "numMA." which is prefixed to column names
dim(data3)
data3 <- na.omit(data3)  # get rid of NA if any
X <- as.matrix(data3)

# get the TEST data into form for datamining
blockmatrix4 <- datablock2ma(datablocklist=datablock(sdfTest)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix4) # Splits to numeric and character matrix 
data4 <- numchar[1]
data4<-as.data.frame(data4)
names(data4) <- gsub("numMA.", "", names(data4)) # gets rid of "numMA." which is prefixed to column names
dim(data4)
data4 <- na.omit(data4)  # get rid of NA if any
Y <- as.matrix(data4)

xtrain <- X[,2:ncol(X)]; # do not use 1st X value, (the Activity) in the training data
ytrain <- X[,1] # 1st column is the Activity value i.e. the training label

xtest <- Y[,2:ncol(Y)]; # do not use 1st column of Y value, (the Activity) in the test data
ytest <- Y[,1] # 1st column is the Activity value i.e. the test label

# AMER'S VARIABLES DECIDED BY LITERATURE
xtrain<- as.data.frame(xtrain)
xtest<- as.data.frame(xtest)
ytrain <- as.matrix(ytrain)
ytest <- as.matrix(ytest)
colnames(ytrain) <- c("PIC50")

#reg <- lm(ytrain~xtrain$logP.o.w.) # carefully note "~.,"
#summary(reg)


lmpred <- predict.lm(reg,new=xtest) #pass test data thru and see what it predicts Activity
cbind(lmpred,ytest) # compare predictions with true Activity values
plot(cbind(ytest,lmpred))
pd<-cbind(ytest,lmpred)
plotres(reg)

# Partial Least Squares model fitting
ytrain <- as.data.frame(ytrain)
xtrain <- as.data.frame(xtrain)
ytest <- as.data.frame(ytest)
xtest <- as.data.frame(xtest)
newdata<-cbind(ytrain,xtrain)
newtest <-cbind(ytest,xtest)
plsfit <- plsr(PIC50~ALogP+Molecular_Weight+Num_AromaticRings+
                     Num_RotatableBonds+Num_H_Donors+Num_H_Acceptors+
                     Molecular_FractionalPolarSurfaceArea,data=newdata, vaidation="CV")

plsfit <-  plsr(ytrain~.,data=newdata,validation="CV")
  
summary(plsfit)

plspredict <- predict(plsfit, ncomp = 50, newdata =xtest )
plseval <- data.frame(obs=ytest,pred=plspredict,row.names = NULL)
colnames(plseval)<-c("obs","pred")
defaultSummary(plseval)

cbind(plspredict,ytest,ytest-plspredict)

plot(RMSEP(plsfit), legendpos = "topright")
plot(plsfit, ncomp = 4, asp = 1, line = TRUE)
validationplot(plsfit,val.type="MSEP")

# save data in text file format
mydata<-cbind(ytrain,xtrain)
write.table(mydata, "amerdata.tsv", sep="\t",row.names=FALSE) 



# ===============   IMPROVED MODEL BY EXCLUDING COMPOUNDS =========================
# how to setup X and Y to exclude certain compounds
Y1<- (as.data.frame(Y))
X1<- (as.data.frame(X))

excluded <- c(6, 30, 31, 36, 57, 58, 59, 75, 157) # compounds to be excluded
Y1 <- Y1[-excluded, ]
X1 <- X[ -excluded, ]
X1 <- as.data.frame(X1)
theEND <- length(Y1)
theSTART <-  nrow(X1) - 15  # use last 15 compounds are the test set

# use only independent variables with significant p-vals (predictors)
data=X1[1:theSTART-1,sigind]
sigind <- c(12,23,24,25,29,31,38,48,240,266,313)
rownames(X1) <- NULL
reg <- lm(Y1[1:theSTART-1]~.,data=X1[1:theSTART-1,sigind]) # carefully note "~.," means use all variables
summary(reg)

lmpred <- predict.lm(reg,new=X1[theSTART:theEND,sigind]) #pass test data thru and see how it predicts the Activity
cbind(lmpred,Y1[theSTART:theEND]) # compare predictions with true Activity values
residual <- (lmpred - Y1[theSTART:theEND])
cbind(lmpred,Y1[theSTART:theEND],residual)
plotres(reg)

robust <- rlm(as.matrix(ytrain)~.,data=xtrain) # do robust regression
summary(robust)
rsquare(robust)

# polynomial regression
Y1<-as.numeric(as.matrix(Y1))
fitpoly <- lm(Y1~poly(X1$vsurf_CW6,4)+poly(X1$GCUT_SLOGP_3,4)+poly(X1$GCUT_PEOE_3,4),data=X1)

fitpoly <- lm(Y1~I(X1$vsurf_CW6^2)+
                I(X1$vsurf_CW6^4)+
                I(X1$vsurf_CW6)+
                I(X1$GCUT_SLOGP_3)+
                #I(X1$GCUT_SLOGP_3^4)+
                I(X1$GCUT_PEOE_3^4),data=X1)

summary(fitpoly)

# save data in text file format
colnames(ytrain) <- c("Activity")
mydata<-cbind(ytrain,xtrain)
write.table(mydata, "amerdata.tsv", sep="\t",row.names=FALSE) 

# does the data work in reverse?
newdata <- read.table("amerdata.tsv",sep="\t",header=TRUE)
newdata[1:162,2:12]<- (newdata[1:162,2:12])
newreg <- lm(newdata[1:150,1]~.,data=newdata[1:150,2:9]) # carefully note "~.," means use all variables
newpred <- predict.lm(newreg,new=newdata[140:150,2:12]) 
residual <- (newpred - newdata[140:150,1])
cbind(newpred,newdata[140:150,1],residual)
summary(newreg)




# Now try clustering compounds to determine similarity
apset <- sdf2ap(sdfTrain)
dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="distmat.rda", quiet=TRUE) 
load("distmat.rda") 
hc <- hclust(as.dist(distmat), method="single") 
hc[["labels"]] <- cid(apset) # Assign correct item labels 
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=T) 

heatmap.2(1-distmat, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), 
          col=colorpanel(40, "darkblue", "green", "white"),density.info="none", trace="none") 


# Will PCA reveal anything about the compounds?
pcaSDF <- princomp(xtrain, cor = TRUE)
#pcaSDF <- prcomp(X1, cor = TRUE)
plot(pcaSDF, type = "l")
biplot(pcaSDF)
summary(pcaSDF)
sdfHC <- hclust(dist(pcaSDF$scores), method = "ward.D2")
sdfClusters <- cutree(sdfHC, k = 8)  # enforce number of clusters

compDf <- data.frame(pcaSDF$scores, "cluster" = factor(sdfClusters)) # add cluster to data frame of scores
compDf <- transform(compDf, cluster_name = paste("Cluster",sdfClusters))

p1 <- ggplot(compDf,aes(x=Comp.1, y=Comp.2)) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(aes(color = cluster), alpha = 0.55, size = 5) +
  xlab("PC1") +
  ylab("PC2") + 
  xlim(-5, 6) + 
  ggtitle("PCA clusters from hierarchical clustering of SDF compound data") 
p1 + geom_text(aes(y = Comp.2 + 0.25, label = rownames(compDf)))


# Try different clustering models

# Determine number of clusters
wss <- (nrow(xtrain)-1)*sum(apply(xtrain,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(xtrain,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares") 

# Model based approaches assume a variety of data models and apply maximum likelihood estimation and Bayes 
# criteria to identify the most likely model and number of clusters. Specifically, the Mclust( ) function 
# in the mclust package selects the optimal model according to BIC for EM initialized by hierarchical 
# clustering for parameterized Gaussian mixture models. 
library(mclust)
fit <- Mclust(X1[1:theSTART-1,sigind])
plot(fit) # plot results
summary(fit) # display the best model

# The pvclust( ) function in the pvclust package provides p-values for hierarchical clustering based 
# on multiscale bootstrap resampling. Clusters that are highly supported by the data will have large p values
library(pvclust)
fit <- pvclust(t(X1[1:theSTART-1,sigind]), method.hclust="ward",method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95) 

# Ward Hierarchical Clustering
d <- dist(X1[1:theSTART-1,sigind], method = "euclidean") # distance matrix
fitward <- hclust(d, method="ward")
plot(fitward) # display dendogram
groups <- cutree(fitward, k=7) # cut tree into 7 clusters
# draw dendogram with red borders around the 7 clusters
rect.hclust(fitward, k=7, border="red") 

# K-Means Clustering with 7 clusters
fitkm <- kmeans(X1[1:theSTART-1,sigind], 7)
# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster)
clusplot(X1[1:theSTART-1,sigind], fitkm$cluster, color=TRUE, shade=TRUE,labels=2, lines=0)

# Centroid Plot against 1st two discriminant functions
library(fpc)
plotcluster(X1[1:theSTART-1,sigind], fitkm$cluster) 

# comparing 2 cluster solutions
library(fpc)
cluster.stats(d, fitkm$cluster, fitward$cluster) 



# ----------------------------------------------------------------------------
xtrain<-as.matrix(xtrain)
ytrain<-(as.matrix(ytrain))
xtest<-as.matrix(xtest)
ytest<-(as.matrix(ytest))

#-------- glmnet with k-fold cross-validation for three models (Lasso, Elastic and Ridge) -------
foldid=sample(1:20,size=length(ytrain),replace=TRUE)
cv1=cv.glmnet(xtrain,ytrain,nfold=10,alpha=1,type.measure = "mae")
cv.5=cv.glmnet(xtrain,ytrain,foldid=foldid,alpha=.5,type.measure = "mse")
cv0=cv.glmnet(xtrain,ytrain,foldid=foldid,alpha=0.0,type.measure = "mse")

par(mfrow=c(2,2))
plot(cv1,xlab="log(Lambda) for Alpha=1");plot(cv.5,xlab="log(Lambda) for Alpha=0.5");plot(cv0,xlab="log(Lambda) for Alpha=0")
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name,xlim=c(-6,4),ylim=c(2,4))
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
legend("topright",legend=c("alpha= 1","alpha= 0.5","alpha=0"),pch=19,col=c("red","grey","blue"),cex = 0.75,pt.cex = .8,y.intersp=0.3)

plotres(cv1) # plot residuals
plotres(cv0) # plot residuals
plotres(cv.5) # plot residuals


# the residual version as per Teixeria et al paper 2009
cvpredicted<- predict(cv1,newx=xtest,type="response",s=cv1$lambda.min)
residual <- (cvpredicted - ytest)
results <- data.frame(predicted=cvpredicted,actual=ytest,residual=residual)
names(results) <- c("predicted", "actual", "residual")
results

rsquare(NULL,ytest,cvpredicted)
get_r2_ss(ytest,cvpredicted,w = rep(1, length(ytest)))



# ------------- FEATURE SELECTION USING CARET --------------
ctrl <- rfeControl(functions = lmFuncs, method = "repeatedcv",repeats = 5,verbose = FALSE)
subsets <- c(1,2,4,5,10,15,20,25)
lmProfile <- rfe(xtrain, ytrain,sizes = subsets, rfeControl = ctrl)

#------------ SOME FUNCTIONS FOR CALCULATING R2 ---------------------   
rsquare = function(object=NULL, y, fitted.y){
  # get y and fitted.y from a lm or glm object
  if (! is.null(object)) {
    fitted.y = fitted(object)
    if (class(fitted.y) == "numeric")
      y = object$model[[1]]
    else
      y = object$model[,1]
  }
  
  # compute coefficient of determination
  if (class(fitted.y) == "numeric") {
    return(cor(y, fitted.y)^2)
  } else {
    R2 = double(ncol(y))
    for (ic in 1:ncol(y))
      R2[ic] = cor(y[,ic], fitted.y[,ic])^2
    return(R2)  }
}

get_r2_cor <- function(y, y_pred, w) {
  # Calculate R2 using the correlation coefficient method
  xy = cbind(y, y_pred)
  return(boot::corr(d=xy, w=w) ^ 2)
}

get_r2_ss <- function(y, y_pred, w) {
  # Calculate R2 using the Sum of Squares method
  # https://en.wikipedia.org/wiki/Coefficient_of_determination#Definitions
  ss_residual = sum(w * (y - y_pred) ^ 2)
  ss_total = sum(w * (y - weighted.mean(y, w)) ^ 2)
  return(1 - ss_residual / ss_total)
}

get_r2_likehood <- function(model, model_intercept, n) {
  # Calculate R2 using the generalized (likelihood) method
  # https://en.wikipedia.org/wiki/Coefficient_of_determination#Generalized_R2
  L_0 = exp(as.numeric(logLik(model_intercept)))
  L_null = exp(as.numeric(logLik(model)))
  return(1 - (L_0 / L_null) ^ (2 / n))
}

simulate <- function(weighted = T, n = 50, seed = 0) {
  # Randomly generate data, perform regression, and return the r-squared values
  # produced by various methods.
  
  # Simulate x (the predictor), y (the outcome), and w (the observation weights)
  set.seed(seed)
  x = runif(n)
  y = runif(n)
  if (weighted) {w = runif(n)} else {w = rep(1, n)}
  
  # Fit linear regression models and compute predictions
  model_intercept = lm(y ~ 1, weight=w)
  model = lm(y ~ x, weight=w)
  y_pred = predict(model)
  
  # Calculate and return the four R2 measures
  return(c(
    r2_cor = get_r2_cor(y, y_pred, w),
    r2_lm = summary(model)$r.squared,
    r2_likelihood = get_r2_likehood(model, model_intercept, n),
    r2_ss = get_r2_ss(y, y_pred, w)
  ))
}

# remove uninformative variables i.e. all or majority are zeros
clean_zeros <- function(thedata) {
  
  limit <- nrow(thedata)*.90
  crappy <- xtrain[, colSums(xtrain != 0) > limit] 
  
  #thedata <- thedata[,which(!apply(thedata==0,2,all))]
  
  #thedata <- thedata[, which(apply(thedata, 2, function(col) !any(table(col) > limit)))] 
  
  return(thedata)
}

