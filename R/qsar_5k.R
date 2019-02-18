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

# REMOVE superfluous data
rm(data1,data5,blockmatrix1,blockmatrix5,sdf187,sdf5k1,sdf5k2,X)


########## try models on the compounds
nn_pred <- compute(nn_model,X1)
rf_pred <- predict(rf_model,X1)
svm_pred <-predict(svm_model,X1)
rbf_pred <- predict(rbf_model,X1)
pls_pred <- predict(pls_model,(X1[,1:numberPC]))
lm_pred <- predict(lm_model, as.data.frame(X1[,1:numberPC]))

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
pls_temp <- pls_temp %>% 
  dplyr::arrange(desc(PIC50))
head(pls_temp,10)

lm_temp <- data.frame(Compounds=names(lm_pred),PIC50=lm_pred)
lm_temp <- lm_temp %>% 
  dplyr::arrange(desc(PIC50))
head(lm_temp,10)

## COMPOUNDS in common between models
# VENN: ensure models output in common format to allow comparisions with Venn diagram however
# a large amount of processing is required for the mlp, since neuralnet package is very different.
library(VennDiagram)
library(gplots)
topcomp <- 10

candidates_rf <-  head(rf_temp$Compounds,topcomp)
candidates_svm <- head(svm_temp$Compounds,topcomp)
candidates_rbf <- head(rbf_temp$Compounds,topcomp)
candidates_nn <- head(nn_temp$Compounds,topcomp)

# make a nice well structured data frame, for some reason svm and rbf only produce 1264 candidates from 1280 test data  
comp_model <- data.frame(rf =as.vector(candidates_rf),
                         svm=as.vector(candidates_svm), 
                         rbf=as.vector(candidates_rbf),
                         mlp=as.vector(candidates_nn))

comp_model[] <- lapply(comp_model, as.character)

commoncomps <- Reduce(intersect, list(candidates_nn,
                                         candidates_rf,
                                         candidates_rbf,
                                         candidates_svm))  # compounds identified by all four classifiers


######################################################################
plot.new()
venn.plot <- venn.diagram(list(comp_model$rf,comp_model$svm,comp_model$rbf,comp_model$mlp), 
                          NULL, 
                          fill=c("red", "blue","green","yellow"), 
                          alpha=c(0.5,0.5,0.5,0.5), 
                          cex = 2, 
                          cat.fontface=2, 
                          margins =c(10,10),
                          cat.cex=2,
                          category.names=c("RF", "SVM","RBF","MLP"))
grid.draw(venn.plot)



