# qsar_5k.R
# Preprocess Marks and Amers SDF files, chemical compound data
# Jan 2017
# Get the 5,000 compound file into good shape for analysis

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

# ------------- USE THE TRAINED NEURALNET -----------------------
xtest  <- as.data.frame((X1[,1:15]))  # we use "only" the first 15 PC's
pred5k <- compute(net,xtest)
pred5k <- pred5k$net.result

write.table(pred5k,"C:/common_laptop/R-files/QSAR/new5K.csv",sep=",",row.names=TRUE)


# read in amers results from his model PIC50
amerpred <- read.SDFset("C:/R-files/QSAR/5000predicted_by_PLS.sdf")

blockmatrix6 <- datablock2ma(datablocklist=datablock(amerpred)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrix6) # Splits to numeric and character matrix 
data6 <- numchar[1]
data6 <-as.data.frame(data6)
names(data6) <- gsub("numMA.", "", names(data6)) # gets rid of "numMA." which is prefixed to column names
dim(data6)
data6 <- na.omit(data6)  # get rid of NA if any
AX <- as.matrix(data6)


# -------------- COMPARE BETWEEN NEURAL NET AND PLS MODELS ------------
#compare_results <-cbind(names(AX),pred5k,AX)  # does not work anymore 1/12/18

compare_results <-cbind(pred5k,AX)  # show NN and Amers side by side

colnames(compare_results)<-c("neuralnet","AmerModel")
diff_models <- as.numeric(compare_results[,2]) - as.numeric(compare_results[,1])
diff_models <-abs(diff_models)
amerdata<- as.numeric(compare_results[,2])
nndata <- as.numeric(compare_results[,1])
compounds <- rownames(AX)
alldata <- cbind(compounds,diff_models,amerdata,pred5k)
alldata <- data.frame(compounds,diff_models,amerdata,nndata) # REMEMBER THIS FORMAT!!!
rownames(alldata)<-NULL
write.table(alldata,"C:/R-files/QSAR/compare5K-old.csv",sep=",",row.names=FALSE)



#bestresults <- arrange(alldata,desc(nndata)) # If they only want one model then sort by maximum values of IC50,
                                                # in this example its sorted by the neuralnet

# use the filter function in dplyr to extract compounds with low model differences 
# (between neural and amers) # and high IC50 values - its a tradeoff. Write the 
# list to file or an error message if we cant find any with # the settings in filter() function.
bestresults <- filter(alldata, diff_models < .1, pred5k > 9)
if(nrow(bestresults) > 0){
  write.table(bestresults,"C:/R-files/QSAR/bestresults.csv",sep=",",row.names=FALSE)}else{
   write.table("Sorry I couldn't find any compounds with these settings for diff_models and IC50","C:/R-files/QSAR/bestresults.csv",sep=",",row.names=FALSE)}


