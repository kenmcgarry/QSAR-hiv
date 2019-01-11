# amer_pca.R
library(ChemmineR)
library(dplyr)
library(neuralnet)
library(data.table)
library(qcc)
library("factoextra")
library("FactoMineR")
library(pcaMethods)
library(xtable)

setwd("C:/common_laptop/R-files/QSAR")  # point to where file is located
# load("C:/common_laptop/R-files/QSAR/amer_pca_data.RData")

sdf187 <- read.SDFset("QSAR_187.sdf")    # load in SDF 187 selected compounds
#sdf5k <- read.SDFset("5000compounds.sdf")   # load in the old 5K compounds
#sdf5k <- read.SDFset("2nd-5000compounds.sdf") # load in the new 5K compounds

valid <- validSDF(sdf187); sdf187 <- sdf187[valid] # remove invalid data if we have any.
length(sdf187)

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
n <- nrow(X1)

summary(pcaSDF)
screeplot(pcaSDF, type="lines",col=3)
#The loadings for the principal components are stored in:
pcaSDF$rotation # with princomp(): pca$loadings
plot(pcaSDF)

#biplot of first two principal components
biplot(pcaSDF,cex=0.8)
abline(h = 0, v = 0, lty = 2, col = 8)

pcavar <- pcaSDF$sdev^2  # variances
pareto.chart(pcavar, ylab="Variances")  # plot pareto chart
#library(qcc)
#library("factoextra")
#library("FactoMineR")
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

res.pca <- PCA(X[,2:ncol(X)], graph = FALSE)
print (res.pca) 
eigenvalues <- res.pca$eig # get eigenvalues
head(eigenvalues[, 1:2])
fviz_screeplot (res.pca, ncp=10) # % variance explained by each PC
individual_contributions <- fviz_pca_var(res.pca) # Variables factor map
# Contribution of variables on PC1 and PC2
individual_contributions <- fviz_contrib(res.pca, choice = "var", axes = 1:2)  # set axes=1 to consider only PC1


plot(head(individual_contributions$data[, 2], 10), type="b", pch=20, xaxt="n", ylab="Percentage % Contribution", xlab="Observation Row Number", main="Top 10 - Contribution of individuals on PC1 & PC2")
axis(1, at=1:10, labels=head(individual_contributions$data[, 1], 10))


# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 20)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 20)
#The total contribution to PC1 and PC2
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 20)
# The most important (or, contributing) variables can be highlighted 
# on the correlation plot 
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Change the transparency by contrib values
fviz_pca_var(res.pca, alpha.var = "contrib")

# graph of compounds
#fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Avoid text overlapping (slow if many points))

# biplot, a biplot can be interpreted as follow:
# an individual that is on the same side of a given variable has a high value 
# for this variable;
# an individual that is on the opposite side of a given variable has a low value 
# for this variable.

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")  # Individuals color

# not too sure the usefulness of this command?
res.pca <- PCA(X[,2:ncol(X)], ind.sup = 24:27, 
               quanti.sup = 11:12, quali.sup = 13, graph=FALSE)
fviz_pca_var(res.pca)


## Visualize variable with cos2 >= 0.6
fviz_pca_var(res.pca, select.var = list(cos2 = 0.6))
# Top 5 active variables with the highest cos2
fviz_pca_var(res.pca, select.var= list(cos2 = 5))

# Select by names
name <- list(name = c("SlogP", "logP.o.w.", "lip_acc","logS"))
fviz_pca_var(res.pca, select.var = name)
# top x contributing individuals and variables
fviz_pca_biplot(res.pca, select.ind = list(contrib = 20), 
                select.var = list(contrib = 20),
                ggtheme = theme_minimal())

# SCREE PLOT
# http://www.sthda.com/english/wiki/print.php?id=202
fviz_screeplot(res.pca, ncp=15,main="")
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2],15)
sum(eigenvalues[1:15,2])

####
#source("https://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")






