# qsar_models_crossvalidated.R
# https://rstudio-pubs-static.s3.amazonaws.com/214482_7b51dfba415d4b6a8bb304083ce64de1.html
# https://github.com/julianhatwell/Statistical_Learning_Basics

# In response to reviewers,  modified 14/5/19 for cross-validation of RF, MLP, SVM by 
# using the CARET tuning mechanism "train() function" for optimizing model parameters
# 
# -------------- NEURAL NETWORK (MLP) MODEL ---------------------
library(caret)
library(neuralnet)
library(e1071)
library(randomForest)
library(plotmo)
library(xlsx)
library(dplyr)
library(e1071)
library(RSNNS)
library(kerndwd)
library(kernlab)



n <- colnames(trainingdata)
form <- as.formula(paste("PIC50 ~", paste(n[!n %in% "PIC50"], collapse = " + ")))

fitControl <- trainControl( # k-fold CV.
  method = "cv",
  verboseIter = TRUE,
  number = 5)

tunGrid <- expand.grid(data.frame(layer1 = 20, # optimum parameter range
                                  layer2 = 22,layer3=0)) 

tunGrid <- expand.grid(data.frame(layer1 = 18:25, layer2 = 18:25, layer3 = 0)) # takes 2 hours

# CV trained neural network
nn_cv <- caret::train(form, data = trainingdata,
                 method = 'neuralnet', 
                 trControl = fitControl,
                 linear.output=TRUE,
                 act.fct = "logistic",
                 threshold=0.001,
                 lifesign="full",
                 preProcess = c("center", "scale"),
                 stepmax = 1e+06,
                 algorithm="rprop+",
                 tuneGrid  = tunGrid)
nn_cv

nn_pred <- neuralnet::compute(nn_cv$finalModel,testdata[,1:numberPC])
results <- cbind(nn_pred$net.result,ytestB,ytestB-nn_pred$net.result)
colnames(results)<-c("neuralnet","PIC50","Residual")
head(results,20)

# calculate R2, rmse and mae
rm(xy)
xy<-cbind(nn_pred$net.result,ytestB,nn_pred$net.result-ytestB)
cat("\n RSQUARE",boot::corr(d=xy) ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(nn_pred$net.result)))

# Roy and Trophsa measures
colnames(xy) <- c("PIC50","Pred","Residual")
tropsha(xy)
roy(xy)

plot(xy,col='red',main='Actual vs predicted NN',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

# do the residual plot using plotres()
plotres(nn_cv,which=1:4)

write.xlsx(results, "mydata.xlsx") 

#########################################################################################
# -------------- RANDOM FOREST (RF) MODEL ---------------------
cat("\nTraining Random Forest model")

tunegrid <- expand.grid(mtry=c(5:15)) # only parameter that can be modified using train()

rf_cv <- caret::train(form, data=trainingdata, 
                   method="rf", tuneGrid=tunegrid,ntree=500,
                   trControl=fitControl)
print(rf_cv)
plot(rf_cv)

rf_pred <- predict(rf_cv, testdata[,1:numberPC])
results <- cbind(rf_pred,ytestB,ytestB-rf_pred)
colnames(results)<-c("RandomForest","PIC50","Residual")
head(results,20)
rm(xy)

xy<-cbind(ytestB,rf_pred)
cat("\n RSQUARE",boot::corr(d=xy) ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(rf_pred)))

# Roy and Trophsa measures
colnames(xy) <- c("PIC50","Pred")
tropsha(xy)
roy(xy)

plot(xy,col='red',main='Actual vs predicted Random Forest, (PCA=15 comps)',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='RF',pch=18,col='red', bty='n')
plotres(rf_cv)

write.xlsx(head(results,20), "c:mydata.xlsx") 

#########################################################################################
# -------------- SUPPORT VECTOR MACHINE (SVM) MODEL ---------------------
# https://rpubs.com/ezgi/classification

n <- colnames(trainingdata)
formzero <- as.formula(paste("PIC50 ~", paste(n[!n %in% "PIC50"], collapse = " + ")))
formzero <- paste(formzero, " +0")

fitControl2 <- trainControl( # k-fold CV.
  method = "repeatedcv",
  verboseIter = TRUE,
  repeats=10,
  number = 10)

fitControl3 <- trainControl(method = "LOOCV",verboseIter = TRUE)

grid_radial <- expand.grid(sigma = c(0.60,0.85,0.1,0.15,0.2,0.3),
                               C = c(1,2,3,4,5,6))

svm_cv <- caret::train(form, data = trainingdata, method = "svmRadial",#"svmLinear","svmRadial"
                    trControl=fitControl3,
                    tuneGrid=grid_radial,
                    metric="Rsquared",
                    tuneLength = 10,
                    savePred=TRUE,
                    preProcess = c("center", "scale"))#,
                    #tuneLength = 10)
svm_cv
#Predict using SVM regression
svm_pred <- predict(svm_cv, testdata[,1:numberPC])

rm(xy)
results<-cbind(ytestB,svm_pred)
head(results,20)
boot::corr(d=results)^2
rmse(as.vector(unlist(ytestB))-as.vector(unlist(svm_pred)))

# Roy and Trophsa measures
xy<-cbind(ytestB,svm_pred)
colnames(xy) <- c("PIC50","Pred")
tropsha(xy)

plot(xy,col='red',main='Actual vs predicted SVM, (PCA=15 comps)',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='SVM',pch=18,col='red', bty='n')
plotres(svm_cv)

results <- cbind(svm_pred,ytestB,ytestB-svm_pred)
colnames(results)<-c("SVM","PIC50","Residual")
head(results)

write.xlsx(head(results,20), "mydata.xlsx") 

#########################################################################################
# -------------- RADIAL BASIS FUNCTION (RBF) MODEL ---------------------
cat("\nTraining RBF model")
fitControlRBF <- trainControl( # k-fold CV.
  method = "cv",
  verboseIter = TRUE,
  number = 5)

rbf_cv <- caret::train(form, data = trainingdata, 
      method = 'krlsRadial', 
      trControl = fitControlRBF,
      tuneGrid = data.frame(.lambda = NA, 
                            .sigma  = c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50,100,500)),
      preProc = c("center","scale"))


rbf_cv <- caret::train(form, data = trainingdata, 
                       method = 'gaussprPoly',
                       trControl = fitControl,
                       #tuneGrid = data.frame(.sigma  = c(0.0001, 0.001, 0.01, 0.1, 1, 10 ,20)),
                       preProc = c("center", "scale"))

summary(rbf_cv)
rbf_pred <- predict(rbf_cv,testdata[,1:numberPC])
rm(xy)
xy <- cbind(ytestB,rbf_pred)
boot::corr(d=xy) ^ 2
rmse(as.vector(unlist(ytestB))-as.vector(unlist(rbf_pred)))
# Roy and Trophsa measures
colnames(xy) <- c("PIC50","Pred")
tropsha(xy)

plot(xy,col='red',main='Actual vs predicted RBF',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='RBF',pch=18,col='red', bty='n')
plotres(rbf_cv)

results <- cbind(rbf_pred,ytestB,ytestB-rbf_pred)
colnames(results)<-c("RBF","PIC50","Residual")
head(results)

write.xlsx(xy, "mydata.xlsx")

#########################################################################################
# -------------- PARTIAL LEAST SQUARES (PLS) MODEL ---------------------

pls_cv <- caret::train(form,data = trainingdata,
                          method = "pls",
                          preProc = c("center", "scale"),
                          tuneLength = 25,  
                          trControl = fitControl)

summary(pls_cv)
pls_pred <-  predict(pls_cv, testdata[,1:numberPC])
rm(xy)
xy <- cbind(ytestB,pls_pred)#[,1,numberPC])
boot::corr(d=xy) ^ 2
rmse(as.vector(unlist(ytestB))-as.vector(unlist(pls_pred)))

# Roy and Trophsa measures
colnames(xy) <- c("PIC50","Pred")
tropsha(xy)

plot(xy,col='red',main='Actual vs predicted PLS ',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='PLS',pch=18,col='red', bty='n')
plotres(pls_cv)

results <- cbind(pls_pred,ytestB,ytestB-pls_pred)
colnames(results)<-c("PLS","PIC50","Residual")
head(results)
write.xlsx(head(xy,20), "mydata.xlsx") 

#########################################################################################
# -------------- LINEAR REGRESSION (LM) MODEL ---------------------

lm_cv <- caret::train(form, data = trainingdata, 
               method = "lm",
               trControl = fitControl, 
               metric="Rsquared")
summary(lm_cv)
lm_pred <- predict(lm_cv, testdata[,1:numberPC])

rm(xy)
xy<-cbind(ytestB,lm_pred)
boot::corr(d=xy) ^ 2
rmse(as.vector(unlist(ytestB))-as.vector(unlist(lm_pred)))

# Roy and Trophsa measures
colnames(xy) <- c("PIC50","Pred")
tropsha(xy)

plot(xy,col='red',main='Actual vs predicted Linear Regression, (PCA=15 comps) ',pch=18,cex=0.7)
legend('bottomright',legend='Linear',pch=18,col='red', bty='n')
plotres(lm_cv$finalModel)


results <- cbind(lm_pred,ytestB,ytestB-lm_pred)
colnames(results)<-c("LinearRegress","PIC50","Residual")
head(results,20)


write.xlsx(head(xy,20), "mydata.xlsx") 




