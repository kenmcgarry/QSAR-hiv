# qsar_models.R
# Modified 9/2/19 for work comparing RF, MLP, SVM and Linear Regression models

# -------------- NEURAL NETWORK MODEL ---------------------
cat("\nTraining NN model")
n <- colnames(trainingdata)
form <- as.formula(paste("PIC50 ~", paste(n[!n %in% "PIC50"], collapse = " + ")))
nn_model <- neuralnet(form,data=trainingdata, 
                      hidden=c(25,25), 
                      linear.output=TRUE,
                      threshold=0.001,
                      stepmax = 1e+06,
                     lifesign="full",algorithm="rprop+")

nn_pred <- compute(nn_model,testdata[,1:numberPC])
results<-cbind(nn_pred$net.result,ytestB,ytestB-nn_pred$net.result)
colnames(results)<-c("neuralnet","PIC50","Residual")
#results

# calculate R2, rmse and mae
xy<-cbind(ytestB,nn_pred$net.result)
boot::corr(d=xy ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(nn_pred$net.result)))

plot(xy,col='red',main='Actual vs predicted NN, (PCA=15 comps) R2=.84,RMSE=1.0',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

# -------------- SUPPORT VECTOR MACHINE (SVM) MODEL ---------------------
cat("\nTraining SVM model")
svm_model <- svm(form,data=trainingdata,type = "eps-regression", scale = TRUE, center=TRUE)

#Predict using SVM regression
svm_pred <- predict(svm_model, testdata[,1:numberPC])
rm(xy)
xy<-cbind(ytestB,svm_pred)
boot::corr(d=xy ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(svm_pred)))
plot(xy,col='red',main='Actual vs predicted SVM, (PCA=15 comps) R2=,RMSE=3.3',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='SVM',pch=18,col='red', bty='n')

# -------------- RANDOM FOREST (RF) MODEL ---------------------
cat("\nTraining Random Forest model")
rf_model=randomForest(form,data=trainingdata,
                      importance=TRUE,
                      proximity = TRUE,
                      #keep.forest=TRUE,scale=TRUE,center=TRUE,
                      mtry=10,
                      ntree=1500)

rf_pred <- predict(rf_model, testdata[,1:numberPC])
rm(xy)
xy<-cbind(ytestB,rf_pred)
boot::corr(d=xy ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(rf_pred)))
plot(xy,col='red',main='Actual vs predicted Random Forest, (PCA=15 comps) R2=,RMSE=2.5',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='RF',pch=18,col='red', bty='n')


# -------------- RADIAL BASIS FUNCTION (RBF) MODEL ---------------------
cat("\nTraining RBF model")
rbf_model <- rbf(trainingdata[,1:numberPC], trainingdata[,numberPC+1], size=80, maxit=10000, 
                 initFuncParams=c(0, 1, 0, 0.01, 0.01), 
                 learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)
rbf_pred <- predict(rbf_model,testdata[,1:numberPC])
rm(xy)
xy <- cbind(ytestB,rbf_pred)
boot::corr(d=xy ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(rbf_pred)))
plot(xy,col='red',main='Actual vs predicted RBF, (PCA=15 comps) R2=,RMSE=0.01',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='RBF',pch=18,col='red', bty='n')


# -------------- PARTIAL LEAST SQUARES (PLS) MODEL ---------------------
cat("\nTraining PLS model")
pls_model <- plsr(form, data=trainingdata, scale=TRUE, center=TRUE)
pls_pred <-  predict(pls_model, testdata[,1:numberPC])
rm(xy)
xy <- cbind(ytestB,pls_pred[,1,numberPC])
boot::corr(d=xy ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(pls_pred)))
plot(xy,col='red',main='Actual vs predicted PLS, (PCA=15 comps) R2=,RMSE=0.01',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='PLS',pch=18,col='red', bty='n')


# -------------- LINEAR REGRESSION (LM) MODEL ---------------------
cat("\nTraining Linear Regression model")
lm_model <- lm(form, data=trainingdata,scale=TRUE,center=TRUE)
lm_pred <- predict(lm_model, testdata[,1:numberPC])
rm(xy)
xy<-cbind(ytestB,lm_pred)
boot::corr(d=xy ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(lm_pred)))
plot(xy,col='red',main='Actual vs predicted Linear Regression, (PCA=15 comps) R2=,RMSE=0.01',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='Linear',pch=18,col='red', bty='n')



