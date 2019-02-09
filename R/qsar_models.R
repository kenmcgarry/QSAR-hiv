# qsar_models.R
# Modified 10/1/19 for work comparing RF, MLP, SVM and Linear Regression models

# -------------- NEURAL NETWORK MODEL ---------------------
library("neuralnet")

n <- colnames(trainingdata)
f <- as.formula(paste("PIC50 ~", paste(n[!n %in% "PIC50"], collapse = " + ")))
net <- neuralnet(f,data=trainingdata, hidden=c(30,30), linear.output=T,threshold=0.001,
                 lifesign="minimal",algorithm="rprop+")

prednn <- compute(net,xtest)
results<-cbind(prednn$net.result,ytestB,ytestB-prednn$net.result)
colnames(results)<-c("neuralnet","PIC50","Residual")
results

# calculate R2, rmse and mae
xy<-cbind(ytestB,prednn$net.result)
boot::corr(d=xy ^ 2)
rmse(as.vector(unlist(ytestB))-as.vector(unlist(prednn$net.result)))

plot(xy,col='red',main='Actual vs predicted NN, (PCA=15 comps) R2=.84,RMSE=1.0',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

# -------------- SUPPORT VECTOR MACHINE (SVM) MODEL ---------------------



# -------------- RANDOM FOREST (RF) MODEL ---------------------

rf_model=randomForest(PIC50 ~., data=trainingdata,
                      importance=TRUE,
                      proximity = TRUE,
                      keep.forest=TRUE,
                      mtry=15,
                      ntree=2500)

print(rf_model)
round(importance(rf_model), 2)
varImpPlot(rf_model,main="",type=2,color="black",pch=16) 
rf_model
plot(rf_model)

# use random forest on old 5K and new 5K compound files
# rf_pred <- predict(rf_model, xtest)










