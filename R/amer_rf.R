# amer_rf.R
# Can only work if the data used to train and test the NN is loaded
library(randomForest)
library(caret)

load("C:/common_laptop/R-files/QSAR/cv_rf_nn.RData")

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

rf_pred <- predict(rf_model, xtest)


# Now use 10-fold cross-validation
#cvdata <- rbind(trainingdata,testdata)
#cvtest <- rbind(ytrain,ytest)

fitcontrol <- trainControl(method="repeatedcv",
                           number=10, repeats = 3, 
                           verbose = TRUE, savePredictions = TRUE)

rf_fit <- train(x=cvdata[,1:15], y=cvdata[,16],
           data=cvdata,
           trControl=fitcontrol,
           method="rf",
           proximity = TRUE,
           importance=TRUE,
           ntree=5000)

print(rf_fit)

# test RF on 5000 compounds
rf_pred <- predict(rf_fit, xtest)


# train NN with cv
nncontrol <- trainControl(method="cv", number=10)
nn_fit <- train(cvdata[,1:15],cvdata[,16], method="neuralnet", 
               algorithm = "rprop+", 
               learningrate = 0.25,
               stepmax = 1e+06,
               act.fct = 'logistic',
               threshold = 0.01, 
               lifesign="full",
               trControl=nncontrol,
               tuneGrid = expand.grid(.layer1=c(10:30), .layer2=c(10:30),.layer3=c(0)))
               
#pred5k <- compute(nn_fit,xtest)  # only for neuralnet without Cv
#pred5k <- pred5k$net.result

pred5k <- predict(nn_fit, newdata = xtest)  # neuralnet trained using CV and Caret and full 187 compounds
# joinf rf pred with nn pred

joint <- cbind(names(rf_pred),rf_pred,pred5k)
colnames(joint)<-c("Compound","RF","NN")
write.table(joint,"C:/R-files/QSAR/old5K.csv",sep=",",row.names=FALSE)

save(rf_fit, nn_fit, cvdata, xtest, file = "cv_rf_nn.RData")




