#fit187 = glmnet(xtrain, ytrain, alpha = 1, nlambda = 150, standardize=TRUE,family="gaussian")
# fit # fit contains four bits of information:
# 1. num of lambdas[1,],2. num of variables used(Df),3. % of variation accounted for, 4. Lambda
#
#predicted <- predict(fit187,s=c(0.01155),newx=xtest, type = "response")
#cbind(predicted,ytest)


# identify the key independent variables using coef() getting all nonzero coeffs.
#setlambda <- 0.01155
#vals <- coef(fit187, s = setlambda)[which(coef(fit187, s = setlambda) != 0)]
#cvnames <- colnames(X)[which(coef(fit187, s = setlambda) != 0)]
#independents <- which(coef(fit187, s = setlambda) != 0)
#cbind(cvnames,vals)
