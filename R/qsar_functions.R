# qsar_functions.R
#------------ SOME FUNCTIONS FOR CALCULATING R2 ---------------------   

rmse <- function(error) {
  cat("\n RMSE", sqrt(mean(error^2))) # RMSE
  cat("\n MAE", mean(abs(error))) # MAE
  
}


R2 <- function(x, y) {
  #res <- caret::postResample(x, y)
  #rsq <- res[2]
  
  avr_y_actual <- mean(ytestB$PIC50)
  ss_total     <- sum((ytestB$PIC50 - avr_y_actual)^2)
  ss_regression <- sum((svm_pred - avr_y_actual)^2)
  ss_residuals <- sum((ytestB$PIC50 - svm_pred)^2)
  r2 <- 1 - ss_residuals / ss_total
  
  return(r2)
  #1 - sum((ytestB$PIC50-svm_pred)^2)/sum((ytestB$PIC50-mean(ytestB$PIC50))^2)
}

rto.estimates <- function(x, y) {
  b1 <- sum(x * y) / sum(x^2)
  ssr <- b1^2 * sum(x^2)
  sse <- sum(y^2) - ssr
  mse <- sse / (length(x) - 1)
  msr <- ssr / 1
  res.std.err <- sqrt(mse)
  f.stat <- msr / mse
  std.error <- sqrt(mse / sum(x^2))
  
  r2 <- ssr / (sse + ssr)
  
  adj.r2 <- r2 - (1 - r2) * (2 - 1) / (length(x) - 1)
  
  res <- data.frame(rbind(b1, res.std.err, f.stat, std.error, r2, adj.r2))
  rownames(res) <- c('b1', 'Residual Standard Error', 'F-Statistic', 'b1 Standard Error', 
                     'r-squared', 'Adjusted r-squared')
  colnames(res) <- 'Estimates'
  
  print(format(res, scientific = FALSE, digits = 3))
}


# As suggested by the reviewers some further measures used in QSAR maybe more appropriate
# to judge validity of the model as R2 and RSME do not always provide a good estimate.

# Roy measure: 
roy <- function(values){
  r2 <- boot::corr(d=xy) ^ 2
  cat("\nRoy ", r2) # Roy
  
}


# Tropsha measure:  # Tropsha-Golbraikh
tropsha <- function(values){
  Ytr <- mean(trainingdata$PIC50)  # mean value of the dependent variable using training data
  Yte <- mean(testdata$PIC50)
  tp <- 1-(sum((values$PIC50 - values$Pred)^2)/sum((values$PIC50 - Ytr)^2))
  cat("\nTropsha ", tp,"\n")
  return(tp)
}

q2 <- function(values){
  Ytr <- mean(trainingdata$PIC50)  # mean value of the dependent variable using training data
  Yte <- mean(testdata$PIC50)
  tp <- 1-(sum((values$PIC50 - values$Pred)^2)/sum((values$PIC50 - Ytr)^2))
  cat("\nQ2 ", tp,"\n")
  return(tp)
}



