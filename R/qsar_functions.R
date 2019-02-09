# qsar_functions.R
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

simulate1 <- function(weighted = T, n = 50, seed = 0) {
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


rmse <- function(error) {
  cat("\n RMSE", sqrt(mean(error^2))) # RMSE
  cat("\n MAE", mean(abs(error))) # MAE
  
}



