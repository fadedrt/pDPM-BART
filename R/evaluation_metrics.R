prediction <- function(x, y, test_x, npost, trees) {
  y_mean <- mean(y)
  y_sd   <- sd(y)
  
  prediction_train <- matrix(NA, ncol = npost, nrow = nrow(x))
  prediction_test  <- matrix(NA, ncol = npost, nrow = nrow(test_x))
  
  # Generate predictions using the tree ensembles
  for (i in 1:npost) {
    prediction_train[, i] <- get_predictions(trees[[i]], x, single_tree = FALSE)
    prediction_test[, i]  <- get_predictions(trees[[i]], test_x, single_tree = FALSE)
  }
  
  # Rescale predictions back to the original response scale
  f_hat_train <- prediction_train * y_sd + y_mean
  f_hat_test  <- prediction_test * y_sd + y_mean
  
  return(list(
    f_hat_train      = f_hat_train, 
    f_hat_test       = f_hat_test, 
    f_hat_train_mean = rowMeans(f_hat_train), 
    f_hat_test_mean  = rowMeans(f_hat_test)
  ))
}
