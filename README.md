# pDPM-BART

This repository contains code to implement the methodology described in the paper *"Bayesian Additive Regression Trees with Power Dirichlet Process Mixture Errors"*. 

The Bayesian tree structure and framework in this code are adapted from the primary functions of [CSP-BART](https://github.com/ebprado/CSP-BART/tree/main). 

Optimization Note: To enhance computational efficiency, we have specifically accelerated three primary bottleneck functions by re-implementing them with optimized C++ backends. The corresponding source code for these accelerated functions can be found in R/pdpmbart_backend.R.

These weights are **not model parameters** and are **not used during posterior sampling**. They are computed only for evaluating predictive performance and are stored in `Pi_store` for convenience.
## Example

```r
# Load required libraries
library(zeallot) # For multiple assignment (%<-%)
library(Rcpp)    # Seamless R and C++ integration for computational acceleration

# Function to simulate Friedman's dataset with Gaussian errors
sim_fried <- function(N, P, sigma) {
  X <- matrix(runif(N * P), nrow = N)  
  mu <- 10 * sin(pi * X[, 1] * X[, 2]) + 
        20 * (X[, 3] - 0.5)^2 + 
        10 * X[, 4] + 
        5 * X[, 5]
  
  Y <- mu + rnorm(N, 0, sigma)  
  return(list(X = X, Y = Y, mu = mu))
}

# -------------------------
# Simulate training and test datasets
# -------------------------
set.seed(1)
c(x, y, mu) %<-% sim_fried(250, 5, sigma = 1)

c(test_x, test_y, test_mu) %<-% sim_fried(100, 5, sigma = 1)

# -------------------------
# Fit pdpmBART model
# -------------------------
pdpmbart_fit <- bart_mixnormal_pdpmbart(
  y, x,
  ntrees = 50,
  nburn  = 1000,
  npost  = 2000
)

# -------------------------
# Generate predictions
# -------------------------
c(f_hat_train, f_hat_test, f_hat_train_mean, f_hat_test_mean) %<-% 
  prediction(x, y, test_x, npost = 2000, trees = pdpmbart_fit$trees)

# -------------------------
# Evaluation metrics
# -------------------------
cat("pdpmBART RMSE:", round(sqrt(mean((f_hat_test_mean - test_y)^2)), 2), "\n")
cat("pdpmBART MAE:", round(mean(abs(f_hat_test_mean - test_y)), 2), "\n")
cat("pdpmBART MAPE:", round(mean(abs(f_hat_test_mean - test_y) / abs(test_y)) * 100, 2), "%\n")

cat("pdpmBART LPML:",
    calculate_lpml_pdpmbart(
      y = y,
      f_hat_train = pdpmbart_fit$bart_store,
      Mu_store = pdpmbart_fit$Mu_store,
      Sigma2_store = pdpmbart_fit$Sigma2_store,
      Pi_store = pdpmbart_fit$Pi_store
    ), "\n")

dpm_waic_res <- calculate_waic_pdpmbart(pdpmbart_fit, y)
cat("pdpmBART WAIC:", round(dpm_waic_res$WAIC, 2), "\n")
