# pDPM-BART

This repository contains code to implement the methodology described in the paper *"Bayesian Additive Regression Trees with Power Dirichlet Process Mixture Errors"*. 

The Bayesian tree structure and framework in this code are adapted from the primary functions of [CSP-BART](https://github.com/ebprado/CSP-BART/tree/main). 

Optimization Note: To enhance computational efficiency, we have specifically accelerated three primary bottleneck functions by re-implementing them with optimized C++ backends. The corresponding source code for these accelerated functions can be found in R/pdpmbart_backend.R.

## Example

```r
# Load required libraries
library(zeallot) # For multiple assignment (%<-%)
library(Rcpp)    # Seamless R and C++ integration for computational acceleration

# Function to simulate Friedman's dataset with mixture-of-normals errors
sim_fried <- function(N, P) {
  X <- matrix(runif(N * P), nrow = N)  # Random uniform predictors
  mu <- 10 * sin(pi * X[, 1] * X[, 2]) + 
        20 * (X[, 3] - 0.5)^2 + 
        10 * X[, 4] + 
        5 * X[, 5]  # True mean function
  
  # Generate mixture-of-normals residuals
  z <- sample(1:3, size = N, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  mu_k <- c(0, 2, -2)
  sigma_k <- c(1, 0.5, 1.5)
  
  eps <- rnorm(N, mean = mu_k[z], sd = sigma_k[z])
  
  Y <- mu + eps
  return(list(X = X, Y = Y, mu = mu))
}

# -------------------------
# Simulate training and test datasets
# -------------------------
set.seed(1)
c(x, y, mu) %<-% sim_fried(250, 5)

c(test_x, test_y, test_mu) %<-% sim_fried(100, 5)

# Fit pDPM-BART model
pdpmbart_fit <- bart_mixnormal_infinite_power_rr(y, x, ntrees = 50, nburn = 2500, npost = 2500)

# Generate predictions
c(f_hat_train, f_hat_test, f_hat_train_mean, f_hat_test_mean) %<-% 
  prediction(x, y, test_x, npost = 2500, trees = pdpmbart_fit$trees)

# Generate residual shift based on posterior mixture components
residual_shift_samples <- sapply(1:length(pdpmbart_fit$Mu_store), function(s) {
  mu_s <- pdpmbart_fit$Mu_store[[s]]
  pi_s <- pdpmbart_fit$Pi_store[[s]]
  sum(pi_s * mu_s)
})

r_shift <- mean(residual_shift_samples)

# Adjust predictions
final_predictions <- f_hat_test_mean + r_shift

# -------------------------
# Evaluation metrics
# -------------------------
cat("pDPM-BART RMSE:", round(sqrt(mean((final_predictions - test_y)^2)), 2), "\n")
cat("pDPM-BART MAE:", round(mean(abs(final_predictions - test_y)), 2), "\n")
cat("pDPM-BART MAPE:", round(mean(abs(final_predictions - test_y) / abs(test_y)) * 100, 2), "%\n")
cat("pDPM-BART LPML:", calculate_lpml_pdpmbart(
  y = y, 
  f_hat_train = pdpmbart_fit$bart_store,
  Mu_store = pdpmbart_fit$Mu_store, 
  Sigma2_store = pdpmbart_fit$Sigma2_store, 
  Pi_store = pdpmbart_fit$Pi_store
), "\n")
