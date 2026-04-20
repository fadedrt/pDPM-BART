# =================================================================
# Prediction Function for BART-based Models
# =================================================================
prediction <- function(x, y, test_x, npost, trees) {
  y_mean <- mean(y)
  y_sd   <- sd(y)
  
  prediction_train <- matrix(NA, ncol = npost, nrow = nrow(x))
  prediction_test  <- matrix(NA, ncol = npost, nrow = nrow(test_x))
  
  # Generate predictions from posterior tree ensembles
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


# =================================================================
# WAIC Computation for pDPM-BART Model (Numerically Stable Version)
# =================================================================
calculate_waic_pdpmbart <- function(fit, y) {
  # 1. Safely extract dimensions from posterior samples
  npost <- length(fit$Mu_store)
  n <- length(y)
  
  # Defensive check for valid posterior samples
  if (npost == 0 || is.null(npost)) {
    stop("Error: No posterior samples found in fit object. Please check fit$Mu_store.")
  }
  
  # 2. Initialize log-likelihood matrix (npost × n)
  log_lik_matrix <- matrix(NA, nrow = npost, ncol = n)
  
  cat("Computing WAIC (stable and vectorized version)...\n")
  
  # 3. Loop over posterior samples to compute pointwise log-likelihood
  for (s in 1:npost) {
    # Extract parameters from iteration s (already rescaled to original y space)
    mu_s    <- fit$Mu_store[[s]]
    sigma_s <- sqrt(fit$Sigma2_store[[s]])
    pi_s    <- fit$Pi_store[[s]]
    K_s     <- fit$K_store[s]
    f_s     <- fit$bart_store[s, ]
    
    # Residuals under current posterior draw
    resids <- y - f_s
    
    # Matrix of log mixture component likelihoods
    log_mix_mat <- matrix(-Inf, nrow = n, ncol = K_s)
    
    for (k in 1:K_s) {
      if (pi_s[k] > 0) {
        # Vectorized Gaussian log-density computation
        log_mix_mat[, k] <- log(pi_s[k]) + dnorm(resids, mean = mu_s[k], sd = sigma_s[k], log = TRUE)
      }
    }
    
    # Log-Sum-Exp trick for numerical stability
    max_vals <- apply(log_mix_mat, 1, max)
    log_lik_matrix[s, ] <- max_vals + log(rowSums(exp(log_mix_mat - max_vals)))
  }
  
  # 4. Compute WAIC components
  
  # (1) lppd: log pointwise predictive density (goodness-of-fit, higher is better)
  lppd_i <- apply(log_lik_matrix, 2, function(x) {
    max_x <- max(x)
    max_x + log(mean(exp(x - max_x)))
  })
  lppd <- sum(lppd_i)
  
  # (2) p_waic: effective number of parameters (complexity penalty, lower is better)
  p_waic_i <- apply(log_lik_matrix, 2, var)
  p_waic <- sum(p_waic_i)
  
  # (3) WAIC: final criterion (lower is better)
  waic <- -2 * (lppd - p_waic)
  
  cat("--------------------------\n")
  cat("lppd (fit, higher is better):   ", round(lppd, 3), "\n")
  cat("p_waic (penalty, lower is better):   ", round(p_waic, 3), "\n")
  cat("WAIC (overall, lower is better):     ", round(waic, 3), "\n")
  cat("--------------------------\n")
  
  return(list(
    WAIC   = waic, 
    lppd   = lppd, 
    p_waic = p_waic
  ))
}


# =================================================================
# LPML Computation for pDPM-BART Model
# Based on Conditional Predictive Ordinates (CPO)
# =================================================================
calculate_lpml_pdpmbart <- function(y, f_hat_train, Mu_store, Sigma2_store, Pi_store) {
  npost <- nrow(f_hat_train)
  n <- length(y)
  log_cpo <- numeric(n)
  
  for (i in 1:n) {
    log_pred <- numeric(npost)
    
    for (s in 1:npost) {
      # Residual under posterior draw s
      mu_resid <- y[i] - f_hat_train[s, i]
      
      mu_s     <- Mu_store[[s]]
      sigma_s  <- sqrt(Sigma2_store[[s]])
      pi_s     <- Pi_store[[s]]
      
      # Log-density under mixture components
      log_mix_comps <- log(pi_s) + dnorm(mu_resid, mean = mu_s, sd = sigma_s, log = TRUE)
      
      # Log-Sum-Exp aggregation for mixture density
      max_val <- max(log_mix_comps)
      log_pred[s] <- max_val + log(sum(exp(log_mix_comps - max_val)))
    }
    
    # Harmonic mean estimator for CPO (numerically stabilized)
    u <- -log_pred
    max_u <- max(u)
    log_cpo[i] <- -(max_u + log(mean(exp(u - max_u))))
  }
  
  # LPML is the sum of log CPO values
  return(sum(log_cpo))
}
