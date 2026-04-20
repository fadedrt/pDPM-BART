#' MCMC Updates for pDPM-BART Residual Model Parameters
#' 
#' This script implements Metropolis-Hastings and Gibbs updates for the 
#' powered Dirichlet process mixture (pDPM) residual model used in BART.
#' It includes pseudo-likelihood evaluation under pCRP, parameter updates 
#' for alpha and r, and conjugate updates for Gaussian mixture components.
#' 
#' @param z Cluster assignments
#' @param r Power parameter in pCRP
#' @param alpha Concentration parameter
#' @param x Residual observations
#' @param Mu Cluster means
#' @param Sigma2 Cluster variances
#' @param mu_0 Prior mean
#' @param nu Prior precision scaling
#' @param a Shape parameter for inverse-gamma prior
#' @param b Rate parameter for inverse-gamma prior
#' @param K Number of clusters
#' 
#' @export


# =================================================================
# Core Engine: Leave-one-out pseudo-likelihood under pCRP
# Ensures consistency with cluster update mechanism in C++ layer
# =================================================================
log_pcrp_pseudo_lik <- function(z, r, alpha) {
  n <- length(z)
  log_prob <- 0
  
  # Initialize cluster counts
  n_k <- rep(0, max(z))
  tbl <- table(z)
  n_k[as.integer(names(tbl))] <- as.integer(tbl)
  
  # Leave-one-out evaluation for each observation
  for (i in 1:n) {
    cluster_i <- z[i]
    n_k[cluster_i] <- n_k[cluster_i] - 1 # remove current observation
    
    sum_n_r <- sum(n_k[n_k > 0]^r)
    denominator <- alpha + sum_n_r
    
    if (n_k[cluster_i] == 0) {
      numerator <- alpha          
    } else {
      numerator <- n_k[cluster_i]^r 
    }
    
    log_prob <- log_prob + log(numerator) - log(denominator)
    n_k[cluster_i] <- n_k[cluster_i] + 1 # restore
  }
  return(log_prob)
}

# =================================================================
# Metropolis-Hastings update for alpha (using pseudo-likelihood)
# =================================================================
sample_alpha_mh <- function(alpha_current, z, r, a_alpha, b_alpha, mh_proposal_sd = 0.5) {
  log_alpha_current <- log(alpha_current)
  log_alpha_prop <- rnorm(1, mean = log_alpha_current, sd = mh_proposal_sd)
  alpha_prop <- exp(log_alpha_prop)
  
  log_prior_current <- dgamma(alpha_current, shape = a_alpha, rate = b_alpha, log = TRUE)
  log_prior_prop    <- dgamma(alpha_prop,    shape = a_alpha, rate = b_alpha, log = TRUE)
  
  # Evaluate pseudo-likelihood on current partition
  log_lik_current <- log_pcrp_pseudo_lik(z, r, alpha_current)
  log_lik_prop    <- log_pcrp_pseudo_lik(z, r, alpha_prop)
  
  log_jacobian <- log_alpha_prop - log_alpha_current
  log_acc <- (log_lik_prop + log_prior_prop + log_jacobian) - (log_lik_current + log_prior_current)
  
  if (!is.na(log_acc) && log(runif(1)) < log_acc) return(alpha_prop)
  return(alpha_current)
}

# =================================================================
# Metropolis-Hastings update for power parameter r (shifted half-normal prior)
# =================================================================
sample_r_mh <- function(r_current, z, alpha, r_prior_sigma, mh_proposal_sd) {
  
  # Symmetric random walk proposal
  r_prop <- rnorm(1, mean = r_current, sd = mh_proposal_sd)
  
  # Enforce constraint: r > 1
  if (r_prop <= 1) return(r_current)
  
  # Prior: r - 1 ~ Normal(0, sigma)
  log_prior_current <- dnorm(r_current - 1, mean = 0, sd = r_prior_sigma, log = TRUE)
  log_prior_prop    <- dnorm(r_prop - 1,    mean = 0, sd = r_prior_sigma, log = TRUE)
  
  # Pseudo-likelihood evaluation
  log_lik_current <- log_pcrp_pseudo_lik(z, r_current, alpha)
  log_lik_prop    <- log_pcrp_pseudo_lik(z, r_prop, alpha)
  
  # MH acceptance ratio
  log_acc <- (log_lik_prop + log_prior_prop) - (log_lik_current + log_prior_current)
  
  if (!is.na(log_acc) && log(runif(1)) < log_acc) return(r_prop)
  
  return(r_current)
}



# =================================================================
# Update Gaussian mixture parameters (Mu and Sigma2)
# using Normal-Inverse-Gamma conjugacy
# =================================================================
update_sigmamu <- function(x, z, Mu, Sigma2, mu_0, nu, a, b, K) {
  Mu_sample <- numeric(K)
  Sigma2_sample <- numeric(K)
  
  # Ensure prior parameters are valid
  a <- max(a, 1e-6)  
  b <- max(b, 1e-6)  
  nu <- max(nu, 1e-6) 
  
  for (k in 1:K) {
    idx_k <- which(z == k)
    N_k <- length(idx_k)  
    
    # Case 1: Empty cluster (draw from prior)
    if (N_k == 0) {
      Sigma2_sample[k] <- 1 / rgamma(1, shape = a, rate = b)
      Sigma2_sample[k] <- pmin(pmax(Sigma2_sample[k], 1e-8), 1e8)
      Mu_sample[k] <- rnorm(1, mean = mu_0, sd = sqrt(Sigma2_sample[k]/nu))
      next
    }
    
    # Case 2: Non-empty cluster
    x_k <- x[idx_k]
    x_k[is.na(x_k) | is.infinite(x_k)] <- 0
    
    # Sample statistics
    x_bar <- mean(x_k)
    sum_xk <- sum(x_k)
    sse <- sum((x_k - x_bar)^2)
    
    # Posterior update for Sigma2 (marginalizing Mu)
    shape_sigma <- a + N_k / 2
    rate_sigma <- b + sse / 2 + (nu * N_k * (x_bar - mu_0)^2) / (2 * (N_k + nu))
    
    shape_sigma <- max(shape_sigma, 1e-6)
    rate_sigma <- max(rate_sigma, 1e-6)
    
    Sigma2_sample[k] <- 1 / rgamma(1, shape = shape_sigma, rate = rate_sigma)
    Sigma2_sample[k] <- pmin(pmax(Sigma2_sample[k], 1e-8), 1e8)
    
    # Posterior update for Mu
    mean_mu <- (sum_xk + nu * mu_0) / (N_k + nu)
    var_mu <- Sigma2_sample[k] / (N_k + nu)
    
    var_mu <- max(var_mu, 1e-8)
    sd_mu <- sqrt(var_mu)
    
    Mu_sample[k] <- rnorm(1, mean = mean_mu, sd = sd_mu)
  }
  
  return(list(Mu = Mu_sample, Sigma2 = Sigma2_sample))
}


# =================================================================
# Update mixture weights (empirical proportions)
# =================================================================
update_Pi_inf <- function(z, K) {
  counts <- tabulate(z, nbins = K)
  
  if (length(counts) != K) {
    stop("Cluster labels in z must be integers from 1 to K.")
  }
  
  Pi_sample <- counts / sum(counts)
  
  return(Pi_sample)
}


# =================================================================
# Escobar & West (1995) Gibbs update for alpha (DPM-BART)
# =================================================================
update_alphagibbs <- function(K, a_alpha, b_alpha, alpha, tau, x) { 
  n <- length(x)
  if (n == 0) stop("Sample size n cannot be zero.")
  
  if (length(tau) > 1) {
    log_tau <- sum(log(tau))  
  } else {
    log_tau <- log(tau)       
  }
  
  denominator <- n * (b_alpha - log_tau)
  if (denominator <= 0) denominator <- 1e-8  
  
  numerator <- (a_alpha + K - 1)
  prob1 <- numerator / (numerator + denominator)
  
  comp <- sample(c(1, 2), size = 1, prob = c(prob1, 1 - prob1))
  
  rate_param <- b_alpha - log_tau
  if (rate_param <= 0) rate_param <- 1e-8  
  
  if (comp == 1) {
    alpha_new <- rgamma(1, shape = a_alpha + K, rate = rate_param)
  } else {
    alpha_new <- rgamma(1, shape = a_alpha + K - 1, rate = rate_param)
  }
  alpha_new <- max(alpha_new, 1e-8)
  
  tau_new <- rbeta(n = 1, shape1 = alpha + 1, shape2 = n)
  
  return(list(alpha_new = alpha_new, tau_new = tau_new))
}


# =================================================================
# Update Dirichlet weights for variable selection probabilities
# =================================================================
update_s_  <- function(var_count, p, alpha_s) {
  shape <- alpha_s / p + var_count
  temp  <- rgamma(length(shape), shape, rate = 1)
  temp / sum(temp)
}
