
dpmbart <- function(
    y,
    x,
    sparse           = FALSE,
    ntrees           = 50,
    node_min_size    = 5,
    alpha_tree       = 0.95,
    beta_tree        = 2,
    init_K           = 3, 
    a                = 3,    
    nu               = 0.1,  
    a_alpha          = 1,    
    b_alpha          = 2,    
    nburn            = 1000,
    npost            = 500,
    nthin            = 1
) {
  
  sigma2_mu <- 1 / (4 * ntrees)
  aux_identify_factor_variables <- NULL
  common_variables              <- NULL
  
  total_iterations <- nburn + npost * nthin
  n <- length(y)
  p <- ncol(x)
  
  cat("Starting standard DPM-BART (r = 1) sampling...\n")
  cat("Using Escobar & West (1995) data augmentation for alpha update\n")
  cat("Prior calibrated via OLS residual variance to match P(sigma < sigma_hat) = 0.90\n\n")
  
  # --- Storage containers ---
  tree_store      <- vector('list', npost)
  y_hat_store     <- matrix(NA, ncol = n, nrow = npost)
  bart_store      <- matrix(NA, ncol = n, nrow = npost)
  var_count       <- rep(0, p)
  var_count_store <- matrix(0, ncol = p, nrow = npost)
  s_prob_store    <- matrix(0, ncol = p, nrow = npost)
  
  Sigma2_store    <- list()
  z_store         <- matrix(NA, ncol = n, nrow = npost)
  Mu_store        <- list()
  alpha_store     <- rep(NA, npost)
  Pi_store        <- list()
  K_store         <- rep(NA, npost)
  
  tree_fits_cache <- matrix(0, ncol = ntrees, nrow = n)
  
  # --- Standardize response ---
  y_mean <- mean(y)
  y_sd   <- sd(y)
  if (y_sd < 1e-8) y_sd <- 1e-6
  y_scale <- (y - y_mean) / y_sd
  
  # =====================================================================
  # Data-Driven Prior Calibration
  # =====================================================================
  lm_fit    <- lm(y_scale ~ as.matrix(x))
  sigma_hat <- sd(lm_fit$residuals)
  q_01      <- qchisq(0.1, df = 2 * a)
  b_std     <- (q_01 * sigma_hat^2) / 2
  # =====================================================================
  
  mu_0_std <- 0  
  s <- rep(1 / p, p)
  
  # --- Initialization ---
  alpha_current <- 1 / init_K         
  z             <- sample(1:init_K, size = n, replace = TRUE)  
  K             <- init_K                                      
  
  # Initialize auxiliary variable tau (Escobar & West scheme)
  tau_current   <- rbeta(1, shape1 = alpha_current + 1, shape2 = n) 
  
  residual_init           <- rep(0, n)  
  residuals_by_group_init <- split(residual_init, z)
  
  Mu <- sapply(residuals_by_group_init, function(grp) {
    (sum(grp) + nu * mu_0_std) / (length(grp) + nu)
  })
  
  Sigma2 <- sapply(seq_along(residuals_by_group_init), function(k) {
    grp <- residuals_by_group_init[[k]]
    n_k <- length(grp)
    sse <- sum((grp - mean(grp))^2)
    (2 * b_std + sse + (nu * n_k * (mean(grp) - mu_0_std)^2) / (nu + n_k)) / (2 * (a + n_k / 2))
  })
  
  curr_trees <- create_stump(num_trees = ntrees, y = y_scale, X = x)
  new_trees  <- curr_trees
  yhat_bart  <- get_predictions(curr_trees, x, single_tree = (ntrees == 1))
  
  pb <- utils::txtProgressBar(min = 1, max = total_iterations, style = 3, width = 60, title = 'Running DPM-BART...')
  
  # --- Main MCMC loop ---
  for (i in seq_len(total_iterations)) {
    utils::setTxtProgressBar(pb, i)
    
    # 1. BART tree updates
    for (j in seq_len(ntrees)) {
      current_partial_residuals <- y_scale - yhat_bart + tree_fits_cache[, j]
      move_type <- sample_move(curr_trees[[j]], i, nburn)
      
      new_trees[[j]] <- update_tree(
        y = current_partial_residuals, X = x, type = move_type, 
        curr_tree = curr_trees[[j]], node_min_size = node_min_size, 
        s = s, common_vars = common_variables, aux_factor_var = aux_identify_factor_variables
      )
      
      l_old <- tree_full_conditional_mixnormal(tree = curr_trees[[j]], R = current_partial_residuals, Mu = Mu, Sigma2 = Sigma2, z = z, sigma2_mu = sigma2_mu, common_vars = common_variables, aux_factor_var = aux_identify_factor_variables) + get_tree_prior(curr_trees[[j]], alpha_tree, beta_tree, common_variables)
      l_new <- tree_full_conditional_mixnormal(tree = new_trees[[j]], R = current_partial_residuals, Mu = Mu, Sigma2 = Sigma2, z = z, sigma2_mu = sigma2_mu, common_vars = common_variables, aux_factor_var = aux_identify_factor_variables) + get_tree_prior(new_trees[[j]], alpha_tree, beta_tree, common_variables)
      
      ratio <- if (isTRUE(new_trees[[j]]$ForceStump)) 1 else l_new - l_old
      if (ratio > 0 || ratio > log(runif(1))) {
        curr_trees[[j]] <- new_trees[[j]]
        if (move_type == 'change') {
          var_count[curr_trees[[j]]$var[1]] <- var_count[curr_trees[[j]]$var[1]] + 1
          var_count[curr_trees[[j]]$var[2]] <- var_count[curr_trees[[j]]$var[2]] - 1
        }
        if (move_type == 'grow')  var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] + 1 
        if (move_type == 'prune') var_count[curr_trees[[j]]$var] <- var_count[curr_trees[[j]]$var] - 1 
      }
      
      curr_trees[[j]] <- simulate_mu_mixnormal(tree = curr_trees[[j]], R = current_partial_residuals, Mu = Mu, Sigma2 = Sigma2, z = z, sigma2_mu = sigma2_mu, common_vars = common_variables, aux_factor_var = aux_identify_factor_variables)
      current_fit <- get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart   <- yhat_bart - tree_fits_cache[, j] + current_fit
      tree_fits_cache[, j] <- current_fit
    }
    
    y_hat <- yhat_bart
    current_residuals <- y_scale - y_hat  
    
    # 2. Dirichlet Process update (r explicitly fixed at 1.0)
    z_update_res <- update_z_inf_pCRP(x = current_residuals, K = K, z = z, alpha_current = alpha_current, 
                                      nu = nu, a = a, b = b_std, mu_0 = mu_0_std, r = 1.0)
    z <- z_update_res$z_update
    K <- z_update_res$final_K
    
    # Speed optimization: Use C-level tabulate to clean up empty clusters instantly
    counts <- tabulate(z)
    non_empty <- which(counts > 0)
    if (length(non_empty) < K) {
      mapping <- integer(max(z))
      mapping[non_empty] <- seq_along(non_empty)
      z <- mapping[z]
      K <- length(non_empty)
    }
    
    # Update cluster-specific parameters
    sigmamu_res <- update_sigmamu(x = current_residuals, z = z, Mu = Mu, Sigma2 = Sigma2, 
                                  mu_0 = mu_0_std, nu = nu, a = a, b = b_std, K = K)
    Mu     <- sigmamu_res$Mu
    Sigma2 <- sigmamu_res$Sigma2
    
    # Compute mixture weights via empirical proportions
    Pi     <- update_Pi_inf(z = z, K = K)
    
    # Update alpha using Escobar & West (1995) method
    alphatau_res <- update_alphagibbs(
      K       = K,
      a_alpha = a_alpha,
      b_alpha = b_alpha,
      alpha   = alpha_current,
      tau     = tau_current,
      x       = current_residuals 
    )
    alpha_current <- alphatau_res$alpha_new
    tau_current   <- alphatau_res$tau_new
    
    # 3. Variable selection (optional sparsity)
    if (isTRUE(sparse) && i > floor(total_iterations * 0.1)) {
      s <- update_s(var_count, p, 1)
    }
    
    # 4. Posterior storage
    if ((i > nburn) && ((i - nburn) %% nthin == 0)) {
      idx <- (i - nburn) / nthin
      
      tree_store[[idx]]      <- curr_trees
      y_hat_store[idx, ]     <- y_hat * y_sd + y_mean  
      bart_store[idx, ]      <- yhat_bart * y_sd + y_mean  
      var_count_store[idx, ] <- var_count
      s_prob_store[idx, ]    <- s
      
      alpha_store[idx]       <- alpha_current
      Mu_store[[idx]]        <- Mu * y_sd  
      Sigma2_store[[idx]]    <- Sigma2 * (y_sd^2)  
      z_store[idx, ]         <- z
      Pi_store[[idx]]        <- Pi
      K_store[idx]           <- K
    }
  } 
  
  close(pb) 
  
  cat("\nStandard DPM-BART sampling completed!\n")
  
  results <- list(
    trees            = tree_store,
    y_hat            = y_hat_store, 
    bart_store       = bart_store,
    var_count_store  = var_count_store,
    s_prob_store     = s_prob_store,
    alpha_store      = alpha_store,
    Mu_store         = Mu_store,
    Sigma2_store     = Sigma2_store,
    z_store          = z_store,
    Pi_store         = Pi_store,
    K_store          = K_store,
    model_params     = list(
      ntrees         = ntrees, 
      init_K         = init_K, 
      nburn          = nburn, 
      npost          = npost, 
      nthin          = nthin, 
      y_mean         = y_mean, 
      y_sd           = y_sd
    )
  )
  
  class(results) <- "bart_mixnormal_dp"
  return(results)
}

