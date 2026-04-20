# =================================================================
# 辅助函数：Escobar & West (1995) 经典的 alpha Gibbs 抽样
# =================================================================
update_alphagibbs <- function(K, a_alpha, b_alpha, alpha, tau, x) { 
  n <- length(x)
  if (n == 0) stop("x不能为空，样本量n不能为0！")
  
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
# 主函数：标准 DPM-BART 基线模型 (r = 1 固定)
# =================================================================
bart_mixnormal_dp <- function(
    y,
    x,
    sparse           = FALSE,
    ntrees           = 50,
    node_min_size    = 5,
    alpha_tree       = 0.95,
    beta_tree        = 2,
    init_K           = 3, 
    fraction         = 0.3,  
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
  
  cat("开始标准 DPM-BART (r=1) 抽样...\n")
  cat("已启用 Escobar & West (1995) 数据增强法更新 alpha\n")
  cat("期望簇内方差占比 (fraction):", fraction * 100, "%\n\n")
  
  # --- 存储空间 ---
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
  
  y_mean <- mean(y)
  y_sd   <- sd(y)
  if (y_sd < 1e-8) y_sd <- 1e-6
  y_scale <- (y - y_mean) / y_sd
  
  mu_0_std <- 0  
  b_std    <- fraction * 1.0 * (a - 1) 
  s <- rep(1 / p, p)
  
  # --- 初始化 ---
  alpha_current <- 1 / init_K         
  z             <- sample(1:init_K, size = n, replace = TRUE)  
  K             <- init_K                                      
  
  # 【修改 1】：初始化辅助变量 tau
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
  
  # --- 主循环 ---
  for (i in seq_len(total_iterations)) {
    utils::setTxtProgressBar(pb, i)
    
    # 1. BART Tree 更新
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
    
    # 2. DP 更新 (强制 r = 1.0)
    z_update_res <- update_z_inf_pCRP(x = current_residuals, K = K, z = z, alpha_current = alpha_current, 
                                      nu = nu, a = a, b = b_std, mu_0 = mu_0_std, r = 1.0)
    z <- z_update_res$z_update
    K <- z_update_res$final_K
    
    z_counts <- table(z)
    non_empty_labels <- as.integer(names(z_counts))
    if (length(non_empty_labels) < K) {
      z <- match(z, non_empty_labels)
      K <- length(non_empty_labels)
    }
    
    sigmamu_res <- update_sigmamu(x = current_residuals, z = z, Mu = Mu, Sigma2 = Sigma2, 
                                  mu_0 = mu_0_std, nu = nu, a = a, b = b_std, K = K)
    Mu     <- sigmamu_res$Mu
    Sigma2 <- sigmamu_res$Sigma2
    
    # 【更新处】：改为直接基于经验比例计算 Pi (去除了 alpha_current 参数)
    Pi     <- update_Pi_inf(z = z, K = K)
    
    # 【修改 2】：使用经典的 Escobar & West (1995) 方法更新 alpha
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
    
    # 3. 变量选择
    if (isTRUE(sparse) && i > floor(total_iterations * 0.1)) {
      s <- update_s(var_count, p, 1)
    }
    
    # 4. 后验收集
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
  
  cat("\n标准 DPM-BART 抽样完成！\n")
  
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





