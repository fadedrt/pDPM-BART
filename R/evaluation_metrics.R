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


calculate_waic_pdpmbart <- function(fit, y) {
  # 1. 安全提取基础维度
  # 直接从存储的矩阵或列表长度获取，不依赖额外的 meta info
  npost <- length(fit$Mu_store)
  n <- length(y)
  
  # 防御性检查
  if (npost == 0 || is.null(npost)) {
    stop("错误：无法从 fit 对象中找到后验样本，请检查 fit$Mu_store 是否为空。")
  }
  
  # 2. 初始化对数似然矩阵 (npost 行, n 列)
  log_lik_matrix <- matrix(NA, nrow = npost, ncol = n)
  
  cat("正在计算 WAIC (极速安全版)...\n")
  
  # 3. 遍历每次 MCMC 抽样，计算逐点对数似然
  for (s in 1:npost) {
    # 提取第 s 次迭代的参数 (这些参数在主函数中已被还原到 y 的原始尺度)
    mu_s    <- fit$Mu_store[[s]]
    sigma_s <- sqrt(fit$Sigma2_store[[s]])
    pi_s    <- fit$Pi_store[[s]]
    K_s     <- fit$K_store[s]
    f_s     <- fit$bart_store[s, ]
    
    # 计算当前迭代下的残差 (向量化，长度为 n)
    resids <- y - f_s
    
    # 构建临时矩阵，存放每个样本在每个组分下的对数似然
    log_mix_mat <- matrix(-Inf, nrow = n, ncol = K_s)
    
    for (k in 1:K_s) {
      if (pi_s[k] > 0) {
        # 向量化计算 dnorm，一次性算完 n 个样本在第 k 个簇下的似然
        log_mix_mat[, k] <- log(pi_s[k]) + dnorm(resids, mean = mu_s[k], sd = sigma_s[k], log = TRUE)
      }
    }
    
    # Log-Sum-Exp (LSE) 技巧整合各组分的密度，防止指数运算下溢
    max_vals <- apply(log_mix_mat, 1, max)
    log_lik_matrix[s, ] <- max_vals + log(rowSums(exp(log_mix_mat - max_vals)))
  }
  
  # 4. 计算 WAIC 的核心组件
  
  # (1) lppd: 对数逐点预测密度的和 (拟合优度，越大越好)
  # 使用 LSE 技巧在对数尺度上求算术平均：log(mean(exp(x))) = log(sum(exp(x))) - log(S)
  lppd_i <- apply(log_lik_matrix, 2, function(x) {
    max_x <- max(x)
    max_x + log(mean(exp(x - max_x)))
  })
  lppd <- sum(lppd_i)
  
  # (2) p_waic: 对数似然的后验方差和 (复杂度惩罚，越小越好)
  p_waic_i <- apply(log_lik_matrix, 2, var)
  p_waic <- sum(p_waic_i)
  
  # (3) 最终 WAIC (越小越好)
  waic <- -2 * (lppd - p_waic)
  
  cat("--------------------------\n")
  cat("lppd (拟合优度, 越大越好):   ", round(lppd, 3), "\n")
  cat("p_waic (惩罚项, 越小越好):   ", round(p_waic, 3), "\n")
  cat("WAIC 总值 (越小越好):        ", round(waic, 3), "\n")
  cat("--------------------------\n")
  
  return(list(
    WAIC   = waic, 
    lppd   = lppd, 
    p_waic = p_waic
  ))
}

calculate_lpml_pdpmbart <- function(y, f_hat_train, Mu_store, Sigma2_store, Pi_store) {
  npost <- nrow(f_hat_train)
  n <- length(y)
  log_cpo <- numeric(n)
  
  for (i in 1:n) {
    log_pred <- numeric(npost)
    
    for (s in 1:npost) {
      mu_resid <- y[i] - f_hat_train[s, i]
      mu_s     <- Mu_store[[s]]
      sigma_s  <- sqrt(Sigma2_store[[s]])
      pi_s     <- Pi_store[[s]]
      
      # 计算当前样本在各簇下的 log-density
      log_mix_comps <- log(pi_s) + dnorm(mu_resid, mean = mu_s, sd = sigma_s, log = TRUE)
      
      # Log-Sum-Exp 技巧整合混合分布密度
      max_val <- max(log_mix_comps)
      log_pred[s] <- max_val + log(sum(exp(log_mix_comps - max_val)))
    }
    
    # 标准 CPO 调和平均估计 (Harmonic Mean Estimator)
    u <- -log_pred
    max_u <- max(u)
    log_cpo[i] <- -(max_u + log(mean(exp(u - max_u))))
  }
  
  return(sum(log_cpo))
}
