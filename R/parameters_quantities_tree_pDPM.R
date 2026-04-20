# =================================================================
# 2. R 层面的严谨先验与 MH 抽样（基于精确联合似然与潜排列 pi）
# =================================================================
# 之前被误删的 pCRP 伪似然版本 (仅计算分子)
# 之前被误删的 pCRP 伪似然版本 (仅计算分子)
# 完美修正版：基于生成序列的精确 pCRP 对数似然评估
# =================================================================
# 核心引擎：计算完全序列伪似然 (带分母惩罚，拒绝 r 无限膨胀)
# =================================================================
# =================================================================
# 核心引擎：计算完全序列伪似然 (带分母惩罚，拒绝 r 无限膨胀)
# =================================================================
# =================================================================
# 修正版 1：使用“留一法”伪似然计算 Z 的对数概率
# 这样能保证 R 层的参数更新 (alpha, r) 与 C++ 层的 Z 更新逻辑完全统一
# =================================================================
log_pcrp_pseudo_lik <- function(z, r, alpha) {
  n <- length(z)
  log_prob <- 0
  
  # 获取最终各簇人数
  n_k <- rep(0, max(z))
  tbl <- table(z)
  n_k[as.integer(names(tbl))] <- as.integer(tbl)
  
  # 对每个人进行“留一法”评估
  for (i in 1:n) {
    cluster_i <- z[i]
    n_k[cluster_i] <- n_k[cluster_i] - 1 # 排除自己
    
    sum_n_r <- sum(n_k[n_k > 0]^r)
    denominator <- alpha + sum_n_r
    
    if (n_k[cluster_i] == 0) {
      numerator <- alpha          
    } else {
      numerator <- n_k[cluster_i]^r 
    }
    
    log_prob <- log_prob + log(numerator) - log(denominator)
    n_k[cluster_i] <- n_k[cluster_i] + 1 # 恢复
  }
  return(log_prob)
}
# 辅助函数 1：严格计算 pCRP 下聚类划分 Z 的对数先验概率 (依赖历史入座轨迹)
# =================================================================
# 修正版 1：使用“留一法”伪似然计算 Z 的对数概率
# 这样能保证 R 层的参数更新 (alpha, r) 与 C++ 层的 Z 更新逻辑完全统一
# =================================================================
# 辅助函数 2：使用 Metropolis-Hastings 更新 alpha (纯伪似然，无 pi)
sample_alpha_mh <- function(alpha_current, z, r, a_alpha, b_alpha, mh_proposal_sd = 0.5) {
  log_alpha_current <- log(alpha_current)
  log_alpha_prop <- rnorm(1, mean = log_alpha_current, sd = mh_proposal_sd)
  alpha_prop <- exp(log_alpha_prop)
  
  log_prior_current <- dgamma(alpha_current, shape = a_alpha, rate = b_alpha, log = TRUE)
  log_prior_prop    <- dgamma(alpha_prop,    shape = a_alpha, rate = b_alpha, log = TRUE)
  
  # 直接在当前划分 z 上计算伪似然
  log_lik_current <- log_pcrp_pseudo_lik(z, r, alpha_current)
  log_lik_prop    <- log_pcrp_pseudo_lik(z, r, alpha_prop)
  
  log_jacobian <- log_alpha_prop - log_alpha_current
  log_acc <- (log_lik_prop + log_prior_prop + log_jacobian) - (log_lik_current + log_prior_current)
  
  if (!is.na(log_acc) && log(runif(1)) < log_acc) return(alpha_prop)
  return(alpha_current)
}

# 辅助函数 3：使用 Metropolis-Hastings 更新幂次 r (半正态先验，纯伪似然)
sample_r_mh <- function(r_current, z, alpha, r_prior_sigma, mh_proposal_sd) {
  
  # 1. 对称随机游走提议新的 r 值
  r_prop <- rnorm(1, mean = r_current, sd = mh_proposal_sd)
  
  # 2. 物理硬墙：强制 r > 1。如果提议的值 <= 1，直接拒绝（因为先验概率为 0）
  if (r_prop <= 1) return(r_current)
  
  # 3. 计算对数先验概率 (Shifted Half-Normal: r - 1 ~ Normal(0, sigma))
  # 常数项在 MH 比例中会自动抵消，因此直接用 dnorm 是完全数学等价的
  log_prior_current <- dnorm(r_current - 1, mean = 0, sd = r_prior_sigma, log = TRUE)
  log_prior_prop    <- dnorm(r_prop - 1,    mean = 0, sd = r_prior_sigma, log = TRUE)
  
  # 4. 直接在当前划分 z 上计算伪似然
  log_lik_current <- log_pcrp_pseudo_lik(z, r_current, alpha)
  log_lik_prop    <- log_pcrp_pseudo_lik(z, r_prop, alpha)
  
  # 5. 计算 Metropolis-Hastings 接受率
  log_acc <- (log_lik_prop + log_prior_prop) - (log_lik_current + log_prior_current)
  
  # 6. 接受 / 拒绝 判定
  if (!is.na(log_acc) && log(runif(1)) < log_acc) return(r_prop)
  
  return(r_current)
}





update_sigmamu <- function(x, z, Mu, Sigma2, mu_0, nu, a, b, K) {
  Mu_sample <- numeric(K)
  Sigma2_sample <- numeric(K)
  
  # 先验参数合法性检查（兜底）
  a <- max(a, 1e-6)  
  b <- max(b, 1e-6)  
  nu <- max(nu, 1e-6) 
  
  for (k in 1:K) {
    idx_k <- which(z == k)
    N_k <- length(idx_k)  
    
    # 1. 处理空族群（N_k=0）：使用先验兜底
    if (N_k == 0) {
      Sigma2_sample[k] <- 1 / rgamma(1, shape = a, rate = b)
      Sigma2_sample[k] <- pmin(pmax(Sigma2_sample[k], 1e-8), 1e8) # 防御
      Mu_sample[k] <- rnorm(1, mean = mu_0, sd = sqrt(Sigma2_sample[k]/nu))
      next
    }
    
    # 2. 非空族群：数值稳定性处理
    x_k <- x[idx_k]
    x_k[is.na(x_k) | is.infinite(x_k)] <- 0
    
    # ------------------ 严格的数学修正部分 ------------------
    # 计算样本统计量
    x_bar <- mean(x_k)
    sum_xk <- sum(x_k)
    sse <- sum((x_k - x_bar)^2) # 基于样本均值的组内平方和
    
    # 3. 更新 Sigma2（边际后验：已积分掉 Mu）
    shape_sigma <- a + N_k / 2
    # 严格遵循 NIG 的 rate 后验公式：b + SSE/2 + 先验惩罚项
    rate_sigma <- b + sse / 2 + (nu * N_k * (x_bar - mu_0)^2) / (2 * (N_k + nu))
    
    # 兜底：确保shape和rate>0
    shape_sigma <- max(shape_sigma, 1e-6)
    rate_sigma <- max(rate_sigma, 1e-6)
    
    Sigma2_sample[k] <- 1 / rgamma(1, shape = shape_sigma, rate = rate_sigma)
    # 限制Sigma2范围，避免极端值导致后续Mu的方差爆炸
    Sigma2_sample[k] <- pmin(pmax(Sigma2_sample[k], 1e-8), 1e8)
    
    # 4. 更新 Mu（条件后验：依赖于刚刚抽出的新 Sigma2）
    mean_mu <- (sum_xk + nu * mu_0) / (N_k + nu)
    var_mu <- Sigma2_sample[k] / (N_k + nu)
    
    # 兜底：确保方差>0
    var_mu <- max(var_mu, 1e-8)
    sd_mu <- sqrt(var_mu)
    
    Mu_sample[k] <- rnorm(1, mean = mean_mu, sd = sd_mu)
    # ------------------------------------------------------
  }
  
  return(list(Mu = Mu_sample, Sigma2 = Sigma2_sample))
}
update_Pi_inf <- function(z, K) {
  # 1. 统计每个族群的数量（确保标签在1~K范围内）
  counts <- tabulate(z, nbins = K)
  
  # 验证计数长度是否为K（避免z中标签超出1~K导致的错误）
  if (length(counts) != K) {
    stop("z中的标签必须是1到", K, "之间的整数！")
  }
  
  # 2. 直接计算经验比例
  # 因为此时的 K 已经过滤掉了空簇，所以 sum(counts) 等于总样本数 n
  Pi_sample <- counts / sum(counts)
  
  return(Pi_sample)
}

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


update_s_  <- function(var_count, p, alpha_s) {
  shape <- alpha_s / p + var_count
  temp  <- rgamma(length(shape), shape, rate = 1)
  temp / sum(temp)
}
