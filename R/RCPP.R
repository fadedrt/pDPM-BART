library(Rcpp)

# 直接在 R 脚本中编译共用的 C++ 加速核心
cppFunction('
List compute_mixnormal_sums_cpp(NumericVector R, IntegerVector z, 
                                NumericVector Mu, NumericVector Sigma2, 
                                IntegerVector node_indices, IntegerVector which_terminal) {
    int n = R.size();
    int n_terminals = which_terminal.size();
    
    // 寻找最大节点编号以分配内存 (通常很小，几十到几百)
    int max_node = 0;
    for(int i = 0; i < n; i++) {
        if(node_indices[i] > max_node) {
            max_node = node_indices[i];
        }
    }
    
    // 初始化求和数组，全填 0
    std::vector<double> s1(max_node + 1, 0.0);
    std::vector<double> s2(max_node + 1, 0.0);
    std::vector<double> s3(max_node + 1, 0.0);
    
    // 单次循环完成所有观测值的分配和计算 (O(N) 极速遍历)
    for(int i = 0; i < n; i++) {
        int cluster = z[i] - 1; // R (1-based) 转换为 C++ (0-based)
        double mu_i = Mu[cluster];
        double sig2_i = Sigma2[cluster];
        
        double diff = R[i] - mu_i;
        double inv_2sig2 = 1.0 / (2.0 * sig2_i);
        
        int node = node_indices[i];
        
        s1[node] += diff * inv_2sig2;
        s2[node] += (diff * diff) * inv_2sig2;
        s3[node] += inv_2sig2;
    }
    
    // 提取叶节点对应的结果并返回
    NumericVector out_sum1(n_terminals);
    NumericVector out_sum2(n_terminals);
    NumericVector out_sum3(n_terminals);
    
    for(int k = 0; k < n_terminals; k++) {
        int node = which_terminal[k];
        if(node <= max_node) {
            out_sum1[k] = s1[node];
            out_sum2[k] = s2[node];
            out_sum3[k] = s3[node];
        } else {
            // 如果某个叶节点没有分到任何数据，求和自动为0
            out_sum1[k] = 0.0;
            out_sum2[k] = 0.0;
            out_sum3[k] = 0.0;
        }
    }
    
    return List::create(Named("sum1") = out_sum1,
                        Named("sum2") = out_sum2,
                        Named("sum3") = out_sum3);
}
')

tree_full_conditional_mixnormal <- function(tree, R, Mu, Sigma2, z, sigma2_mu, common_vars, aux_factor_var) {
  
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  
  if(nrow(tree$tree_matrix) != 1) {
    terminal_ancestors = get_ancestors(tree) 
    aux_table = table(terminal_ancestors[,1], terminal_ancestors[,2]) 
    which_terminal_no_double_split = NULL
    for (k in 1:nrow(aux_table)){
      split_var_ancestors = names(aux_table[k,])[aux_table[k,] != 0]
      invalid_ancestors = any(unlist(lapply(aux_factor_var, function(x) all(split_var_ancestors %in% x))))
      if (invalid_ancestors == TRUE) {which_terminal_no_double_split[k] = rownames(aux_table)[k]}
    }
    which_terminal_no_double_split = as.numeric(which_terminal_no_double_split[!is.na(which_terminal_no_double_split)])
  } else{
    which_terminal_no_double_split = 0 
  }
  
  sigma2_mu_aux = rep(sigma2_mu, length(which_terminal))
  sigma2_mu_aux[which(which_terminal %in% which_terminal_no_double_split)] = 1e-200
  
  # 【调用 C++ 加速端】
  sums <- compute_mixnormal_sums_cpp(R, z, Mu, Sigma2, tree$node_indices, which_terminal)
  
  sum1 <- sums$sum1
  sum2 <- sums$sum2
  sum3 <- sums$sum3
  
  log_post = -0.5 * sum(log(1 + 2 * sum3 * sigma2_mu_aux)) + sum((sum1^2)/((1/(2*sigma2_mu_aux)) + sum3) - sum2)
  
  return(log_post)
}



simulate_mu_mixnormal <- function(tree, R, Mu, Sigma2, z, sigma2_mu, common_vars, aux_factor_var) {
  
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  
  if(nrow(tree$tree_matrix) != 1) {
    terminal_ancestors = get_ancestors(tree) 
    aux_table = table(terminal_ancestors[,1], terminal_ancestors[,2]) 
    which_terminal_no_double_split = NULL
    for (k in 1:nrow(aux_table)){
      split_var_ancestors = names(aux_table[k,])[aux_table[k,] != 0]
      invalid_ancestors = any(unlist(lapply(aux_factor_var, function(x) all(split_var_ancestors %in% x))))
      if (invalid_ancestors == TRUE) {which_terminal_no_double_split[k] = rownames(aux_table)[k]}
    }
    which_terminal_no_double_split = as.numeric(which_terminal_no_double_split[!is.na(which_terminal_no_double_split)]) 
  } else{
    which_terminal_no_double_split = 0 
  }
  
  sigma2_mu_aux = rep(sigma2_mu, length(which_terminal))
  sigma2_mu_aux[which(which_terminal %in% which_terminal_no_double_split)] = 1e-200
  
  nj = as.numeric(tree$tree_matrix[which_terminal,'node_size'])
  
  # 【调用同样的 C++ 加速端】
  sums <- compute_mixnormal_sums_cpp(R, z, Mu, Sigma2, tree$node_indices, which_terminal)
  
  sum1 <- sums$sum1
  sum3 <- sums$sum3
  
  mu   <- rnorm(length(nj),
                mean = sum1/(sum3+1/(2*sigma2_mu_aux)),
                sd = sqrt(1/(2*sum3+1/sigma2_mu_aux)))
  
  tree$tree_matrix[,'mu'] = NA
  tree$tree_matrix[which_terminal,'mu'] = mu
  tree$tree_matrix[which_terminal_no_double_split, 'mu'] = 0 
  
  return(tree)
}


cppFunction('
List fill_tree_details_cpp(NumericMatrix tree_matrix, NumericMatrix X) {
    int n_obs = X.nrow();
    int n_nodes = tree_matrix.nrow();
    
    NumericMatrix new_tree_matrix = clone(tree_matrix);
    IntegerVector node_indices(n_obs, 1); 
    
    for(int i = 1; i < n_nodes; ++i) {
        int curr_parent = new_tree_matrix(i, 3) - 1; 
        int split_var = new_tree_matrix(curr_parent, 4) - 1; 
        double split_val = new_tree_matrix(curr_parent, 5);
        bool is_left = (new_tree_matrix(curr_parent, 1) - 1 == i);
        
        int node_size = 0;
        for(int obs = 0; obs < n_obs; ++obs) {
            if(node_indices[obs] == curr_parent + 1) {
                double x_val = X(obs, split_var);
                if(is_left) {
                    if(x_val < split_val) {
                        node_size++;
                        node_indices[obs] = i + 1;
                    }
                } else {
                    if(x_val >= split_val) {
                        node_size++;
                        node_indices[obs] = i + 1;
                    }
                }
            }
        }
        new_tree_matrix(i, 7) = node_size; 
    }
    
    return List::create(Named("tree_matrix") = new_tree_matrix,
                        Named("node_indices") = node_indices);
}
')

# 用包装函数替换原来的纯 R 循环
fill_tree_details = function(curr_tree, X) {
  # 调用刚才编译的 C++ 后端
  res <- fill_tree_details_cpp(curr_tree$tree_matrix, X)
  
  return(list(tree_matrix = res$tree_matrix,
              node_indices = res$node_indices))
}




# =================================================================
# 1. C++ 核心加速代码：随机扫描的留一法 Gibbs 抽样 (完全使用你提供的版本)
# =================================================================
cppFunction('
List update_z_inf_pCRP_pure_cpp(NumericVector x, int K, IntegerVector z, 
                                double alpha_current, double nu, double a, 
                                double b, double mu_0, double r) {
    int n = x.size();
    
    // 1. 复制状态并转换为 0-based
    std::vector<int> z_update(n);
    for(int i = 0; i < n; i++) {
        z_update[i] = z[i] - 1; 
    }
    
    std::vector<int> N_k(K, 0);
    std::vector<double> sum_x(K, 0.0);
    std::vector<double> sum_x2(K, 0.0);
    
    for(int i = 0; i < n; i++) {
        int cl = z_update[i];
        if(cl >= 0 && cl < K) {
            N_k[cl]++;
            sum_x[cl] += x[i];
            sum_x2[cl] += x[i] * x[i];
        }
    }
    
    double nu_mu0 = nu * mu_0;
    double nu_mu0_sq = nu * std::pow(mu_0, 2);
    
    // 【修复 1】：新族群的常数项，正确公式为 sqrt( nu / (1 + nu) )
    double sqrt_const_log = 0.5 * std::log(nu / (1.0 + nu));
    
    // 2. 调用 R 内置的 sample 保证打乱顺序的随机种子行为一致
    IntegerVector idx_seq = seq_len(n) - 1;
    IntegerVector shuffle_idx = Rcpp::sample(idx_seq, n, false);
    
    // 3. 开始遍历更新
    for(int idx = 0; idx < n; idx++) {
        int i = shuffle_idx[idx];
        int z_i_old = z_update[i];
        
        // 临时移除
        if (z_i_old >= 0 && z_i_old < K) {
            N_k[z_i_old]--;
            sum_x[z_i_old] -= x[i];
            sum_x2[z_i_old] -= x[i] * x[i];
        }
        
        double sum_N_r = 0.0;
        for(int k = 0; k < K; k++) {
            if(N_k[k] > 0) sum_N_r += std::pow(N_k[k], r);
        }
        double denominator = sum_N_r + alpha_current;
        
        std::vector<double> log_p(K + 1, -INFINITY);
        
        // 计算已有族群对数概率
        for(int k = 0; k < K; k++) {
            if(N_k[k] == 0) continue;
            
            double weight_k_log = std::log(std::pow(N_k[k], r) / denominator);
            
            int N_k_new = N_k[k] + 1;
            double sum_x_new = sum_x[k] + x[i];
            double sum_x2_new = sum_x2[k] + x[i] * x[i];
            
            // 【修复 2】：移除错误的 -1.0 和 -1.5，严格遵守 a + N/2
            double alpha_old = std::max(a + N_k[k] / 2.0, 1e-6);
            double alpha_new = std::max(a + N_k_new / 2.0, 1e-6);
            
            double term_beta_old = sum_x2[k] + nu_mu0_sq - std::pow(sum_x[k] + nu_mu0, 2) / (N_k[k] + nu);
            double beta_old = b + 0.5 * term_beta_old;
            
            double term_beta_new = sum_x2_new + nu_mu0_sq - std::pow(sum_x_new + nu_mu0, 2) / (N_k_new + nu);
            double beta_new = b + 0.5 * term_beta_new;
            
            // 【修复 3】：修复分子分母颠倒的问题 (N_old + nu) / (N_new + nu)
            double log_lik_k = -0.5 * std::log(2.0 * M_PI) + 
                               0.5 * std::log((N_k[k] + nu) / (N_k_new + nu)) +
                               std::lgamma(alpha_new) - std::lgamma(alpha_old) + 
                               alpha_old * std::log(beta_old) - alpha_new * std::log(beta_new);
                                  
            log_p[k] = weight_k_log + log_lik_k;
        }
        
        // 计算新族群对数概率
        double weight_new_log = std::log(alpha_current / denominator);
        
        // 【修复 4】：严格遵守 a + 1/2
        double alpha0 = std::max(a + 0.5, 1e-6);
        double term_beta0 = x[i]*x[i] + nu_mu0_sq - std::pow(x[i] + nu_mu0, 2) / (1.0 + nu);
        double beta0 = b + 0.5 * term_beta0;
        
        // 使用正确的 sqrt_const_log
        double log_lik_new = -0.5 * std::log(2.0 * M_PI) + sqrt_const_log +
                             std::lgamma(alpha0) - std::lgamma(a) + 
                             a * std::log(b) - alpha0 * std::log(beta0);
                                  
        log_p[K] = weight_new_log + log_lik_new;
        
        // 【防溢出终极武器：Log-Sum-Exp Trick】
        double max_log_p = -INFINITY;
        for(int k = 0; k <= K; k++) {
            if(log_p[k] > max_log_p) {
                max_log_p = log_p[k];
            }
        }
        
        std::vector<double> p_z(K + 1, 0.0);
        double sum_pz = 0.0;
        
        for(int k = 0; k <= K; k++) {
            if(log_p[k] != -INFINITY) {
                p_z[k] = std::exp(log_p[k] - max_log_p);
                sum_pz += p_z[k];
            }
        }
        
        // 归一化
        NumericVector p_z_norm(K + 1);
        if (sum_pz < 1e-32) {
            for(int k = 0; k <= K; k++) p_z_norm[k] = 1.0 / (K + 1.0);
        } else {
            for(int k = 0; k <= K; k++) p_z_norm[k] = p_z[k] / sum_pz;
        }
        
        // 抽样
        IntegerVector choices = seq_len(K + 1);
        IntegerVector sampled = Rcpp::sample(choices, 1, false, p_z_norm);
        int z_i_new = sampled[0] - 1; 
        
        z_update[i] = z_i_new;
        
        // 状态分配
        if (z_i_new == K) {
            N_k.push_back(1);
            sum_x.push_back(x[i]);
            sum_x2.push_back(x[i] * x[i]);
            K++;
        } else {
            N_k[z_i_new]++;
            sum_x[z_i_new] += x[i];
            sum_x2[z_i_new] += x[i] * x[i];
        }
    }
    
    // 4. 返回 1-based 的结果给 R
    IntegerVector z_out(n);
    for(int i = 0; i < n; i++) {
        z_out[i] = z_update[i] + 1;
    }
    
    return List::create(Named("z_update") = z_out,
                        Named("final_K") = K);
}
')

update_z_inf_pCRP <- function(x, K, z, alpha_current, nu, a, b, mu_0, r) {
  # 彻底移除了溢出惩罚限制，完全依赖纯粹的贝叶斯先验约束
  res <- update_z_inf_pCRP_pure_cpp(
    x = x, K = K, z = z, alpha_current = alpha_current, 
    nu = nu, a = a, b = b, mu_0 = mu_0, r = r
  )
  return(list(z_update = res$z_update, final_K = res$final_K))
}

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






