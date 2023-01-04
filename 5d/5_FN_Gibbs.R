###############################################################
# Gibbs Sampling Implementation on 5-dimensional Half Normal Distribution
# original author: Jin
###############################################################
library(MASS)
library(clusterGeneration)
library(base)
library(dplyr)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(purrr)
library(tidyr)
library(reshape2)
library(latex2exp)
library(ggpubr)
library(hexbin)
###############################################################
## indicator function of random vector X
###############################################################
Indicator <- function(x){
  return(ifelse(FALSE %in% ifelse(x >= 0, TRUE, FALSE), 0, 1))
}
###############################################################
## symbolic function
###############################################################
s_f <- function(q, n){
  sign <- matrix(0, n, length(q))
  for(i in 1:n){
    sign[i,] <- sample(c(-1,1), length(q), replace = TRUE)
  }
  sign <- unique(sign)
  return(sign)
}
###############################################################
## symbolic matrix
###############################################################
symbolic_as_diag <- function(vector){
  return(diag(vector))
}
###############################################################
## Folded Normal Sampling Funcion
###############################################################
FN_Sampling <- function(u, sig){
  sample <- abs(rnorm(1, mean = u, sd = sqrt(sig)))
  return(sample)
}
###############################################################
## Mixing Folded Normal Sampling Funcion
###############################################################
Mixing_Sample <- function(N,d,com,pos,x,Sigma){
  ##################################################
  # 输入
  ##################################################
  ## 输入参数
  # N <- 10000 # 样本量
  # d <- 5 # 维度5d
  # com <- nrow(symbolic) # 组合2^5=32
  # pos # 剔除第几行第几列,目前pos = 1
  # x <- c(0.5, 0.3, 1, 2) # xj # 已剔除pos项的向量
  # Sigma # 协方差矩阵
  ##################################################
  # 相关矩阵
  ##################################################
  ## 协方差矩阵的det
  det <- det(Sigma)
  ## 协方差矩阵的伴随矩阵
  Sigma_ac <- det * solve(Sigma) 
  ## 剔除第一行第一列的协方差矩阵
  Sigma_e <- Sigma[-pos, -pos]
  det_e <- det(Sigma_e)
  ## 其伴随矩阵
  Sigma_e_ac <- det_e * solve(Sigma_e)
  ##################################################
  # 混合分布随机数
  ##################################################
  ## 混合分布的系数
  eps <- rep(0,com)
  for(k in 1:com){
    for(i in 1:(d-1)){
      for(j in 1:(d-1)){
        eps[k] <- eps[k] - 0.5 / det_e * symbolic[k,-pos][i] * x[i] * symbolic[k,-pos][j] * x[j] * Sigma_e_ac[i,j]
      }
    }
  }
  ## 方式一直接取exp，但容易出现nan值
  # normal_con <- exp(eps) / sum(exp(eps))
  # ## 方式二 求倒数
  eps_matrix <- matrix(rep(eps,com),byrow = T,com,com)
  for(i in 1:com){
    eps_matrix[i,] <- eps_matrix[i,] - eps_matrix[i,i]
  }
  eps_matrix <- exp(eps_matrix)
  normal_con <- 1/apply(eps_matrix, 1, sum)
  normal_con[is.na(normal_con)] <- 0
  # print(normal_con)
  ## 混合分布的参数
  u <- (1 / det_e) * symbolic[,-pos] %*% matrix(x * Sigma_ac[pos,][-pos], ncol = 1)
  sig <- det / det_e 
  ## 选择混合分布之一
  choice <- sample(c(1:nrow(symbolic)),1,prob = normal_con, replace = T)
  ## 生成随机数
  sample <- FN_Sampling(u[choice],sig)
  return(sample)
}
###############################################################
# Gibbs Sampling Function
###############################################################
Gibbs <- function(N,x_init,Sig){
  x_seq <- matrix(0, N, 5)
  x_seq1 <- x_seq2 <- x_seq3 <- x_seq4 <- rep(0,4)
  x_seq[1,] <- x_init
  for (i in 2:N){
    ## x1|
    x_seq1 <- x_seq[i - 1, -1]
    x_seq[i, 1] <- Mixing_Sample(N = 1, d = 5, com = nrow(symbolic),
                                 pos = 1, x = x_seq1, Sigma = Sig)
    ## x2
    x_seq2 <- c(x_seq[i, 1],x_seq1[-1])
    x_seq[i, 2] <- Mixing_Sample(N = 1, d = 5, com = nrow(symbolic),
                                 pos = 2, x = x_seq2, Sigma = Sig)
    ## x3
    x_seq3 <- c(x_seq[i, 1:2],x_seq1[-(1:2)])
    x_seq[i, 3] <- Mixing_Sample(N = 1, d = 5, com = nrow(symbolic),
                                 pos = 3, x = x_seq3, Sigma = Sig)
    ## x4
    x_seq4 <- c(x_seq[i, 1:3],x_seq1[-(1:3)])
    x_seq[i, 4] <- Mixing_Sample(N = 1, d = 5, com = nrow(symbolic),
                                 pos = 4, x = x_seq3, Sigma = Sig)
    ## x5
    x_seq5 <- x_seq[i, 1:4]
    x_seq[i, 5] <- Mixing_Sample(N = 1, d = 5, com = nrow(symbolic),
                                 pos = 5, x = x_seq4, Sigma = Sig)
    
  }
  return(x_seq)
}
###############################################################
## correlation
###############################################################
est_corr <- function(sig1, sig2, rho){
  intergrand <- function(x){
    x^2 * 2 * dnorm(x) * pnorm(x * rho/sqrt(1 - rho^2))
  }
  intergal <- integrate(intergrand, 0, Inf)$value - 0.5
  cov <- 2 * sqrt(sig1) * sqrt(sig2) * {(1 - rho^2) * sqrt(1 - rho^2)/pi + rho * intergal} - 2 * sqrt(sig1) * sqrt(sig2)/pi
  sig1 <- sig1 * (1 - 2/pi)
  sig2 <- sig2 * (1 - 2/pi)
  cor <- cov/sqrt(sig1 * sig2)
  return(cor)
}
# est_corr(sig1 = 1, sig2 = 1, rho = 0.95)
###############################################################
## input 
###############################################################
set.seed(1010)
eigen <- rexp(5,0.5)
symbolic <- s_f(1:5, 1000)
com <- nrow(symbolic)
Sigma <- genPositiveDefMat(dim = 5,
                           covMethod = "eigen",
                           eigenvalue = eigen)$Sigma
# R <- cor(data.frame(x1=rnorm(10,0,1),x2=rnorm(10,0,1),
#                     x3=rnorm(10,1.2),x4=rnorm(10,0,2),
#                     x5=rnorm(10,2,3))) # get correlations
# D <- diag(sqrt(runif(5,0,8)))
# Sigma <- D %*% R %*% D
N <- 10000
Chain_n <- 4
x_initial <- list(first = 0, second = 0, third = 0, fourth = 0)
x_init <- x_initial[[1]] <-matrix(c(0.5, 0.5, 1.5, 1.5, 1.0), ncol = 5)
x_initial[[2]] <- matrix(c(1.0, 1.5,0.5,1.0, 0.5), ncol = 5)
x_initial[[3]] <- matrix(c(1.5, 1.0, 1.0, 0.5, 1.5), ncol = 5)
x_initial[[4]] <- matrix(c(1.0, 0.5,1.0, 1.0, 1.5), ncol = 5)
###############################################################
## Sampling
###############################################################
set.seed(1010)
sample_p1 <- sample_p2 <- sample_p3 <- sample_p4 <- sample_p5 <- 0
Chain_all <- list(Chain1 = 0, Chain2 = 0, Chain3 = 0, Chain4 = 0)
t_start <- Sys.time()
for(j in 1:Chain_n){
  result <- Gibbs(N,x_init,Sigma)
  Chain_all[[j]] <- result
  if(j == 1){
    sample_p1 <- result[,1]
    sample_p2 <- result[,2]
    sample_p3 <- result[,3]
    sample_p4 <- result[,4]
    sample_p5 <- result[,5]
  }
  else{
    sample_p1 <- c(sample_p1, result[,1])
    sample_p2 <- c(sample_p2, result[,2])
    sample_p3 <- c(sample_p3, result[,3])
    sample_p4 <- c(sample_p4, result[,4])
    sample_p5 <- c(sample_p5, result[,5])
  }
  if(j != Chain_n){
    x_init = x_initial[[j+1]]
  }
}
t_end <- Sys.time()
t_end - t_start # 4.424277 mins 4.612946 mins
###############################################################
## data(in array)
###############################################################
dim2 <- paste0("Chain:", 1:4)
# dim3 <- paste0("theta", 1:5)
dim3 <- paste0("theta", "[",1:5,"]")
sample_array <- array(c(sample_p1,sample_p2,sample_p3,sample_p4,sample_p5),c(N,4,5),dimnames = list(NULL,dim2,dim3))
names(dimnames(sample_array)) <- c("Iteration", "Chain", "Parameter")
###############################################################
## ergodic mean plot
###############################################################
###############################################################
### ergodic data
################################
Chain_Para1_erg <-
  sample_array[,,1] %>%
  as.data.frame() %>%
  mutate(Iteration=seq_len(N)) %>%
  melt(id.vars = "Iteration", measure.vars = paste0("Chain:",1:4),variable.name = "Chain",value.name = "Para") %>%
  group_by(Chain) %>%
  mutate(ergodic_mean = cumsum(Para)/seq_len(N),
         Theta = "theta1") %>%
  as.data.frame()
Chain_Para2_erg <-
  sample_array[,,2] %>%
  as.data.frame() %>%
  mutate(Iteration=seq_len(N)) %>%
  melt(id.vars = "Iteration", measure.vars = paste0("Chain:",1:4),variable.name = "Chain",value.name = "Para") %>%
  group_by(Chain) %>%
  mutate(ergodic_mean = cumsum(Para)/seq_len(N),
         Theta = "theta2") %>%
  as.data.frame()
Chain_Para3_erg <-
  sample_array[,,3] %>%
  as.data.frame() %>%
  mutate(Iteration=seq_len(N)) %>%
  melt(id.vars = "Iteration", measure.vars = paste0("Chain:",1:4),variable.name = "Chain",value.name = "Para") %>%
  group_by(Chain) %>%
  mutate(ergodic_mean = cumsum(Para)/seq_len(N),
         Theta = "theta3") %>%
  as.data.frame()
Chain_Para4_erg <-
  sample_array[,,4] %>%
  as.data.frame() %>%
  mutate(Iteration=seq_len(N)) %>%
  melt(id.vars = "Iteration", measure.vars = paste0("Chain:",1:4),variable.name = "Chain",value.name = "Para") %>%
  group_by(Chain) %>%
  mutate(ergodic_mean = cumsum(Para)/seq_len(N),
         Theta = "theta4") %>%
  as.data.frame()
Chain_Para5_erg <-
  sample_array[,,5] %>%
  as.data.frame() %>%
  mutate(Iteration=seq_len(N)) %>%
  melt(id.vars = "Iteration", measure.vars = paste0("Chain:",1:4),variable.name = "Chain",value.name = "Para") %>%
  group_by(Chain) %>%
  mutate(ergodic_mean = cumsum(Para)/seq_len(N),
         Theta = "theta5") %>%
  as.data.frame()
Chain_erg <- 
  Chain_Para1_erg %>% 
  rbind(Chain_Para2_erg) %>%
  rbind(Chain_Para3_erg) %>%
  rbind(Chain_Para4_erg) %>%
  rbind(Chain_Para5_erg) %>%
  mutate(Theta_c = recode(Theta,
                          "theta1" = "theta[1]",
                          "theta2" = "theta[2]",
                          "theta3" = "theta[3]",
                          "theta4" = "theta[4]",
                          "theta5" = "theta[5]")
  )
###############################################################
### implement ergodic mean plot
################################
Chain_erg_Gibbs <-
  Chain_erg %>%
  group_by(Chain) %>%
  # top_n(n = -10000, wt = Iteration) %>%
  ggplot(aes(x=Iteration, y=ergodic_mean, colour=Chain, group=Chain)) +
  facet_grid(rows = vars(Theta_c), labeller = label_parsed) +
  labs(title = "吉布斯采样") +
  geom_line() +
  theme_minimal(base_size = 14) +
  # geom_hline(aes(yintercept = sqrt(2/pi)), linetype = 5, size = 0.3) +
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.key = element_rect(
          color = "white", 
          fill = "white"),
        strip.text.y = element_text(size = 15,angle = 0),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold")) +
  theme(axis.title.x = element_text(vjust = -0.1),
        axis.title.y = element_text(margin = margin(r = 15))) +
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(35, "pt")) +
  theme(legend.title.align = 0.5) +
  ylab("遍历均值") +
  xlab("迭代次数") +
  scale_color_discrete(labels = c(expression(paste(theta^(0)==symbol(J)[1])),
                                  expression(paste(theta^(0)==symbol(J)[2])),
                                  expression(paste(theta^(0)==symbol(J)[3])),
                                  expression(paste(theta^(0)==symbol(J)[4])))) 
###############################################################
## trace_plot
##############################################################
color_scheme_set("mix-blue-red")
traceplot_Gibbs <- 
mcmc_trace(sample_array, n_warmup = 2500,
           facet_args = list(ncol = 1, strip.position = "left", labeller = ggplot2::label_parsed)) +
  theme(axis.text.y = element_text(angle = 0),
        legend.position = "top",
        legend.spacing.x = unit(1.0, 'cm'),
        legend.key.size = unit(20, "pt"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  scale_color_discrete(name = "Chain", labels = c(expression(paste(theta^(0)==symbol(J)[1])),
                                  expression(paste(theta^(0)==symbol(J)[2])),
                                  expression(paste(theta^(0)==symbol(J)[3])),
                                  expression(paste(theta^(0)==symbol(J)[4])))) 
###############################################################
## PSRF: R^hat (top 2500)
###############################################################
top_n <- 2500
top_sub <- c(1:top_n, (N + 1):(N + top_n), (2 * N + 1):(2 * N + top_n),(3 * N + 1):(3 * N + top_n))
sample_p1_top <- sample_p1[top_sub]
sample_p2_top <- sample_p2[top_sub]
sample_p3_top <- sample_p3[top_sub]
sample_p4_top <- sample_p4[top_sub]
sample_p5_top <- sample_p5[top_sub]
dim2 <- paste0("chain:", 1:4)
# dim3 <- paste0("theta", 1:2)
dim3 <- c("theta[1]", "theta[2]","theta[3]", "theta[4]","theta[5]")
sample_array_top <- array(c(sample_p1_top,sample_p2_top,sample_p3_top,sample_p4_top,sample_p5_top),c(top_n,4,5),dimnames = list(NULL,dim2,dim3))
names(dimnames(sample_array_top)) <- c("Iteration", "Chain", "Parameter")
sample_cbind_para_top <- cbind(sample_array_top[,1,],sample_array_top[,2,],sample_array_top[,3,],sample_array_top[,4,])
sample_rbind_para_top <- rbind(sample_array_top[,1,],sample_array_top[,2,],sample_array_top[,3,],sample_array_top[,4,])
J = 4
L = top_n
col_index <- matrix(0,length(x_init),Chain_n)
for(j in 1:length(x_init)){
  col_index[j,] <- which(colnames(sample_cbind_para_top) == colnames(sample_rbind_para_top)[j])
}
for(j in 1:length(x_init)){
  x_j_mean <- apply(sample_cbind_para_top, 2, mean)[col_index[j,]]
  x_mean <- apply(sample_rbind_para_top, 2, mean)[j]
  B <- L/(J-1) * sum((x_j_mean - x_mean)^2)
  sj2 <- 1/(L-1) * apply((sample_cbind_para_top[,col_index[j,]] - x_j_mean)^2,2,sum)
  W <- 1/J * sum(sj2)
  R <- ((L-1)/L * W + 1/L * B + B/(L*J))/W
  print(R)
}
###############################################################
## (in array) after eliminating the Burn-in period data(The previous 5000)
###############################################################
Elim_n <- 2500
Elim_sub <- c(1:Elim_n, (N + 1):(N + Elim_n), (2 * N + 1):(2 * N + Elim_n),(3 * N + 1):(3 * N + Elim_n))
sample_p1_elim <- sample_p1[-Elim_sub]
sample_p2_elim <- sample_p2[-Elim_sub]
sample_p3_elim <- sample_p3[-Elim_sub]
sample_p4_elim <- sample_p4[-Elim_sub]
sample_p5_elim <- sample_p5[-Elim_sub]
sample_array_elim <- array(c(sample_p1_elim,sample_p2_elim,sample_p3_elim,sample_p4_elim,sample_p5_elim),
                           c(N-Elim_n,4,5),dimnames = list(NULL,dim2,dim3))
names(dimnames(sample_array_elim)) <- c("Iteration", "Chain", "Parameter")
###############################################################
## autocorrelation
###############################################################
color_scheme_set("green")
acf_Gibbs <- 
  mcmc_acf(sample_array_elim, lags = 20) + 
  hline_at(0.5, linetype = 2, size = 0.15, color = "gray") +
  xlab("滞后阶数") +
  ylab("自相关函数") +
  theme(axis.title.y = element_text(margin = margin(r = 10), size = 10),
        axis.title.x = element_text(vjust = -0.1, size = 10)) +
  facet_grid(Chain~Parameter,labeller = label_parsed)
###############################################################
## effective sample size
###############################################################
ess_Vari <- function(Chain) {
  C_N <- length(Chain)
  V <- map_dbl(seq_len(C_N - 1),
               function(t) {
                 mean(diff(Chain, lag = t) ^ 2, na.rm = TRUE)
               }
  )
  return(V)
}
ess <- function(Chain_all,Chain_n, para_n){
  V = 0
  para_all = c(0)[-1]
  for(i in 1:Chain_n){
    V = V + ess_Vari(Chain_all[[i]][(Elim_n+1):N,para_n])
    para_all = c(para_all,Chain_all[[i]][(Elim_n+1):N,para_n])
  }
  rho <- head_while(1 - V / 2 / var(para_all), ~ . > 0)
  (N-Elim_n)*Chain_n / (1 + sum(rho))
}
ess_para1 <- ess(Chain_all,4, 1) # the ess of 4 chains(theta1): 30000
ess_para2 <- ess(Chain_all,4, 2) # the ess of 4 chains(theta2): 30000
ess_para3 <- ess(Chain_all,4, 3) # the ess of 4 chains(theta3): 30000
ess_para4 <- ess(Chain_all,4, 4) # the ess of 4 chains(theta4): 30000
ess_para5 <- ess(Chain_all,4, 5) # the ess of 4 chains(theta5): 30000
###############################################################
## scatter
###############################################################
color_scheme_set("viridisA") 
scatter_density_Gibbs <- 
  mcmc_pairs(sample_array_elim, 
             diag_fun = "dens",
             off_diag_fun = "hex")
###############################################################
## maginal histogram
###############################################################
color_scheme_set("teal")
maginal_hist_Gibbs <- 
mcmc_hist_by_chain(sample_array_elim, regex_pars = "theta", binwidth = 0.25,
                   facet_args = list(labeller = ggplot2::label_parsed)) +
  facet_text(size = 15) 
###############################################################
## Compare MCMC estimates to "true" parameter values 
###############################################################
###############################################################
sample_rbind_para <- rbind(sample_array_elim[,1,],sample_array_elim[,2,],sample_array_elim[,3,],sample_array_elim[,4,])
## mean
est_mean <- apply(sample_rbind_para , 2, mean)
true_mean <- sqrt(diag(Sigma)) * sqrt(2/pi)
## var
est_var <- apply(sample_rbind_para, 2, var)
true_var <- diag(Sigma) * (1 - 2/pi)
## rho
D <- diag(sqrt(diag(Sigma)))
original_R <- solve(D) %*% Sigma %*% solve(D)
true_R <- matrix(0,5,5)
for(i in 1:5){
  for(j in 1:5){
    if(i == j){
      true_R[i,j] <- 1
    }
    else{
      true_R[i,j] <- est_corr(sig1 = Sigma[i,i], sig2 = Sigma[j,j], rho = original_R[i,j])
    }
  }
}
# true_R <- R
est_R <- cor(sample_rbind_para)
true_R
###############################################################
## MCSE 
###############################################################
### theta1
sd(sample_rbind_para[,1])/ess_para1
### theta2
sd(sample_rbind_para[,2])/ess_para2
### theta3
sd(sample_rbind_para[,3])/ess_para3
### theta4
sd(sample_rbind_para[,4])/ess_para4
### theta5
sd(sample_rbind_para[,5])/ess_para5
###############################################################
## PSRF: R^hat
###############################################################
sample_cbind_para <- cbind(sample_array_elim[,1,],sample_array_elim[,2,],sample_array_elim[,3,],sample_array_elim[,4,])
J = 4
L = N-Elim_n
col_index <- matrix(0,length(x_init),Chain_n)
for(j in 1:length(x_init)){
  col_index[j,] <- which(colnames(sample_cbind_para) == colnames(sample_rbind_para)[j])
}
for(j in 1:length(x_init)){
  x_j_mean <- apply(sample_cbind_para, 2, mean)[col_index[j,]]
  x_mean <- apply(sample_rbind_para, 2, mean)[j]
  B <- L/(J-1) * sum((x_j_mean - x_mean)^2)
  sj2 <- 1/(L-1) * apply((sample_cbind_para[,col_index[j,]] - x_j_mean)^2,2,sum)
  W <- 1/J * sum(sj2)
  R <- ((L-1)/L * W + 1/L * B + B/(L*J))/W
  print(R)
} # 0.9999423 1.000091
###############################################################
## REM 是在所有维度上近似变量期望的误差总结
###############################################################
REM <- sum(abs(true_mean - est_mean))/sum(abs(true_mean)) # 0.003681218

