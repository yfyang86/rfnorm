###############################################################
# HMC and Gibbs Implementation on 100-dimensional Independent Half Normal Distribution
# original author: Jin
###############################################################
library(MASS)
library(ggplot2)
library(patchwork)
###############################################################
## indicator function of random vector X
###############################################################
Indicator <- function(x){
  return(ifelse(FALSE %in% ifelse(x >= 0, TRUE, FALSE), 0, 1))
}
################################################################################################################
## Hamiltonian Monte Carlo                                                        111111
################################################################################################################
## Hamiltonian Monte Carlo Function
###############################################################
HMC = function (U, grad_U,L, current_q)
{
  q = current_q
  p = matrix(rnorm(length(q),0,1), ncol = 1)
  current_p = p
  # epsilon <- runif(1, 0.0104,0.0156)
  epsilon <- 0.1
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept = 0
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) # ensure position vairable big than zero
  {
    accept = 1
    return (list(abs(q), accept))  # accept
  }
  else
  {
    accept = 0
    return (list(current_q, accept))  # reject
  }
}
###############################################################
# target 100-dimensional Independent foleded normal distribution
###############################################################
## potential energy
###############################################################
U_P <- function(q)
{ 
  value <- 0.5 * (matrix(q, nrow = 1) %*% inverse %*% matrix(q, ncol = 1))
  return((value - log((1/2 * pi) ^ (-length(q) / 2) * det(sigma) ^ (-1 / 2))))
}
###############################################################
## gradient 
###############################################################
dU <- function(q)
{
  return(inverse %*% matrix(q, ncol = 1)) 
}
###############################################################
## input 
###############################################################
set.seed(1010)
q_init <- runif(100, 0, 2)
sigma <- diag(seq(0.01, 1, 0.01))
inverse <- solve(sigma)
# epsilon uniformly from interval(0.0104,0.0156) L=50效果几乎相同
L = 100
N = 10000
###############################################################
# simulation bivariate half normal distribution by HMC
###############################################################
## simulation
###############################################################
q_HMC <- matrix(NA, nrow = length(q_init), ncol = N) # position matrix
# sample_p1 <- sample_p2 <- 0
# Chain_all <- list(Chain1 = 0, Chain2 = 0, Chain3 = 0, Chain4 = 0)
isAccept <- 0
t1 <- Sys.time()
for (i in 1:N) {
  result = HMC(U = U_P, grad_U = dU, L = L, current_q = q_init)
  q_HMC[,i] = result[[1]]
  isAccept = isAccept + result[[2]]
  q_init = q_HMC[,i]
}
t2 <- Sys.time()
diff_t_hmc <- round(as.numeric(t2-t1),digits = 4)
###############################################################
## compute acceptance rate
###############################################################
isAccept/N
###############################################################
## estiamte mean & var
###############################################################
true_mean_hmc <- sqrt(seq(0.01, 1.00, 0.01)) * sqrt(2/pi)
true_var_hmc <- seq(0.01, 1.00, 0.01) * (1 - 2/pi)
est_mean_hmc <- apply(q_HMC,1,mean)
est_var_hmc <- apply(q_HMC,1,var)
scatter_hmc <- data.frame(cbind(true_mean_hmc, true_var_hmc, est_mean_hmc, est_var_hmc))
###############################################################
## visualization by scatter
###############################################################
scatter_mean_HMC_1 <-
  ggplot(data = scatter_hmc, aes(x = true_mean_hmc, y = est_mean_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8),
        plot.title = element_text(hjust = 0.5, vjust = 2.0, face = "bold")) +
  expand_limits(y = c(0.0, 1.0)) +
  expand_limits(x = c(0.0, 1.0)) +
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  geom_abline(intercept = 0,slope = 1) +
  xlab("分量的真实均值") +
  ylab("分量的样本均值") +
  labs(title = paste0(diff_t_hmc,"秒, ", 100 * isAccept/N,"%"))
scatter_var_HMC_1 <- 
  ggplot(data = scatter_hmc, aes(x = true_var_hmc, y = est_var_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8)) +
  expand_limits(y = c(0.0, 0.6)) +
  expand_limits(x = c(0.0, 0.6)) +
  scale_x_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  scale_y_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("分量的真实方差") +
  ylab("分量的样本方差") 

################################################################################################################
## Hamiltonian Monte Carlo                                                          222222222
################################################################################################################
## Hamiltonian Monte Carlo Function
###############################################################
HMC = function (U, grad_U,L, current_q)
{
  q = current_q
  p = matrix(rnorm(length(q),0,1), ncol = 1)
  current_p = p
  # epsilon <- runif(1, 0.01, 0.02)
  # epsilon <- runif(1, 0.0104,0.0156)
  epsilon <- 0.025
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept = 0
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) # ensure position vairable big than zero
  {
    accept = 1
    return (list(abs(q), accept))  # accept
  }
  else
  {
    accept = 0
    return (list(current_q, accept))  # reject
  }
}
###############################################################
# target 100-dimensional Independent foleded normal distribution
###############################################################
## potential energy
###############################################################
U_P <- function(q)
{ 
  value <- 0.5 * (matrix(q, nrow = 1) %*% inverse %*% matrix(q, ncol = 1))
  return((value - log((1/2 * pi) ^ (-length(q) / 2) * det(sigma) ^ (-1 / 2))))
}
###############################################################
## gradient 
###############################################################
dU <- function(q)
{
  return(inverse %*% matrix(q, ncol = 1)) 
}
###############################################################
## input 
###############################################################
set.seed(1010)
q_init <- runif(100, 0, 2)
sigma <- diag(seq(0.01, 1, 0.01))
inverse <- solve(sigma)
# epsilon uniformly from interval(0.0104,0.0156) L=50效果几乎相同
L = 100
N = 10000
###############################################################
# simulation bivariate half normal distribution by HMC
###############################################################
## simulation
###############################################################
q_HMC <- matrix(NA, nrow = length(q_init), ncol = N) # position matrix
# sample_p1 <- sample_p2 <- 0
# Chain_all <- list(Chain1 = 0, Chain2 = 0, Chain3 = 0, Chain4 = 0)
isAccept <- 0
t1 <- Sys.time()
for (i in 1:N) {
  result = HMC(U = U_P, grad_U = dU, L = L, current_q = q_init)
  q_HMC[,i] = result[[1]]
  isAccept = isAccept + result[[2]]
  q_init = q_HMC[,i]
}
t2 <- Sys.time()
diff_t_hmc <- round(as.numeric(t2-t1),digits = 4)
###############################################################
## compute acceptance rate
###############################################################
isAccept
###############################################################
## estiamte mean & var
###############################################################
true_mean_hmc <- sqrt(seq(0.01, 1.00, 0.01)) * sqrt(2/pi)
true_var_hmc <- seq(0.01, 1.00, 0.01) * (1 - 2/pi)
est_mean_hmc <- apply(q_HMC,1,mean)
est_var_hmc <- apply(q_HMC,1,var)
scatter_hmc <- data.frame(cbind(true_mean_hmc, true_var_hmc, est_mean_hmc, est_var_hmc))
###############################################################
## visualization by scatter
###############################################################
scatter_mean_HMC_2 <-
  ggplot(data = scatter_hmc, aes(x = true_mean_hmc, y = est_mean_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8),
        plot.title = element_text(hjust = 0.5, vjust = 2.0, face = "bold")) +
  expand_limits(y = c(0.0, 1.0)) +
  expand_limits(x = c(0.0, 1.0)) +
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  geom_abline(intercept = 0,slope = 1) +
  xlab("分量的真实均值") +
  ylab("分量的样本均值") +
  labs(title = paste0(diff_t_hmc,"秒, ", 100 * isAccept/N,"%"))
scatter_var_HMC_2 <- 
  ggplot(data = scatter_hmc, aes(x = true_var_hmc, y = est_var_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8)) +
  expand_limits(y = c(0.0, 0.6)) +
  expand_limits(x = c(0.0, 0.6)) +
  scale_x_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  scale_y_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("分量的真实方差") +
  ylab("分量的样本方差") 

################################################################################################################
## Hamiltonian Monte Carlo                                                          3333
################################################################################################################
## Hamiltonian Monte Carlo Function
###############################################################
HMC = function (U, grad_U,L, current_q)
{
  q = current_q
  p = matrix(rnorm(length(q),0,1), ncol = 1)
  current_p = p
  epsilon <- runif(1, 0.01,0.02)
  # epsilon <- runif(1, 0.0104,0.0156)
  # epsilon <- 0.018
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept = 0
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) # ensure position vairable big than zero
  {
    accept = 1
    return (list(abs(q), accept))  # accept
  }
  else
  {
    accept = 0
    return (list(current_q, accept))  # reject
  }
}
###############################################################
# target 100-dimensional Independent foleded normal distribution
###############################################################
## potential energy
###############################################################
U_P <- function(q)
{ 
  value <- 0.5 * (matrix(q, nrow = 1) %*% inverse %*% matrix(q, ncol = 1))
  return((value - log((1/2 * pi) ^ (-length(q) / 2) * det(sigma) ^ (-1 / 2))))
}
###############################################################
## gradient 
###############################################################
dU <- function(q)
{
  return(inverse %*% matrix(q, ncol = 1)) 
}
###############################################################
## input 
###############################################################
set.seed(1010)
q_init <- runif(100, 0, 2)
sigma <- diag(seq(0.01, 1, 0.01))
inverse <- solve(sigma)
# epsilon uniformly from interval(0.0104,0.0156) L=50效果几乎相同
L = 100
N = 10000
###############################################################
# simulation bivariate half normal distribution by HMC
###############################################################
## simulation
###############################################################
q_HMC <- matrix(NA, nrow = length(q_init), ncol = N) # position matrix
# sample_p1 <- sample_p2 <- 0
# Chain_all <- list(Chain1 = 0, Chain2 = 0, Chain3 = 0, Chain4 = 0)
isAccept <- 0
t1 <- Sys.time()
for (i in 1:N) {
  result = HMC(U = U_P, grad_U = dU, L = L, current_q = q_init)
  q_HMC[,i] = result[[1]]
  isAccept = isAccept + result[[2]]
  q_init = q_HMC[,i]
}
t2 <- Sys.time()
diff_t_hmc <- round(as.numeric(t2-t1),digits = 4)
###############################################################
## compute acceptance rate
###############################################################
isAccept
###############################################################
## estiamte mean & var
###############################################################
true_mean_hmc <- sqrt(seq(0.01, 1.00, 0.01)) * sqrt(2/pi)
true_var_hmc <- seq(0.01, 1.00, 0.01) * (1 - 2/pi)
est_mean_hmc <- apply(q_HMC,1,mean)
est_var_hmc <- apply(q_HMC,1,var)
scatter_hmc <- data.frame(cbind(true_mean_hmc, true_var_hmc, est_mean_hmc, est_var_hmc))
###############################################################
## visualization by scatter
###############################################################
scatter_mean_HMC_3 <-
  ggplot(data = scatter_hmc, aes(x = true_mean_hmc, y = est_mean_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8),
        plot.title = element_text(hjust = 0.5, vjust = 2.0, face = "bold")) +
  expand_limits(y = c(0.0, 1.0)) +
  expand_limits(x = c(0.0, 1.0)) +
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  geom_abline(intercept = 0,slope = 1) +
  xlab("分量的真实均值") +
  ylab("分量的样本均值") +
  labs(title = paste0(diff_t_hmc,"秒, ", 100 * isAccept/N,"%"))
scatter_var_HMC_3 <- 
  ggplot(data = scatter_hmc, aes(x = true_var_hmc, y = est_var_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8)) +
  expand_limits(y = c(0.0, 0.6)) +
  expand_limits(x = c(0.0, 0.6)) +
  scale_x_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  scale_y_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("分量的真实方差") +
  ylab("分量的样本方差") 
################################################################################################################
## Hamiltonian Monte Carlo                                                          4444
################################################################################################################
## Hamiltonian Monte Carlo Function
###############################################################
HMC = function (U, grad_U,L, current_q)
{
  q = current_q
  p = matrix(rnorm(length(q),0,1), ncol = 1)
  current_p = p
  epsilon <- runif(1, 0.01,0.02)
  # epsilon <- runif(1, 0.0104,0.0156)
  # epsilon <- 0.018
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept = 0
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) # ensure position vairable big than zero
  {
    accept = 1
    return (list(abs(q), accept))  # accept
  }
  else
  {
    accept = 0
    return (list(current_q, accept))  # reject
  }
}
###############################################################
# target 100-dimensional Independent foleded normal distribution
###############################################################
## potential energy
###############################################################
U_P <- function(q)
{ 
  value <- 0.5 * (matrix(q, nrow = 1) %*% inverse %*% matrix(q, ncol = 1))
  return((value - log((1/2 * pi) ^ (-length(q) / 2) * det(sigma) ^ (-1 / 2))))
}
###############################################################
## gradient 
###############################################################
dU <- function(q)
{
  return(inverse %*% matrix(q, ncol = 1)) 
}
###############################################################
## input 
###############################################################
set.seed(1010)
q_init <- runif(100, 0, 2)
sigma <- diag(seq(0.01, 1, 0.01))
inverse <- solve(sigma)
# epsilon uniformly from interval(0.0104,0.0156) L=50效果几乎相同
L = 25
N = 10000
###############################################################
# simulation bivariate half normal distribution by HMC
###############################################################
## simulation
###############################################################
q_HMC <- matrix(NA, nrow = length(q_init), ncol = N) # position matrix
# sample_p1 <- sample_p2 <- 0
# Chain_all <- list(Chain1 = 0, Chain2 = 0, Chain3 = 0, Chain4 = 0)
isAccept <- 0
t1 <- Sys.time()
for (i in 1:N) {
  result = HMC(U = U_P, grad_U = dU, L = L, current_q = q_init)
  q_HMC[,i] = result[[1]]
  isAccept = isAccept + result[[2]]
  q_init = q_HMC[,i]
}
t2 <- Sys.time()
diff_t_hmc <- round(as.numeric(t2-t1),digits = 4)
###############################################################
## compute acceptance rate
###############################################################
isAccept
###############################################################
## estiamte mean & var
###############################################################
true_mean_hmc <- sqrt(seq(0.01, 1.00, 0.01)) * sqrt(2/pi)
true_var_hmc <- seq(0.01, 1.00, 0.01) * (1 - 2/pi)
est_mean_hmc <- apply(q_HMC,1,mean)
est_var_hmc <- apply(q_HMC,1,var)
scatter_hmc <- data.frame(cbind(true_mean_hmc, true_var_hmc, est_mean_hmc, est_var_hmc))
###############################################################
## visualization by scatter
###############################################################
scatter_mean_HMC_4 <-
  ggplot(data = scatter_hmc, aes(x = true_mean_hmc, y = est_mean_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8),
        plot.title = element_text(hjust = 0.5, vjust = 2.0, face = "bold")) +
  expand_limits(y = c(0.0, 1.0)) +
  expand_limits(x = c(0.0, 1.0)) +
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  geom_abline(intercept = 0,slope = 1) +
  xlab("分量的真实均值") +
  ylab("分量的样本均值") +
  labs(title = paste0(diff_t_hmc,"秒, ", 100 * isAccept/N,"%"))
scatter_var_HMC_4 <- 
  ggplot(data = scatter_hmc, aes(x = true_var_hmc, y = est_var_hmc)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8)) +
  expand_limits(y = c(0.0, 0.6)) +
  expand_limits(x = c(0.0, 0.6)) +
  scale_x_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  scale_y_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("分量的真实方差") +
  ylab("分量的样本方差") 

# scatter_mean_HMC / scatter_var_HMC
################################################################################################################
## Gibbs sampling
################################################################################################################
###############################################################
## One Dimension Folded Normal Sampling Funcion
###############################################################
FN_Sampling <- function(number, mu, sig){
  sample <- abs(rnorm(number, mean = mu, sd = sqrt(sig)))
  density <- function(x, u, s){
    1/sqrt(2 * pi)/sqrt(s) * (exp(-(x - u)^2/(2 * s)) + exp(-(x + u)^2/(2 * s)) * Indicator(x))
  }
  return(list(sample,density(sample, mu, sig)))
}
###############################################################
# Gibbs Sampling Function
###############################################################
Gibbs <- function(N,length,sigma){
  x_seq <- matrix(0, N, length)
  for (i in 1:N){
    for(j in 1:length){
      x_seq[i,j] <- FN_Sampling(1, 0, sigma[j,j])[[1]]
    }
  }
  return(x_seq)
}
###############################################################
# Gibbs Sampling algorithm (100 dimension)
###############################################################
set.seed(1010)
sigma <- diag(seq(0.01, 1, 0.01))
N <- 10000
length <- 100
t3 <- Sys.time()
result <- Gibbs(N, length, sigma)
t4 <- Sys.time()
diff_t_gibbs <- round(as.numeric(t4-t3),digits = 4)
###############################################################
## estiamte mean & var
###############################################################
true_mean_gibbs <- sqrt(seq(0.01, 1.00, 0.01)) * sqrt(2/pi)
true_var_gibbs <- seq(0.01, 1.00, 0.01) * (1 - 2/pi)
est_mean_gibbs <- apply(result,2,mean)
est_var_gibbs <- apply(result,2,var)
scatter_gibbs <- data.frame(cbind(true_mean_gibbs, true_var_gibbs,est_mean_gibbs,est_var_gibbs))
###############################################################
## visualization by scatter
###############################################################
scatter_mean_Gibbs <-
  ggplot(data = scatter_gibbs, aes(x = true_mean_gibbs, y = est_mean_gibbs)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold")) +
  expand_limits(y = c(0.0, 1.0)) +
  expand_limits(x = c(0.0, 1.0)) +
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), labels = c(paste0("0.0"), seq(0.2, 0.8, 0.2), paste0("1.0"))) +
  geom_abline(intercept = 0,slope = 1) +
  xlab("分量的真实均值") +
  ylab("分量的样本均值")
  # labs(title = paste0("吉布斯采样(",diff_t_gibbs,"秒)"))
scatter_var_Gibbs <- 
  ggplot(data = scatter_gibbs, aes(x = true_var_gibbs, y = est_var_gibbs)) + 
  geom_point(size = 1.5, alpha = 0.55, color = "red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r = 10),size = 8),
        axis.title.x = element_text(vjust=-1,size = 8)) +
  expand_limits(y = c(0.0, 0.6)) +
  expand_limits(x = c(0.0, 0.6)) +
  scale_x_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  scale_y_continuous(breaks = seq(0.0, 0.6, 0.1), labels = c(paste0("0.0"), seq(0.1, 0.6, 0.1))) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("分量的真实方差") +
  ylab("分量的样本方差") 
# scatter_mean_Gibbs / scatter_var_Gibbs
###############################################################
## visualization HMC & Gibbs
###############################################################
A <- scatter_mean_Gibbs | scatter_var_Gibbs
B <- scatter_mean_HMC_1 / scatter_var_HMC_1
C <- scatter_mean_HMC_2 / scatter_var_HMC_2
D <- scatter_mean_HMC_3 / scatter_var_HMC_3
E <- scatter_mean_HMC_4 / scatter_var_HMC_4


library(ggpubr)
gibbs <- ggarrange(scatter_mean_Gibbs, scatter_var_Gibbs, nrow = 1)
annotate_figure(gibbs,
                top = text_grob(paste0("吉布斯采样(",diff_t_gibbs,"秒)"), face = "bold", size = 14,hjust = 0.4)
                )
ggarrange(B,C,D,E, ncol = 2, nrow = 2,labels = "AUTO")

# (scatter_mean_Gibbs | scatter_mean_HMC_1) / (scatter_var_Gibbs | scatter_var_HMC_1) /
#   (scatter_mean_HMC_2 | scatter_mean_HMC_3) / (scatter_var_HMC_2 | scatter_var_HMC_3) + 
#   plot_annotation(tag_levels = 'A')

