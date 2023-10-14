
######################################################################################################
############################################FN MLE simulation#########################################
######################################################################################################
## load libraries 
source("utils.R")
library(MASS)
library(mvtnorm)

#################################an example of estimating parameters##################################
## SImulate Data
set.seed(123)
n <- 100
u1 <- 4; u2 <- 6; sig1 <- 1; sig2 <- 2; rho <- 0.2
mean <- c(u1, u2)
sigma2 <- matrix(c(sig1^2, rho * sig1 * sig2, rho * sig1 * sig2, sig2^2),
                 2, 2)
dat <- mvrnorm(n, mean, sigma2)
## Initilization
covmat <- cov(dat)
init <- c(mean(dat[,1]), mean(dat[,2]), chol(covmat)[c(1, 3, 4)])
inputx = abs(dat)

result <- optim(init, loglik.G2cholesky, method = "BFGS")

result$par
## Covariance Matrix
est_helper.G2(result$par, method = "chole")[["sigma2"]]
#        [,1]      [,2]
#[1,] 0.9411666 0.3998349
#[2,] 0.3998349 3.2814058

######################################################################################################
###################################2d FN MLE simulation###############################################
######################################################################################################

#####n=20,30,40,50,60,70,80,90,100; mu=(2.5,2.5),(5,5),(7.5,7.5),(10,10),(12.5,12.5).
#set.seed(123456)
result_df <- data.frame()
for(i in 1:1000){
  n <- 20
  mu=c(2.5,2.5)
  ss <- matrix(c(25,5,5,25),2,2)
  dat <- mvrnorm(n, mu, ss)

  covmat <- cov(dat)
  init <- c(mean(dat[,1]), mean(dat[,2]), covmat[c(1, 3, 4)])
  inputx = abs(dat)
  
  optim(init, loglik.G2, method = "BFGS")$par
  
  
  fit=optim(init, loglik.G2, method = "BFGS",hessian=T)
  
  fisher_info<-solve(fit$hessian)
  prop_sigma<-sqrt(diag(fisher_info))
  prop_sigma<-diag(prop_sigma)
  a=diag(prop_sigma)
  prop_sigma=a
  upper<-fit$par+1.96*prop_sigma
  lower<-fit$par-1.96*prop_sigma
  
  true_para=c(mu,ss[c(1,3,4)])
  
  p=c()
  for(k in 1:5){
    p[k]=(lower[k]<=true_para[k] & true_para[k]<=upper[k])
  }
  result_df=rbind(result_df, p)
}

##calculate coverage rate of parameters
p=c()
for(k in 1:5){
  p[k]=sum(result_df[which(result_df[,k]==TRUE),k])/1000
}
round(p,2)

######################################################################################################
#########################################4d FN MLE simulation#########################################
######################################################################################################

#####n=20,30,40,50,60,70,80,90,100; mu=(2.5,2.5,2.5,2.5),(5,5,5,5),(7.5,7.5,7.5,7.5),(10,10,10,10),(12.5,12.5,12.5,12.5).
#set.seed(123456)
result_df <- data.frame()
for(i in 1:1000){
  n <- 20
  mu=c(2.5,2.5,2.5,2.5)
  ss=matrix(rep(5,16),4,4)
  diag(ss)=25
  
  dat <- mvrnorm(n, mu, ss)
  sss <- cov(dat)
  init = c(mean(dat[,1]), mean(dat[,2]), mean(dat[,3]), mean(dat[,4]),
           sss[c(1, 5, 9, 13, 6, 10, 14, 11, 15, 16)])
  inputx = abs(dat)
  
  optim(init, loglik.G2, method = "BFGS")$par
  
  fit=optim(init, loglik.G2, method = "BFGS",hessian=T)
  
  fisher_info<-solve(fit$hessian)
  prop_sigma<-sqrt(diag(fisher_info))
  prop_sigma<-diag(prop_sigma)
  a=diag(prop_sigma)
  prop_sigma=a
  upper<-fit$par+1.96*prop_sigma
  lower<-fit$par-1.96*prop_sigma
  
  true_para=c(mu,ss[c(1, 5, 9, 13, 6, 10, 14, 11, 15, 16)])
  p=c()
  for(k in 1:14){
    p[k]=(lower[k]<=true_para[k] & true_para[k]<=upper[k])
  }
  result_df=rbind(result_df, p)
}
p=c()
for(k in 1:14){
  p[k]=sum(result_df[which(result_df[,k]==TRUE),k])/1000
}
round(p,2)

####################################################################################################
#######################################real data application########################################
####################################################################################################
source("utils.R")
library(MASS)
library(ggplot2)
library("shape")
library("MASS") 

data <- read.csv("bmi.csv",header=T,sep=" ") ##R package VGAM bmi data
data=data[,c(2,3)]
colnames(data)=NULL
colnames(data)=c("age","bmi")



# kernel density estimate
kde_estimation <- kde2d(data[, 1], data[, 2])
length(kde_estimation$x)
dim(kde_estimation$z)



#estimate mle
dat=data.frame(data$age,data$bmi)

dat=as.matrix(data)

inputx = dat

sss = chol(cov(dat))

init = c(mean(dat[,1]), mean(dat[,2]),
         sss[c(1,3,4)])
fit <- optim(init, loglik.G2cholesky,hessian=T)

muest=fit$par[c(1,2)];sigma=matrix(c(fit$par[3],0,fit$par[4],fit$par[5]),2,2)
fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
prop_sigma<-diag(prop_sigma)
a=diag(prop_sigma)
prop_sigma=a


sigmaest=est_helper.G2(fit$par, method = "chole")[["sigma2"]]

x <- seq(min(data$age), max(data$age), length.out = 25)
y <- seq(min(data$bmi), max(data$bmi), length.out = 25)

mu1=muest[1]

mu2=muest[2]

summary(data)

sgm1=sqrt(sigmaest[1,1])

sgm2=sqrt(sigmaest[2,2])

rou=sigmaest[1,2]/(sgm1*sgm2)



f1=function(x,y){(1.0/(2.0*pi*sgm1*sgm2*sqrt(1-rou^2)))*(exp((-1.0/(2.0*(1-rou^2)))*((((x-mu1)^2)/(sgm1^2))-(2*rou*(x-mu1)*(y-mu2)/(sgm1*sgm2))+(((y-mu2)^2)/sgm2^2)))+exp((-1.0/(2.0*(1-rou^2)))*((((x+mu1)^2)/(sgm1^2))-(2*rou*(x+mu1)*(y-mu2)/(sgm1*sgm2))+(((y-mu2)^2)/sgm2^2)))
                                                         +exp((-1.0/(2.0*(1-rou^2)))*((((x-mu1)^2)/(sgm1^2))-(2*rou*(x-mu1)*(y+mu2)/(sgm1*sgm2))+(((y+mu2)^2)/sgm2^2)))+ +exp((-1.0/(2.0*(1-rou^2)))*((((x+mu1)^2)/(sgm1^2))-(2*rou*(x+mu1)*(y+mu2)/(sgm1*sgm2))+(((y+mu2)^2)/sgm2^2))))
}

z <- outer(x, y, f1) #generate pdf of folded normal


png("combined_plot.png", width =18 , height = 8, units = "in", res = 300)

layout(matrix(c(1, 2), nrow = 1))
par(mgp = c(3, 2, 5))  
persp(x, y, z, theta = 60, phi = 10, expand = 0.6, r = 180, ltheta = 0, shade = 0.5, ticktype = "detailed", xlab = "Age", ylab = "Body Mass Index", zlab = "Density", col = "lightblue", main = "", cex.axis = 1.25,cex.lab=1.75)


par(mgp = c(3, 2, 5))  
persp(kde_estimation$x, kde_estimation$y, kde_estimation$z, theta = 60, phi = 10, expand = 0.6, r = 180, ltheta = 0, shade = 0.5, ticktype = "detailed", xlab = "Age", ylab = "Body Mass Index", zlab = "Density", col = "lightgray", main = "", cex.axis = 1.25,cex.lab=1.75)
dev.off()
