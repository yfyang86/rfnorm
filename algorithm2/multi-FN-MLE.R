##multi-dimensional FN distribution mle 
library(MASS)
library(mvtnorm)

permsign <- function(n) {
    if (n < 2) stop("n>=2")
    if (n == 2) {
        return(
            matrix(
                c(1, 1, 1, -1, -1, 1, -1, -1),
                byrow = TRUE, ncol = 2));
    }
    ex = permsign(n - 1)
    return(rbind(cbind(1, ex), cbind(-1, ex)));
}

prodbyrow <- function(mat, vec) {
        t(t(mat) * vec)
}

permute_f <- function(xvec, mu, sigma2){
    n <- length(mu)
    if (n < 2) stop("dim >= 2")
    sign_vec <- permsign(n)
    m <- nrow(sign_vec)
    pdf_val <- 0.
    i_row <- 1
    while (i_row <= m) {
        pdf_val <- pdf_val + dmvnorm(prodbyrow(xvec, sign_vec[i_row, ]), 
                              mean = mu, sigma = sigma2)
        i_row = i_row + 1
        }
    return(pdf_val);
}

# mu, covmat for 2 dim only
est_helper <- function(params) {
    u1 <- params[1]
    u2 <- params[2]
    s1 <- params[3]
    s2 <- params[4]
    rho <- params[5]
    list(mu = c(u1, u2),
         sigma2 = matrix(
        c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2),
        2, 2))
}
# ll for 2 dim only
loglikdim2 <- function(params) {
    u1 <- params[1]
    u2 <- params[2]
    s1 <- exp(params[3])
    s2 <- exp(params[4])
    rho <- params[5]
    ttmmpp = est_helper(c(u1, u2, s1, s2, rho))
    return(-sum(log(permute_f(matrix(c(x, y), ncol = 2), ttmmpp$mu, ttmmpp$sigma2))));
}


# mu, covmat by row (upper triangle),  dim >= 2
est_helper.G2 <- function(params) {
    pp = length(params)

    df_f <- Vectorize(function(p){p^2/2 + 3*p/2 - pp})
    p = which(df_f(1:100) == 0)
    mu = params[1:p]
    sigma2 = matrix(0, p, p)
    s = p + 1
    for(i in 1:p){
        for(j in i:p){
            if (i == j){
                sigma2[i, j] = params[s]/2
            }else{
                sigma2[i, j] = params[s]
            }
            s = s + 1
        }
    }
    return(list(p = p, mu = mu, sigma2 = sigma2 + t(sigma2)));
}

# log-likelihood; mu, covmat by row (upper triangle),  dim >= 2
loglik.G2 <- function(params, ...) {
    ttmmpp = est_helper.G2(params)
    return(-sum(log(permute_f(matrix(inputx, ncol = ttmmpp$p), ttmmpp$mu, ttmmpp$sigma2))));
}


##????figure final ?��?2ά??3άregion?۵???̬??mle
data1=read.csv(file="data1.csv",header=T,sep=",")
data1=data1[,-1]
data2=read.csv(file="data2.csv",header=T,sep=",")
data2=data2[,-1]
##2άregion
#########2 dim folded normal nlm
f1=function(x,y,u1,u2,s1,s2,rho) 1/(2*pi*s1*s2*sqrt(1-rho^2))*(exp(-0.5/(1-rho^2)*((x-u1)^2/s1^2)-2*rho*(x-u1)*(y-u2)/(s1*s2)+(y-u2)^2/s2^2)+exp(-0.5/(1-rho^2)*((x+u1)^2/s1^2)-2*rho*(x+u1)*(y+u2)/(s1*s2)+(y+u2)^2/s2^2)+exp(-0.5/(1-rho^2)*((x+u1)^2/s1^2)+2*rho*(x+u1)*(y-u2)/(s1*s2)+(y-u2)^2/s2^2)+exp(-0.5/(1-rho^2)*((x-u1)^2/s1^2)+2*rho*(x-u1)*(y+u2)/(s1*s2)+(y+u2)^2/s2^2))
minusl3=function(params,r){
  u1=params[1]
  u2=params[2]
  s1=params[3]
  s2=params[4]
  rho=params[5]
  -sum(log(f1(x,y,u1,u2,s1,s2,rho)))
}

a1=matrix(0,7,5)
a1=as.data.frame(a1)
for(i in 1:7){
  start=c(mean(as.numeric(data1[2*i-1,])),mean(as.numeric(data1[2*i,])),sd(as.numeric(data1[2*i-1,])),sd(as.numeric(data1[2*i,])),cor(as.numeric(data1[2*i-1,]),as.numeric(data1[2*i,])))
  x=data1[2*i-1,];y=data1[2*i,]
  nlm(minusl3,start)
  a1[i,]=nlm(minusl3,start)$estimate
}

a1
write.csv(a1,"a1_par.csv")

##3άregion
a2_par=matrix(0,6,9)
a2_par=as.data.frame(a2_par)
for(i in 1:6){
  c=data2[(3*i-2):(3*i),]
  mu=c(mean(as.numeric(c[1,])),mean(as.numeric(c[2,])),mean(as.numeric(c[3,])))
  inputx = t(c)
  sss=cov(t(c))
  init=c(mu,c(sss[1:3], sss[5:6], sss[9]))
  a2_par[i,]=nlm(loglik.G2, init)$estimate
}
a2_par
write.csv(a2_par,file="a2_par.csv")






