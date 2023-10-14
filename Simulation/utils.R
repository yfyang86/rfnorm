# License: Apache 2
# Author: Yifan Yang <yfyang.86@hotmail.com>

library(MASS)
library(mvtnorm)
library(Rcpp)

src <- '
NumericMatrix prodbyrow_c(NumericMatrix mat, NumericVector vec){
    /**
    example:
    prodbyrow_c(cbind(runif(10), -runif(10)), c(-1, 1))
    */
    NumericMatrix result(mat.nrow(), mat.ncol());
    int i = 0;
    int j = 0;
    for(; i < mat.nrow(); ++i){
        for(; j < mat.ncol(); ++j){
            result(i, j) = (vec(j) < 0 ? -mat(i, j) : mat(i, j));
        }
        j = 0;
    }
    return(result);
}
'

src2 <- '
NumericMatrix permsign_c(int n){
    if(n == 2){
        NumericMatrix result(4, 2);
        result(0, 0) = 1; result(0, 1) = 1;
        result(1, 0) = 1; result(1, 1) = -1;
        result(2, 0) = -1; result(2, 1) = 1;
        result(3, 0) = -1; result(3, 1) = -1;
        return (result);
    }
    NumericMatrix tmp = permsign_c(n-1);
    int v = (2<<(n-2));
    int u = (2<<(n-1));
    NumericMatrix result(u, n);
    for(int i = 0; i < v; i++){
        result(i, 0) = 1;
        for(int j = 1; j < n; j++){
            result(i, j) = tmp(i, j-1);
        }
    }
    
    for(int i = v; i < u; i++){
        result(i, 0) = -1;
        for(int j = 1; j < n; j++){
            result(i, j) = tmp(i-v, j-1);
        }
    }

    return(result);
}
'

cppFunction(src)
cppFunction(src2)


#' Generate combination of {1, -1} of length n
permsign <- function(n) {
    if (n < 2) stop("n>=2")
    if (n == 2) {
        return(
            matrix(
                c(1,  1,
                  1, -1,
                 -1,  1,
                 -1, -1),
                byrow = TRUE, ncol = 2));
    }
    ex <- permsign(n - 1)
    return(rbind(cbind(1, ex), cbind(-1, ex)));
}

#' pointwise mult-op by row
prodbyrow <- function(mat, vec) {
        t(t(mat) * vec)
}

#' PDF of multiple F-N 
permute_f <- function(xvec, mu, sigma2){
    n <- length(mu)
    if (n < 2) stop("dim >= 2")
    sign_vec <- permsign_c(n)
    m <- nrow(sign_vec)
    pdf_val <- 0.
    i_row <- 1
    while (i_row <= m) {
        pdf_val <- pdf_val + dmvnorm(prodbyrow_c(xvec, sign_vec[i_row, ]),
                              mean = mu, sigma = sigma2)
        i_row = i_row + 1
        }
    return(pdf_val);
}

#' mu, covmat for 2 dim only
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

#' ll for 2 dim only
loglikdim2 <- function(params) {
    u1 <- params[1]
    u2 <- params[2]
    s1 <- exp(params[3])
    s2 <- exp(params[4])
    rho <- params[5]
    ttmmpp = est_helper(c(u1, u2, s1, s2, rho))
    return(-sum(log(permute_f(matrix(c(x, y), ncol = 2), ttmmpp$mu, ttmmpp$sigma2))));
}


#' mu, covmat by row (upper triangle), 100 >= dim >= 2
#' 
#' Generally, method == "chole" is more stable (Non-negative Matrix),
#' while method == "tri" only takes garantee that the matrix is symtric.
#' 
#' @param params length(params) = p + (p+1)*p/2
#' params[1:p] stores the mean vector:
#' 
#' -    mu_1 mu_2 ... mu_p
#' 
#' params[-(1:p)] stores the triangle elements of the covariance matrix 
#' (or the Cholesky Decomposition upper triangle matrix v s.t. Cov = v' v)
#'  
#'     sig_{1, 1}, sig_{1, 2}, ..., sig_{1, p}
#'      
#'     0         , sig_{2,2} , ..., sig_{2, p}
#'      
#'     0 ...
#'      
#'     0 ...                      , sig_{p, p}
#' 
#' @param method "tri" or "chole"
est_helper.G2 <- function(params, method = 'tri') {
    pp = length(params)

    df_f <- Vectorize(function(p){p^2/2 + 3*p/2 - pp})
    p = which(df_f(1:100) == 0)
    mu = params[1:p]
    sigma2 = matrix(0, p, p)
    s = p + 1
    for(i in 1:p){
        for(j in i:p){
            if(method == 'tri'){
                if (i == j){
                    sigma2[i, j] = params[s]/2
                }else{
                    sigma2[i, j] = params[s]
                }
            }else{
                # Cholesky
                 sigma2[i, j] = params[s]
            }
            s = s + 1
        }
    }
     if(method == 'tri'){
            sigmamat = sigma2 + t(sigma2)
            }else{
            sigmamat = t(sigma2) %*% (sigma2)
            }
    
    return(list(p = p, mu = mu, sigma2 = sigmamat));
}

# log-likelihood; mu, covmat by row (upper triangle),  dim >= 2
loglik.G2 <- function(params) {
    ttmmpp = est_helper.G2(params)
    return(-sum(log(permute_f(matrix(inputx, ncol = ttmmpp$p), ttmmpp$mu, ttmmpp$sigma2))));
}

# log-likelihood; mu, covmat by row (upper triangle),  dim >= 2
loglik.G2cholesky <- function(params) {
    ttmmpp = est_helper.G2(params, 'chol')
    return(-sum(log(permute_f(matrix(inputx, ncol = ttmmpp$p), ttmmpp$mu, ttmmpp$sigma2))));
}

loglik.G2cholesky <- function(params) {
  ttmmpp = est_helper.G2(params, 'chol')
  return(-sum(log(permute_f(matrix(inputx, ncol = ttmmpp$p), ttmmpp$mu, ttmmpp$sigma2))));
}


#' Lower triangle index of a matrix
#' @param dim_p dimension of the matrix
#' @return lower triangle index of a matrix
f_lower_idx <- function(dim_p){
    return(matrix(1:(dim_p^2),
                    byrow = TRUE, dim_p)[lower.tri(
                                        matrix(0, dim_p, dim_p), diag = TRUE)])
}