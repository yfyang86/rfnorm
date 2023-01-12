#' Multivriate Folded Normal Random Number Generation
#' 
#' Generate multivriate Folded Normal random numbers
#' @param n sample size
#' @param mu mean vector
#' @param Sigma Variance matrix, should be postive, set checkSigma = TRUE to check the postiveness
#' @param method "Gibbs", "HMC", "mvtnorm", default = "Gibbs"
#' @param init initial value (matrix) for Gibbs sampler with p x chain_n dimesions. p is the variable dimension, chain_n is the number of chains.
#' @param gibbs_length iteration length for Gibbs sampler, default = 1000
#' @param gibbs_chain number of chains for Gibbs sampler, default = 4
#' @param checkSigma TRUE/FALSE If TRUE, then mvrfnorm checks whether the variance matrix is positive or not.
#' @return sample array list
#'          - method "Gibbs" or "HMC" or "mvtnorm"
#' @examples
#' ux = matrix(rnorm(4), 2, 2)
#' Sig = ux %*% t(ux)
#' sample_array = mvrfnorm(100,
#' c(0,0),
#' Sig,
#' init=matrix(c(
#'     0.5, 1.5,
#'     0.5, 1.5,
#'     1.5, 0.5,
#'     1.5,1.5
#'     ), ncol = 2),
#' gibbs_chain =  4)
mvrfnorm <- function(n, mu, Sigma, method="Gibbs", 
                    init = matrix(rep(0.5, length(mu)*4), nrow = 4), 
                    gibbs_length=1000, gibbs_chain = 4, 
                    checkSigma = FALSE){
if (sum(method %in% c("Gibbs", "HMC", "mvtnorm")) != 1){
    method <- "Gibbs"
}

p = length(mu)
init = matrix(init, ncol = p)


## p check
if(p < 2){
    stop("Dimension should be larger than one!!!");
}

### n chain
if(dim(init)[1] != gibbs_chain){
    stop('Number of chain does not match!!!');
}

if (length(mu) != dim(init)[2]){
    stop("Length of initial value must be the same as the dimension!!!")
}

## assert Variance-Cov Matrix
if(checkSigma){
    if( sum(svd(Sigma)$d > 0 ) < p){
    stop("Variance Covariance Matrix is not positive definite!!!")
}}

if(method == "mvtnorm") return(
    list(
        method = "mvtnorm",
        value = abs(
            rmvnorm(n, mu, Sigma)
        )
    )
)

if(method == "HMC") {
    stop("Release HMC in the furture.");
} 

#### Gibbs DIM == 2####
if(p != 2){
    stop("Release dim > 2 with Gibbs sampler in the furture.");
}
sig1 <- sqrt(Sigma[1, 1])
sig2 <- sqrt(Sigma[2, 2])
rho <- Sigma[1, 2]/sig1/sig2


j = 1
x_init = matrix(init[j,], ncol = p)
result = rcpp_Gibbs_2d(gibbs_length, x_init, sig1, sig2, rho)
sample_p1 <- result[,1]
sample_p2 <- result[,2]

for (j in 2:gibbs_chain){
    result = rcpp_Gibbs_2d(n, x_init, sig1, sig2, rho)
    sample_p1 <- c(sample_p1, result[,1])
    sample_p2 <- c(sample_p2, result[,2])
    if(j != gibbs_chain) x_init = matrix(init[j,], ncol = p)
}

dim2 <- paste0("chain:", 1:gibbs_chain)
dim3 <- paste0("theta", 1:p)
sample_array <- array(c(sample_p1,sample_p2),c(n, gibbs_chain, p),dimnames = list(NULL,dim2,dim3))
names(dimnames(sample_array)) <- c("Iteration", "Chain", "Parameter")

return(list(method = method, value = sample_array));
}
