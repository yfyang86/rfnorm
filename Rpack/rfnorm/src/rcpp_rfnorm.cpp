#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]


// [[Rcpp::export]]
int rcpp_indicator(NumericVector &x){
    for (auto xx:x){
        if (xx<0) return 0;
    }
    return 1;
}


// [[Rcpp::export]]
List rcpp_FN_Sampling_2d(
    int number,
    double mu,
    double sig
){
    NumericVector sample = rnorm(number, mu, sig);
    NumericVector density;
    for (int i = 0; i<sample.size(); i++){
        sample[i] = fabs(sample[i]);
    }
    for(auto x:sample){
        density.push_back(
            density_double(x, mu, sig)
        );
    }
    List result = List::create( sample, density );
}

NumericVector rcpp_FN_Sampling_2d(double mu, double sig){
    NumericVector sample = rnorm(1, mu, sig);
    double sample_val = fabs(sample(0));
    double density =  density_double(sample_val, mu, sig);
    NumericVector re;
    re.push_back(sample_val);
    re.push_back(density);
    return re;
}

double rcpp_rFN_Sampling_2d(double mu, double sig){
    NumericVector sample = rnorm(1, mu, sig);
    double sample_val = fabs(sample(0));
    return sample_val;
}

//' Folded Normal 2D RNG with Gibbs Sampling
//'
//' Generate Random Number for 2-d Folded Normal with Gibss Sampling
//' @param N sample size
//' @param x_init vector of length 2, initial value
//' @param sig1 variance[1,1]
//' @param sig2 variance[2,2]
//' @param rho rho
//' @return matrix of N*2
// [[Rcpp::export]]
NumericMatrix rcpp_Gibbs_2d(int N, NumericVector x_init, double sig1, double sig2, double rho){
    NumericMatrix x_seq(N, 2);
    double v1 = rho * sig1/sig2;
    double v2 = rho * sig2/sig1;
    double v11 = sig1 * (1 - sq(rho));
    double v21 = sig2 * (1 - sq(rho));
    for (int i = 0; i < x_init.size(); ++i){
        x_seq(0, i) = x_init(i);
    }
    for (int i = 1; i < N; ++i){
        x_seq(i, 0) = rcpp_rFN_Sampling_2d((double)x_seq(i-1,1)*v1, v11);
        x_seq(i, 1) = rcpp_rFN_Sampling_2d((double)x_seq(i,0)*v2, v21);
    }
    return x_seq;
}