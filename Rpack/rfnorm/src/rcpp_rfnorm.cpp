#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

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
    return result;
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



//' Generate the Symbolic Matrix for Folded Normal
//'
//' @param n dimensions, 1 <= n <= 16
//' @return matrix of 2^n*n
// [[Rcpp::export]]
 arma::mat rcpp_symbolic_matrix(int n){
    if (n<0) n=1;
    if (n>16) n=16;
    int s = 2 << (n-1);
    int k = 0;
    arma::mat re(s, n);
    if (n == 1) {
        re(0,0) = 1; re(1, 0) = -1;
        return re;
    }
    if ( n == 2){
              re(0,0) = 1;  re(0,1) = 1;
              re(1,0) = 1;  re(1,1) = -1;
              re(2,0) = -1; re(2,1) = 1;
              re(3,0) = -1; re(3,1) = -1;
        return re;
    }else{
        arma::mat re2 = rcpp_symbolic_matrix(n-1);
        for(int i = 0; i < s; ++i){
            for (int j = 0; j < (n-1); j++){
                k = floor(i/2);
                re(i, j) = re2(k,j);
            }
            re(i, n-1) = 2 * (i%2) - 1;
        }
    }
    return re;
}

void print_arma_mat(arma::mat x){
    std::cout<<" "<<std::endl;
    for(int i =0; i < x.n_rows; ++i){
         for(int j =0; j < x.n_cols; ++j){
             std::cout<< std::setw(6) <<x(i, j)<<"\t";
         }
         std::cout<<" "<<std::endl;
    }
     std::cout<<" "<<std::endl;
}

void print_arma_vec(arma::vec x){
    std::cout<<" "<<std::endl;
    for(int i =0; i < x.n_elem; ++i){
             std::cout<< std::setw(6) <<x(i)<<"\t";
    }
    std::cout<<" "<<std::endl;
}



// [[Rcpp::export]]
double rcpp_mixing_sample(int N, int d, int pos, NumericVector x, arma::mat Sigma){
    pos = pos - 1;
    int com = (2 << (d-1));
    double det = arma::det(Sigma);
    arma::mat symbolic = rcpp_symbolic_matrix(d);
    arma::mat Sigma_ac = arma::inv(Sigma) * det;

    arma::uvec row_sub_ind(d-1);
    double cache_ii = 0;
    int ss = 0;
    for (int i = 0; i < d; ++i){
        if (i != pos) {
            row_sub_ind(ss) = i;
            ss++;
        }
    }

    arma::uvec row_sub_ind_k = {0};
    double det_e = arma::det(Sigma.submat(row_sub_ind, row_sub_ind));
    arma::mat Sigma_e_ac = det_e * arma::inv(Sigma.submat(row_sub_ind, row_sub_ind));

    arma::vec eps(com, arma::fill::zeros);
    for(int k = 0; k < com; ++k){
        for(int i =0; i < (d-1); ++i){
            for(int j =0; j < (d-1); ++j){
                row_sub_ind_k(0) = k;
                arma::mat vvv = symbolic.submat(row_sub_ind_k, row_sub_ind);
                eps(k) = eps(k) - 0.5 / det_e * vvv(0, i) * x(i) * vvv(0, j) * x(j) * Sigma_e_ac(i,j);
            }
        }
    }

    arma::mat eps_matrix(com, com);
    for(int i = 0; i < com ; ++i){
         for(int j = 0; j < com ; ++j){
             eps_matrix(i, j) = eps(j);
    }
    }
    for(int i = 0; i < com ; ++i){
        cache_ii = eps_matrix(i,i);
        for(int j = 0; j < com ; ++j){
            eps_matrix(i, j) = eps_matrix(i,j) - cache_ii;
        }
    }


    eps_matrix = arma::exp(eps_matrix);

    arma::vec normal_con(com, arma::fill::zeros);
    for(int i = 0; i < com ; ++i){
        for(int j = 0; j < com ; ++j){
            normal_con(i) += eps_matrix(i,j);
        }
    if (normal_con(i) < 1e-16){
        normal_con(i) = 0.;
    }else{
        normal_con(i) = 1/normal_con(i);
    }
    }

    double rng = runif(1, 0, 1)(0);
    int choice = 1;
    if( rng >=normal_con(0)){
    for(; choice < com; choice++){
       normal_con(choice) +=  normal_con(choice-1);  
        if(rng < normal_con(choice)) break;
    }
    }else{
    choice = 0;
    }

    double sig = det / det_e ;
    arma::vec Sigma_ac_Vec(d-1, arma::fill::zeros);
    ss = 0 ;
    for (int i=0; i<d; i++){
        if(i != pos){
            Sigma_ac_Vec(ss) = (Sigma_ac(pos, i));
            ss++;
    }}

    for (int i=0; i<(d-1); i++){
        Sigma_ac_Vec(i) = x(i) * Sigma_ac_Vec(i);   
    }

    double u = 0;
    int ii = 0;
    for (int i=0; i<d; i++){
        if(i != pos){
            u += symbolic(choice, i) * Sigma_ac_Vec(ii);
            ii++;
        }
    }
    u = u/det_e;

    return(fabs( rnorm(1, u, sqrt(sig))(0) ));
}




//' Generate Multivariate (p>2) Gibbs Sampler
//'
//' @param N Gibbs runs
//' @param x_init initial value of length p, p > 2
//' @param Sig 
//' @return Gibbs samples matrix of N * p
// [[Rcpp::export]]
arma::mat rcpp_Gibbs_nd(int N, NumericVector x_init, arma::mat Sig){
    int p = x_init.size();
    arma::mat x_seq(N, p, arma::fill::zeros);
    std::vector<NumericVector> x_seq_ar;
    x_seq_ar.resize(p-1);
    int jj;
    for (int i = 0; i < (p-1); i++){
         for (int j = 0; j < (p-1); j++) x_seq_ar[i].push_back(0);
    }
    for (int i = 0; i < p; i++){
        x_seq(0,i) = x_init(i);
    }
    for (int i = 1; i < N; i++){
        for (int ii = 0; ii < (p-1); ii++){
            jj = 0;
            /*
            Gibbs: init 0
            */
            if (ii  == 0){
                for (int j = 0; j < p; j++) {
                if(j != ii) {
                    x_seq_ar[ii][jj] = x_seq(i-1, j+1);
                    jj ++;
                }}
            }
            /*
            Gibbs: 0 <  ... < (p-1)
            */
            if (ii < (p-1) ){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
                for (int j = 0; j < ii; j++) {
                    x_seq_ar[ii][jj] = x_seq(i, j);
                    jj ++;
                }
                for (int j = ii; j < p; j++){
                     x_seq_ar[ii][jj] = x_seq_ar[0](i, j);
                    jj ++;
                }
            }
            /*
            Gibbs: Last p-1
            */
            if (ii == (p-1)){
                 for (int j = 0; j < (p-1); j++){
                     x_seq_ar[p-1][j] = x_seq(i, j);
                 }
            }
            // Gibbs Sampling:
            x_seq(i, ii) = rcpp_mixing_sample(N, p, ii + 1, x_seq_ar[ii], Sig);
        }
    }
return x_seq;
}