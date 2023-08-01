
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

double k_cpp(double x1, double x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
  // return 1;
}

// [[Rcpp::export]]
arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
  arma::mat res(x1.size(), x2.size());
  
  for(int i = 0; (unsigned)i < x1.size(); i++){
    for(int j = 0; (unsigned)j < x2.size(); j++){
      res(i,j) = k_cpp(x1[i],x2[j], a, l);
    }  
  }
  
  return res;
}

double k2_cpp(arma::rowvec x1, arma::rowvec x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-( pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2) ) /(2*pow(l,2)));
}

// [[Rcpp::export]]
arma::mat K2(arma::mat x1, arma::mat x2, double a, double l){
  arma::mat res(x1.n_rows, x2.n_rows);
  
  for(int i = 0; (unsigned)i < x1.n_rows; i++){
    for(int j = 0; (unsigned)j < x2.n_rows; j++){
      res(i,j) = k2_cpp(x1.row(i),x2.row(j), a, l);
    }  
  }
  
  return res;
}

//// GHK FUNCTIONS 

double sim_rnormtrunc(double a, double b, double mu, double sigma){
  
  double alpha = (a - mu) / sigma;
  double beta = (b - mu) / sigma;
  double u = R::runif(0, 1);
  return (sigma * R::qnorm(u * (R::pnorm(beta, 0, 1, 1, 0) - R::pnorm(alpha, 0, 1, 1, 0)) + 
                             R::pnorm(alpha, 0, 1, 1, 0), 0, 1, 1, 0) + mu);
  
}

// [[Rcpp::export]]
arma::vec simulate_eta_cpp(arma::vec y_sign, arma::mat C, arma::vec Xbeta){
  
  int n = y_sign.size();
  arma::vec eta = arma::zeros(n);
  arma::vec bounds(2);
  for (int i = 0; i < n; i++) {
    
    if(y_sign[i] == -1){
      bounds[0] = -100; bounds[1] = 0;
    } else {
      bounds[0] = 0; bounds[1] = 100;
      
    }
    
    for(int l = 0; l <= i; l++){
      bounds[0] -= Xbeta[l];
      bounds[1] -= Xbeta[l];
    }
    for(int l = 0; l < i; l++){
      bounds[0] -= C(i, l) * eta[l];
      bounds[1] -= C(i, l) * eta[l];
      // bounds <- (bounds - (sum(C[i,idxes] * eta[idxes]) + Xbeta[i]) ) / C[i,i]
    }
    bounds = bounds / C(i, i);
    
    eta[i] = sim_rnormtrunc(bounds[0], bounds[1], 0, 1);
    
  }
  
  return eta;
  
}

// [[Rcpp::export]]
double compute_l(arma::vec y_sign, arma::vec Xbeta, arma::mat C){
  
  double prod_l = 1;
  int d = y_sign.size();
  arma::vec eta(d);
  arma::vec bounds(2);
  for(int j = 0; j < d; j++){
    
    if(y_sign[j] == -1){
      bounds[0] = -100; bounds[1] = 0;
    } else {
      bounds[0] = 0; bounds[1] = 100;
    }
    
    bounds[0] -= Xbeta[j];
    bounds[1] -= Xbeta[j];
    for(int l = 0; l < j; l++){
      bounds[0] -= C(j, l) * eta[l];
      bounds[1] -= C(j, l) * eta[l];
    }
    bounds = bounds / C(j, j);
    
    eta[j] = sim_rnormtrunc(bounds[0], bounds[1], 0, 1);
    
    prod_l *= (R::pnorm(bounds[1], 0, 1, 1, 0) - 
                 R::pnorm(bounds[0], 0, 1, 1, 0));
    
  }
  
  return prod_l;
  
}

// [[Rcpp::export]]
arma::vec computeloglik(int B, arma::mat y_sign, arma::mat Xbeta, arma::mat C){
  
  int n = y_sign.n_rows;
  arma::vec loglik = arma::zeros(n);
  for(int i = 0; i < n; i++){
    arma::vec y_sign_i = arma::conv_to<arma::vec>::from(y_sign.row(i));
    arma::vec Xbeta_i = arma::conv_to<arma::vec>::from(Xbeta.row(i));
    for(int l = 0; l < B; l++){
      loglik[i] += (1.0 / B) * compute_l(y_sign_i, Xbeta_i, C);
    }  
  }
  
  return loglik;
}
