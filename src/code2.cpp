// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
  
  // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"


using namespace Rcpp;

// // [[Rcpp::export]]
// arma::vec logistic_CPP(arma::vec x){
//   return (1 / (1 + exp(-x)));
// }
// 
// // [[Rcpp::export]] 
// arma::vec psiModel(arma::vec coeff, 
//                    arma::vec X){
//   
//   arma::mat XX = arma::join_rows(arma::ones(X.n_rows), X, pow(X, 2));
//   
//   arma::vec Xbeta = XX * coeff;
// 
//   // arma::vec psi = logistic_CPP(Xbeta);
//   
//   return Xbeta;
//   
// }
// 
// // [[Rcpp::export]] 
// arma::vec pModel(arma::vec coeff, 
//                    arma::vec X){
//   
//   arma::mat XX = arma::join_rows(arma::ones(X.n_rows), X);
//   
//   arma::vec Xbeta = XX * coeff;
// 
//   // arma::vec psi = logistic_CPP(Xbeta);
//   
//   return Xbeta;
//   
// }
// 
// // [[Rcpp::export]]
// double loglik_cpp(arma::vec pars,
//                   arma::vec y, 
//                   arma::vec M, 
//                   arma::mat X_psi, 
//                   arma::mat X_p,
//                   arma::uvec cov_psi,
//                   arma::uvec cov_p,
//                   arma::vec y_occupied,
//                   arma::vec occ,
//                   arma::vec occ_all,
//                   arma::vec sumM){
//   
//   arma::vec coeff_psi = pars.elem(cov_psi);
//   arma::vec coeff_p = pars.elem(cov_p);
//   
//   // arma::vec Xpsibetapsi = X_psi * coeff_psi;
//   // arma::vec Xpbetap = X_p * coeff_p;
// 
//   arma::vec Xpsibetapsi = psiModel(coeff_psi, X_psi);
//   arma::vec Xpbetap = pModel(coeff_p, X_p);
//     
//   arma::vec psi = logistic_CPP(Xpsibetapsi);
//   arma::vec p = logistic_CPP(Xpbetap);
//   
//   arma::uvec occ_all1 = find(occ_all == 1);  
//   arma::uvec occ0 = find(occ == 0);  
//   arma::uvec occ1 = find(occ == 1);  
//   
//   arma::vec p_occupied = p.elem(occ_all1);
//   arma::vec psi_occ0 = psi.elem(occ0);
//   arma::vec psi_occ1 = psi.elem(occ1);
//   
//   // occupancies for occupied
//   double loglik_psi_occ = sum(log(psi_occ1));
//   
//   // detections for occupied
//   double loglik_p_occ = 0;
//   for(int l = 0; l < y_occupied.size(); l++){
//     loglik_p_occ += R::dbinom(y_occupied[l], 1, p_occupied[l], 1);
//   }
//   
//   // arma::vec p_prod = arma::ones(psi_occ0.size());
//   // int sumM2 = 0;
//   // int l2 = 0;
//   // for(int l = 0; l < occ.size(); l++){
//     //   if(occ[l] == 0){
//       //     for(int k = 0; k < M[l]; k++){
//         //       p_prod[l2] *= (1 - p[sumM2 + k]);
//         //     }
//       //     l2 += 1;
//       //   }
//     //   sumM2 += M[l];
//     // }
//   // 
//     // Rcout << p_prod << std::endl;
//   
//   arma::vec p_prod = arma::ones(psi_occ0.size());
//   for(int l = 0; l < occ0.size(); l++){
//     int i = occ0[l];
//     for(int k = 0; k < M[i]; k++){
//       p_prod[l] *= (1 - p[sumM[i] + k]);
//     }
//   }
//   
//   double loglik_ppsi_nocc = 
//     sum(
//       log(
//         (1 - psi_occ0) + psi_occ0 % p_prod
//       )
//     );
//   
//   
//   // double logprior_psi = R::dnorm(beta_psi[0], 0, 2, 1);
//   // double logprior_p = R::dnorm(beta_p[0], 0, 2, 1);
//   
//   double log_grad = - (1.0 / (2 * 4)) * (sum(coeff_psi % coeff_psi) + 
//                                            sum(coeff_p % coeff_p));
//   
//   // logprior_betapsi <- sum(dnorm(beta_psi[-1], 0, 2, log = T))
//   // logprior_betap <- sum(dnorm(beta_p[-1], 0, 2, log = T))
//   
//   // # - (loglik_psi_occ + loglik_p_occ + loglik_ppsi_nocc +
//     // # logprior_psi + logprior_p + logprior_betapsi + logprior_betap)
//     
//     return(- (loglik_psi_occ + loglik_p_occ + loglik_ppsi_nocc + log_grad));
//   
// }
// 
// // [[Rcpp::export]]
// arma::vec gr_loglik_cpp(arma::vec pars,
//                         arma::vec y, 
//                         arma::vec M, 
//                         arma::mat X_psi, 
//                         arma::mat X_p,
//                         arma::uvec cov_psi,
//                         arma::uvec cov_p,
//                         arma::vec y_occupied,
//                         arma::vec occ,
//                         arma::vec occ_all,
//                         arma::vec sumM){
//   
//   int ncov_psi = X_psi.n_cols;
//   int ncov_p = X_p.n_cols;
//   
//   arma::vec coeff_psi = pars.elem(cov_psi);
//   arma::vec coeff_p = pars.elem(cov_p);
//   
//   // arma::vec Xpsibetapsi = X_psi * coeff_psi;
//   // arma::vec Xpbetap = X_p * coeff_p;
//   
//   arma::vec Xpsibetapsi = psiModel(coeff_psi, X_psi);
//   arma::vec Xpbetap = pModel(coeff_p, X_p);
//   
//   arma::vec psi = logistic_CPP(Xpsibetapsi);
//   arma::vec p = logistic_CPP(Xpbetap);
//   
//   arma::uvec occ_all1 = find(occ_all == 1);  
//   arma::uvec occ0 = find(occ == 0);  
//   arma::uvec occ1 = find(occ == 1);
//   
//   arma::vec p_occupied = p.elem(occ_all1);
//   arma::vec psi_occ0 = psi.elem(occ0);
//   arma::vec psi_occ1 = psi.elem(occ1);
//   
//   arma::vec p_prod = arma::ones(psi_occ0.size());
//   for(int l = 0; l < occ0.size(); l++){
//     int i = occ0[l];
//     for(int k = 0; k < M[i]; k++){
//       p_prod[l] *= (1 - p[sumM[i] + k]);
//     }
//   }
//   
//   arma::vec grad_l = arma::zeros(ncov_psi + ncov_p);
//   
//   // occupancies for occupied
//   for(int l = 0; l < occ1.size(); l++){
//     int i = occ1[l];
//     arma::vec grad_current = arma::conv_to<arma::vec>::from(X_psi.row(i)) * 
//       (exp(- Xpsibetapsi[i]) / (1 + exp(- Xpsibetapsi[i])));
//     grad_l.elem(cov_psi) += grad_current;
//   }
//   
//   // detections for occupied
//   for(int l = 0; l < occ_all1.size(); l++){
//     int i = occ_all1[l];
//     arma::vec grad_current = arma::conv_to<arma::vec>::from(X_p.row(i)) *
//       (y[i] * exp(- Xpbetap[i]) - (1 - y[i])) / (1 + exp(- Xpbetap[i]));
//     grad_l.elem(cov_p) += grad_current;
//   }
//   
//   for(int l = 0; l < occ0.size(); l++){
//     
//     int i = occ0[l];
//     
//     double denom = (p_prod[l] + exp(- Xpsibetapsi[i])) *
//       (1 + exp( Xpsibetapsi[i]));
//     
//     arma::vec grad_current = arma::conv_to<arma::vec>::from(X_psi.row(i)) *
//       (p_prod[l] - 1) / denom;
//     
//     grad_l.elem(cov_psi) += grad_current;
//     
//     double term1 = - psi[i] / (psi[i] * p_prod[l] + (1 - psi[i]));
//     
//     for(int k = 0; k < M[i]; k++){
//       
//       int idx = sumM[i] + k;
//       double currentProd = p_prod[l] / (1 - p[idx]);
//       
//       double term2 = exp(- Xpbetap[idx]);
//       double term3 = pow(1 + exp(- Xpbetap[idx]),2);
//       double term4 = term2 / term3;
//       arma::vec grad_current = arma::conv_to<arma::vec>::from(X_p.row(idx)) *
//         currentProd * term1 * term4;
//       
//       grad_l.elem(cov_p) += grad_current;
//       
//     }
//     
//   }
//   
//   for(int l = 0; l < grad_l.size(); l++){
//     grad_l[l] += - (1.0 / (4)) * pars[l];
//   }
//   
//   return (- grad_l);
// }