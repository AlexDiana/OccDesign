// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logistic_CPP
arma::vec logistic_CPP(arma::vec x);
RcppExport SEXP _OccDesign_logistic_CPP(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_CPP(x));
    return rcpp_result_gen;
END_RCPP
}
// loglik_cpp
double loglik_cpp(arma::vec pars, arma::vec y, arma::vec M, arma::mat X_psi, arma::mat X_p, arma::uvec cov_psi, arma::uvec cov_p, arma::vec y_occupied, arma::vec occ, arma::vec occ_all, arma::vec sumM);
RcppExport SEXP _OccDesign_loglik_cpp(SEXP parsSEXP, SEXP ySEXP, SEXP MSEXP, SEXP X_psiSEXP, SEXP X_pSEXP, SEXP cov_psiSEXP, SEXP cov_pSEXP, SEXP y_occupiedSEXP, SEXP occSEXP, SEXP occ_allSEXP, SEXP sumMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_psi(X_psiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_p(X_pSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cov_psi(cov_psiSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cov_p(cov_pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_occupied(y_occupiedSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type occ(occSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type occ_all(occ_allSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sumM(sumMSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_cpp(pars, y, M, X_psi, X_p, cov_psi, cov_p, y_occupied, occ, occ_all, sumM));
    return rcpp_result_gen;
END_RCPP
}
// gr_loglik_cpp
arma::vec gr_loglik_cpp(arma::vec pars, arma::vec y, arma::vec M, arma::mat X_psi, arma::mat X_p, arma::uvec cov_psi, arma::uvec cov_p, arma::vec y_occupied, arma::vec occ, arma::vec occ_all, arma::vec sumM);
RcppExport SEXP _OccDesign_gr_loglik_cpp(SEXP parsSEXP, SEXP ySEXP, SEXP MSEXP, SEXP X_psiSEXP, SEXP X_pSEXP, SEXP cov_psiSEXP, SEXP cov_pSEXP, SEXP y_occupiedSEXP, SEXP occSEXP, SEXP occ_allSEXP, SEXP sumMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_psi(X_psiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_p(X_pSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cov_psi(cov_psiSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cov_p(cov_pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_occupied(y_occupiedSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type occ(occSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type occ_all(occ_allSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sumM(sumMSEXP);
    rcpp_result_gen = Rcpp::wrap(gr_loglik_cpp(pars, y, M, X_psi, X_p, cov_psi, cov_p, y_occupied, occ, occ_all, sumM));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OccDesign_logistic_CPP", (DL_FUNC) &_OccDesign_logistic_CPP, 1},
    {"_OccDesign_loglik_cpp", (DL_FUNC) &_OccDesign_loglik_cpp, 11},
    {"_OccDesign_gr_loglik_cpp", (DL_FUNC) &_OccDesign_gr_loglik_cpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_OccDesign(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
