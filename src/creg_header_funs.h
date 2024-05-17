#ifndef CREG_HEADER_FUNS
#define CREG_HEADER_FUNS

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

static double const log2pi = std::log(2.0 * M_PI);
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);
arma::vec dmvnrm_arma_fast(arma::mat const &x, arma::mat const &mean,
                           arma::mat const &sigma, bool const logd,
                           int const cores);
arma::colvec dpois_cpp(arma::colvec y, arma::colvec mu_beta,
                       arma::colvec mu_gamma, arma::colvec mu_Beta,
                       arma::colvec mu_Gamma, arma::mat mu_Omega,
                       int const cores);
arma::colvec dnegbin_cpp(arma::colvec y, arma::colvec mu_beta,
                         arma::colvec mu_gamma, arma::colvec mu_Beta,
                         arma::colvec mu_Gamma, arma::mat mu_Omega,
                         double alpha, int const cores);
#endif
