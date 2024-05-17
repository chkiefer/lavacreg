#define STRICT_R_HEADERS
// [[Rcpp::depends("RcppArmadillo")]]
#include "creg_header_funs.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

//' Compute the group-conditional log-likelihood function
//'
//' @param x Model matrix containing the observed variables in order Y, W, Z
//' @param muy Matrix of conditional expectations of Y given Z and/or latent
// covariates ' @param sigmayw Vector containing the size parameter (for
// overdispersion) and measurement error variances ' @param muwz Matrix of
// conditional expectations of W and Z given latent covariate (if existent) and
// given the group ' @param sigmaz (Conditional) covariance matrix of observed
// covariates Z ' @param ghweight Gauss-Hermite weights, if latent covariates
// included ' @param detvarz Determinant of sigmaz (easier to compute in R...) '
//@param dims Integer vector with dimensional information, in order group size,
// integration points, number of Z, and number of W '
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double compute_groupcond_logl(arma::colvec y, arma::mat w, arma::mat z,
                              arma::colvec beta, arma::mat Beta,
                              arma::colvec gamma, arma::mat Gamma,
                              arma::mat Omega, double overdis, arma::colvec nu,
                              arma::mat Lambda, arma::mat Theta,
                              arma::colvec mu_eta, arma::mat Sigma_eta,
                              arma::mat fixeta, arma::vec ghweight,
                              arma::colvec mu_z, arma::mat Sigma_z,
                              arma::mat Sigma_z_lv, int const cores = 1) {

#if defined(_OPENMP)
  omp_set_num_threads(cores);
#endif

  // If size (i.e., overdispersion parameter with size > 0) is provided
  // the indicator is set to FALSE (i.e., NOT Poisson), otherwise Poisson is
  // used
  bool poisson = overdis ? 0 : 1;

  int N = y.n_elem;
  int N_gh = ghweight.n_elem;
  int no_lv = fixeta.n_cols;
  int no_z = z.n_cols;

  if (!N_gh) {
    N_gh = 1;
  }

  // Rcpp::Rcout << "Initialization worked!\n";
  // sleep(1);

  arma::colvec f_w(N * N_gh, arma::fill::ones);
  arma::mat GH;
  if (no_lv) {
    // Adaption of Gauss-Hermite grid
    arma::mat eig_Lam(no_lv, no_lv, arma::fill::zeros);
    arma::vec eigval;
    arma::mat eig_S;
    eig_sym(eigval, eig_S, Sigma_eta);
    eig_Lam.diag() = arma::sqrt(eigval);

    arma::mat A = eig_S * eig_Lam;

    GH = fixeta * A.t();
    GH.each_row() += mu_eta.as_row();

    // Measurement Model part
    arma::mat LamEta = Lambda * GH.t();

    arma::mat mu_w = nu + LamEta.each_col();

    f_w = dmvnrm_arma_fast(w, mu_w.t(), Theta, false, cores);
  } else {
    // Rcpp::Rcout << "I am here!\n";
    // sleep(1);
    arma::vec tmp(N_gh, arma::fill::ones);
    ghweight = tmp;
    // Rcpp::Rcout << "Value of ghweight : " << tmp << "\n";
    // sleep(1);
  }

  // Stochastic covariate part
  arma::colvec f_z(N * N_gh, arma::fill::ones);
  if (no_lv && no_z) {

    arma::mat GHt = GH.t();

    arma::mat tmp =
        Sigma_z_lv * inv_sympd(Sigma_eta) * (GHt.each_col() - mu_eta);

    arma::mat cmu_z = trans(mu_z + tmp.each_col());

    arma::mat cSigma_z =
        Sigma_z - Sigma_z_lv * inv_sympd(Sigma_eta) * Sigma_z_lv.t();
    if (any(cSigma_z.diag() <= 0)) {
      return (R_PosInf);
    }
    f_z = dmvnrm_arma_fast(z, cmu_z, cSigma_z, false, cores);
  } else if (no_z) {
    f_z = dmvnrm_arma_fast(z, mu_z.t(), Sigma_z, false, cores);
  }

  // Regression part
  // beta
  arma::colvec one(N, arma::fill::ones);
  arma::mat zs = arma::join_rows(one, z);
  // arma::colvec mu_y_z = exp(zs * beta);
  arma::colvec mu_beta = exp(zs * beta);

  // Beta
  arma::colvec mu_Beta(N, arma::fill::ones);
  if (Beta.n_rows) {
    arma::mat tmp = arma::diagvec(z * Beta * z.t());
    mu_Beta = exp(tmp);
  }

  // gamma
  arma::colvec mu_gamma(N_gh, arma::fill::ones);
  if (no_lv) {
    mu_gamma = exp(GH * gamma);
  }

  // Gamma
  arma::colvec mu_Gamma(N, arma::fill::ones);
  if (Gamma.n_rows) {
    arma::mat tmp = arma::diagvec(GH * Gamma * GH.t());
    mu_Gamma = exp(tmp);
  }

  // Omega
  arma::mat mu_Omega(N, N_gh, arma::fill::ones);
  if (Omega.n_rows) {
    arma::mat tmp = z * Omega * GH.t();
    mu_Omega = exp(tmp);
  }

  arma::colvec f_y;
  double alpha;
  if (poisson) {
    // Rcpp::Rcout << "Shortly before Poisson...!\n";
    // sleep(1);
    f_y = dpois_cpp(y, mu_beta, mu_gamma, mu_Beta, mu_Gamma, mu_Omega, cores);
  } else {
    alpha = 1 / overdis;
    f_y = dnegbin_cpp(y, mu_beta, mu_gamma, mu_Beta, mu_Gamma, mu_Omega, alpha,
                      cores);
  }

  // Rcpp::Rcout << "Regression part worked!\n";
  // sleep(1);

  double out = 0.0;
  double tmp = 0.0;

  // Rcpp::Rcout << "Value of ghweight : " << ghweight << "\n";
  // // sleep(5);
  // Rcpp::Rcout << "Value of f_w : " << f_w << "\n";
  // // sleep(5);
  // Rcpp::Rcout << "Value of f_z : " << f_z << "\n";
  // // sleep(5);
  // Rcpp::Rcout << "Value of f_y : " << f_y << "\n";
  // sleep(5);

#pragma omp parallel for schedule(static) private(tmp)                         \
    shared(ghweight, f_w, f_y, f_z) reduction(+ : out)
  for (int i = 0; i < N; i++) {
    tmp = 0.0;
    for (int g = 0; g < N_gh; g++) {
      int index = g + N_gh * i;
      tmp += ghweight(g) * f_w(index) * f_z(index) * f_y(index);
      // Rcpp::Rcout << "Value of ghweight : " << ghweight(g) << "\n";
      // Rcpp::Rcout << "Value of f_w : " << f_w(index) << "\n";
      // Rcpp::Rcout << "Value of f_z : " << f_z(index) << "\n";
      // Rcpp::Rcout << "Value of f_y : " << f_y(index) << "\n";
    }
    out += std::log(tmp);
  }

  // Scaling happens outside
  // out = -1 * out / N;

  // Rcpp::Rcout << "The value of out : " << out << "\n";
  return (out);
}
