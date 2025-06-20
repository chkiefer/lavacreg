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
//' covariates ' @param sigmayw Vector containing the size parameter (for
//' overdispersion) and measurement error variances
//' @param muwz Matrix of conditional expectations of W and Z given latent
//' covariate (if existent) and given the group
//' @param sigmaz (Conditional) covariance matrix of observed covariates Z
//' @param ghweight Gauss-Hermite weights, if latent covariates included
//' @param detvarz Determinant of sigmaz (easier to compute in R...)
//' @param dims Integer vector with dimensional information, in order group
//' size, integration points, number of Z, and number of W
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double compute_IRC_groupcond_logl(
    arma::colvec y, arma::mat w, arma::mat z, int N, arma::colvec beta,
    arma::mat Beta, arma::colvec gamma, arma::mat Gamma, arma::mat Omega,
    double overdis, arma::colvec nu, arma::mat Lambda, arma::mat Theta,
    arma::colvec mu_eta, arma::mat Sigma_eta, arma::mat fixeta,
    arma::vec ghweight, arma::colvec mu_z, arma::mat Sigma_z,
    arma::mat Sigma_z_lv, String family, bool fixed_z, bool fst = 1,
    bool sst = 1, int const cores = 1) {

#if defined(_OPENMP)
  omp_set_num_threads(cores);
#endif

  // Set to 1 if verbose output is required (only for debugging)
  bool verbose = 1;

  // int N = y.n_elem;
  int N_gh = ghweight.n_elem;
  int no_lv = fixeta.n_cols;
  int no_z = z.n_cols;
  int no_w = w.n_cols;

  // Set N_gh to 1 if no intpoints are present
  if (N_gh == 0) {
    N_gh = 1;
  }

  // ************************************ //
  // Predictor Model part (or first step) //
  // ************************************ //
  arma::colvec f_wz(N, arma::fill::ones);
  arma::colvec mu_implied;
  arma::mat Sigma_implied;
  arma::mat predobs;
  arma::colvec mu_pred;
  arma::mat Sigma_pred;
  arma::mat Sigma_pred_lv;

  // Distinguish three cases:
  // 1. only LVs (CFA)
  // 2. only stochastic Zs
  // both LVs and stochastic Zs
  if (no_lv && !(no_z && !fixed_z)) {
    // 1. only LVs
    // Compute implied mu and Sigma for indicators
    mu_implied = nu + Lambda * mu_eta;
    Sigma_implied = Lambda * Sigma_eta * Lambda.t() + Theta;

    if (fst == 1) {
      f_wz = dmvnrm_arma_fast(w, mu_implied.t(), Sigma_implied, true, cores);
    }
    //
    predobs = w;
    mu_pred = mu_implied;
    Sigma_pred = Sigma_implied;
    Sigma_pred_lv = Lambda * Sigma_eta;
  } else if (!no_lv && (no_z && !fixed_z)) {
    // 2. Stochastic covariate only
    // direct computation of density
    if (fst == 1) {
      f_wz = dmvnrm_arma_fast(z, mu_z.t(), Sigma_z, true, cores);
    }
    //
  } else if (no_lv && (no_z && !fixed_z)) {
    // 3. LVs and stochastic covariates

    // build joint (implied) mean vector and covariance matrix
    mu_implied = nu + Lambda * mu_eta;
    Sigma_implied = Lambda * Sigma_eta * Lambda.t() + Theta;
    arma::mat Sigma_z_w = Sigma_z_lv * Lambda.t();

    int dim_wz = no_w + no_z;
    arma::mat Sigma_wz(dim_wz, dim_wz, arma::fill::zeros);

    Sigma_wz.submat(0, 0, (no_w - 1), (no_w - 1)) = Sigma_implied;
    Sigma_wz.submat(no_w, no_w, (dim_wz - 1), (dim_wz - 1)) = Sigma_z;
    Sigma_wz.submat(0, no_w, (no_w - 1), (dim_wz - 1)) = Sigma_z_w.t();
    Sigma_wz.submat(no_w, 0, (dim_wz - 1), (no_w - 1)) = Sigma_z_w;
    arma::mat wz = join_rows(w, z);
    arma::colvec mu_wz = join_cols(mu_implied, mu_z);

    // final computations
    if (fst == 1) {
      f_wz = dmvnrm_arma_fast(wz, mu_wz.t(), Sigma_wz, true, cores);
    }
    //
    predobs = wz;
    mu_pred = mu_wz;
    Sigma_pred = Sigma_wz;
    Sigma_pred_lv = join_cols(Lambda * Sigma_eta, Sigma_z_lv);
  }

  // if cfa or first step: return here
  if (sst == 0) {
    return (sum(f_wz));
  }

  // ****************************** //
  // Adaption of Gauss-Hermite grid //
  // ****************************** //
  arma::cube GH_i(N_gh, no_lv, N);
  if (no_lv) {
    arma::mat sig11 = Sigma_eta;
    arma::mat sig22_inv = inv(Sigma_pred);
    arma::mat sig21 = Sigma_pred_lv;
    arma::mat sig_bar = sig11 - sig21.t() * sig22_inv * sig21;

    arma::mat eig_Lam(no_lv, no_lv, arma::fill::zeros);
    arma::vec eigval;
    arma::mat eig_S;
    eig_sym(eigval, eig_S, sig_bar);
    eig_Lam.diag() = arma::sqrt(eigval);

    arma::mat A = eig_S * eig_Lam;
    arma::mat GH = fixeta * A.t();

    // Computation of conditional expectation of eta given W and Z
    arma::mat sig_part =
        sig21.t() * sig22_inv * (predobs.each_row() - mu_pred.t()).t();
    arma::mat mu_bar = mu_eta + sig_part.each_col();
    mu_bar = mu_bar.t();

    for (int i = 0; i < N; i++) {
      arma::mat GH_tmp = GH;
      GH_tmp.each_row() += mu_bar.row(i);
      GH_i.slice(i) = GH_tmp;
    }
  } else {
    arma::vec tmp(N_gh, arma::fill::ones);
    ghweight = tmp;
  }

  // ************************************ //
  // Outcome Model part (or second step) //
  // ************************************ //
  // beta
  arma::colvec one(N, arma::fill::ones);
  arma::mat zs = arma::join_rows(one, z);
  arma::colvec mu_beta = exp(zs * beta);

  // Beta
  arma::colvec mu_Beta(N, arma::fill::ones);
  if (Beta.n_rows) {
    // ersetzen wie unten?
    arma::mat tmp = arma::diagvec(z * Beta * z.t());
    mu_Beta = exp(tmp);
  }

  // gamma
  arma::mat mu_gamma(N, N_gh, arma::fill::ones);
  if (no_lv) {
    for (int i = 0; i < N; i++) {
      arma::mat tmp = GH_i.slice(i) * gamma;
      mu_gamma.row(i) = arma::exp(tmp.t());
    }
  }

  // Gamma
  // is that correct? Does Gamma have to be a matrix?
  arma::mat mu_Gamma(N, N_gh, arma::fill::ones);
  if (Gamma.n_rows) {
    // arma::mat tmp = arma::diagvec(GH * Gamma * GH.t());
    for (int i = 0; i < N; i++) {
      arma::mat tmp = sum((GH_i.slice(i) * Gamma) % GH_i.slice(i), 1);
      mu_Gamma.row(i) = exp(tmp.t());
    }
  }

  // Omega
  arma::mat mu_Omega(N, N_gh, arma::fill::ones);
  if (Omega.n_rows) {
    for (int i = 1; i < N; i++) {
      // still to fix
      arma::mat GH_i_tmp = GH_i.slice(i);
      arma::mat tmp = z * Omega * GH_i_tmp.t();
      mu_Omega = exp(tmp);
    }
  }

  arma::colvec f_y(N * N_gh, arma::fill::ones);
  double alpha;
  if (family == "poisson") {
    f_y =
        dpois_IRC_cpp(y, mu_beta, mu_gamma, mu_Beta, mu_Gamma, mu_Omega, cores);
  } else if (family == "nbinom") {
    alpha = 1 / overdis;
    f_y = dnegbin_IRC_cpp(y, mu_beta, mu_gamma, mu_Beta, mu_Gamma, mu_Omega,
                          alpha, cores);
  } else if (family == "logistic") {
    f_y = dbinom_logis_IRC_cpp(y, mu_beta, mu_gamma, mu_Beta, mu_Gamma,
                               mu_Omega, cores);
  }

  double out = 0.0;
  double tmp = 0.0;

#pragma omp parallel for schedule(static) private(tmp) shared(ghweight, f_y)   \
    reduction(+ : out)
  for (int i = 0; i < N; i++) {
    tmp = 0.0;
    for (int g = 0; g < N_gh; g++) {
      int index = g + N_gh * i;
      tmp += ghweight(g) * f_y(index);
    }
    if (tmp == 0.0) {
      tmp = 1.7e-308; // close to minimum on my laptop
      // Rcpp::Rcout << "Had to correct a value of tmp to : " << tmp << "\n";
    }
    out += std::log(tmp);
  }

  // Scaling happens outside
  // out = -1 * out / N;

  // this is wrong, if first step is omitted!
  if ((no_z && !fixed_z) || no_lv) {
    out += sum(f_wz);
  }

  return (out);
}
