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
    arma::mat Sigma_z_lv, String family, bool fixed_z, bool cfa,
    int const cores = 1) {

#if defined(_OPENMP)
  omp_set_num_threads(cores);
#endif

  bool verbose = 0;

  // If size (i.e., overdispersion parameter with size > 0) is provided
  // the indicator is set to FALSE (i.e., NOT Poisson), otherwise Poisson is
  // used
  // bool poisson = overdis ? 0 : 1;
  // replaced by the family argument!

  // int N = y.n_elem;
  int N_gh = ghweight.n_elem;
  int no_lv = fixeta.n_cols;
  int no_z = z.n_cols;
  int no_w = w.n_cols;

  if (verbose) {
    Rcpp::Rcout << "Value of N : " << N << "\n";
    Rcpp::Rcout << "Value of N_gh : " << N_gh << "\n";
    Rcpp::Rcout << "Value of no_lv : " << no_lv << "\n";
    Rcpp::Rcout << "Value of no_z : " << no_z << "\n";
    Rcpp::Rcout << "Value of no_w : " << no_w << "\n";
  }

  if (!N_gh) {
    N_gh = 1;
  }

  // CFA part (or observed covariates and indicators part)
  arma::colvec f_wz(N, arma::fill::ones);
  arma::colvec mu_implied;
  arma::mat Sigma_implied;
  // if no stochastic covariate or for CFA:
  if (no_lv && (!no_z || fixed_z) || cfa) {
    // Print comment for verbose mode
    if (verbose) {
      Rcpp::Rcout << "First condition (only LVs) valid! \n";
    }

    // Computations
    mu_implied = nu + Lambda * mu_eta;
    Sigma_implied = Lambda * Sigma_eta * Lambda.t() + Theta;
    f_wz = dmvnrm_arma_fast(w, mu_implied.t(), Sigma_implied, true, cores);

    if (verbose) {
      Rcpp::Rcout << "First Value of f_wz : " << f_wz(1) << "\n";
    }
    if (cfa) {
      if (verbose) {
        Rcpp::Rcout << "Value of mu_implied : " << mu_implied << "\n";
        Rcpp::Rcout << "Value of Sigma_implied : " << Sigma_implied << "\n";
      }

      return (sum(f_wz));
    }
  } else if (no_lv && no_z && !fixed_z) {
    // Print comment for verbose mode
    if (verbose) {
      Rcpp::Rcout
          << "Second condition (LVs and stochastic predictors) valid! : \n";
    }

    // hier muss nun extra die gemeinsame Kovarianzmatrix gebaut werden
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

    if (verbose) {
      Rcpp::Rcout << "Value of mu_implied : " << mu_implied << "\n";
      Rcpp::Rcout << "Value of Sigma_implied : " << Sigma_implied << "\n";
      Rcpp::Rcout << "Value of mu_wz : " << mu_wz << "\n";
      Rcpp::Rcout << "Value of Sigma_wz : " << Sigma_wz << "\n";
    }

    // final computations
    f_wz = dmvnrm_arma_fast(wz, mu_wz.t(), Sigma_wz, true, cores);

  } else if (no_z && !fixed_z) {
    // Stochastic covariate only part
    if (verbose) {
      Rcpp::Rcout << "Third condition (stochastic predictors only) valid! : \n";
    }
    f_wz = dmvnrm_arma_fast(z, mu_z.t(), Sigma_z, true, cores);
  }

  if (verbose) {
    Rcpp::Rcout << "First Value of f_wz : " << f_wz(1) << "\n";
  }

  arma::cube GH_i(N_gh, no_lv, N);
  if (no_lv) {
    if (verbose) {
      Rcpp::Rcout << "Adapation of Gauss-Hermite for LVs.\n";
    }
    // Adaption of Gauss-Hermite grid
    arma::mat sig11 = Sigma_eta;
    arma::mat sig22_inv = inv(Sigma_implied);
    arma::mat sig21 = Lambda * sig11;
    arma::mat sig_bar = sig11 - sig21.t() * sig22_inv * sig21;

    if (verbose) {
      Rcpp::Rcout << "Value of Theta : " << Theta << "\n";
      Rcpp::Rcout << "Value of Sigma_eta : " << Sigma_eta << "\n";
      Rcpp::Rcout << "Value of Sigma_implied : " << Sigma_implied << "\n";
      Rcpp::Rcout << "Value of Lambda : " << Lambda << "\n";
      Rcpp::Rcout << "Value of sig_bar : " << sig_bar << "\n";
      Rcpp::Rcout << "Sig_bar worked.\n";
    }

    arma::mat eig_Lam(no_lv, no_lv, arma::fill::zeros);
    arma::vec eigval;
    arma::mat eig_S;
    eig_sym(eigval, eig_S, sig_bar);
    eig_Lam.diag() = arma::sqrt(eigval);

    if (verbose) {
      Rcpp::Rcout << "eig_Lam worked.\n";
    }

    arma::mat A = eig_S * eig_Lam;
    arma::mat GH = fixeta * A.t();

    if (verbose) {
      Rcpp::Rcout << "Value of mu_eta : " << size(mu_eta) << "\n";
      Rcpp::Rcout << "Value of sig21 : " << size(sig21) << "\n";
      Rcpp::Rcout << "Value of sig22_inv : " << size(sig22_inv) << "\n";
      Rcpp::Rcout << "Value of w : " << size(w) << "\n";
      Rcpp::Rcout << "Value of mu_implied : " << size(mu_implied) << "\n";
    }

    // Computation of conditional expectation of eta given W and Z
    arma::mat sig_part =
        sig21.t() * sig22_inv * (w.each_row() - mu_implied.t()).t();
    arma::mat mu_bar = mu_eta + sig_part.each_col();
    mu_bar = mu_bar.t();

    if (verbose) {
      Rcpp::Rcout << "mu_bar worked.\n";
    }

    if (verbose) {
      Rcpp::Rcout << "Value of GH_i : " << size(GH_i) << "\n";
    }
    for (int i = 0; i < N; i++) {
      arma::mat GH_tmp = GH;
      GH_tmp.each_row() += mu_bar.row(i);
      GH_i.slice(i) = GH_tmp;
      // if (verbose) {
      //   Rcpp::Rcout << "Value of GH_i_slice : " << GH_i.slice(i) << "\n";
      //   Rcpp::Rcout << "Value of mu_bar : " << mu_bar.row(i) << "\n";
      //   Rcpp::Rcout << "Value of GH : " << GH << "\n";
      // }
    }
    if (verbose) {
      Rcpp::Rcout << "Value of GH_i : " << size(GH_i) << "\n";
    }

  } else {
    arma::vec tmp(N_gh, arma::fill::ones);
    ghweight = tmp;
  }

  if (verbose) {
    Rcpp::Rcout << "Value of GH_i : " << size(GH_i) << "\n";
  }

  // Regression part

  if (verbose) {
    Rcpp::Rcout << "Regression part starts!\n";
  }
  // beta
  arma::colvec one(N, arma::fill::ones);
  arma::mat zs = arma::join_rows(one, z);
  arma::colvec mu_beta = exp(zs * beta);

  if (verbose) {
    Rcpp::Rcout << "beta worked!\n";
  }

  // Beta
  arma::colvec mu_Beta(N, arma::fill::ones);
  if (Beta.n_rows) {
    // ersetzen wie unten?
    arma::mat tmp = arma::diagvec(z * Beta * z.t());
    mu_Beta = exp(tmp);
  }

  if (verbose) {
    Rcpp::Rcout << "Beta worked!\n";
  }

  // gamma
  arma::mat mu_gamma(N, N_gh, arma::fill::ones);
  if (no_lv) {
    if (verbose) {
      Rcpp::Rcout << "Value of gamma : " << gamma << "\n";
    }
    for (int i = 0; i < N; i++) {
      arma::mat tmp = GH_i.slice(i) * gamma;
      mu_gamma.row(i) = arma::exp(tmp.t());
    }
  }

  if (verbose) {
    // Rcpp::Rcout << "Value of mu_gamma : " << mu_gamma << "\n";
    Rcpp::Rcout << "gamma worked!\n";
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

  if (verbose) {
    Rcpp::Rcout << "Gamma worked!\n";
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
  if (verbose) {
    Rcpp::Rcout << "Omega worked!\n";
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
    if (verbose) {
      Rcpp::Rcout << "First Value of f_y : " << f_y(1) << "\n";
    }
  } else if (family == "logistic") {
    f_y = dbinom_logis_IRC_cpp(y, mu_beta, mu_gamma, mu_Beta, mu_Gamma,
                               mu_Omega, cores);
  }

  if (verbose) {
    Rcpp::Rcout << "Regression part worked!\n";
    Rcpp::Rcout << "size of double:" << sizeof(double) << "\n";
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
  if (verbose) {
    Rcpp::Rcout << "The value of sum(f_y) : " << sum(f_y) << "\n";
    Rcpp::Rcout << "The value of out : " << out << "\n";
  }

  if ((no_z && !fixed_z) || no_lv) {
    out += sum(f_wz);
  }

  if (verbose) {
    Rcpp::Rcout << "The final value of out : " << out << "\n";
  }

  return (out);
}
