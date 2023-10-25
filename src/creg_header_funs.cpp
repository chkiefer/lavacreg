// [[Rcpp::depends("RcppArmadillo")]]
#include "creg_header_funs.h"
#include <RcppArmadillo.h>
#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

// static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat) {
  arma::uword const n = trimat.n_cols;

  for (unsigned j = n; j-- > 0;) {
    double tmp(0.);
    for (unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// https://gallery.rcpp.org/articles/dmvnorm_arma/
arma::vec dmvnrm_arma_fast(arma::mat const &x, arma::mat const &mean,
                           arma::mat const &sigma, bool const logd = false,
                           int const cores = 1) {
  using arma::uword;

#if defined(_OPENMP)
  omp_set_num_threads(cores);
#endif

  uword const n = x.n_rows, xdim = x.n_cols;
  uword const n_gh = mean.n_rows;
  arma::vec out(n * n_gh);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
               constants = -(double)xdim / 2.0 * log2pi,
               other_terms = rootisum + constants;

  arma::rowvec z;
#pragma omp parallel for schedule(static) private(z)
  for (uword i = 0; i < n; i++) {
    for (uword g = 0; g < n_gh; g++) {
      z = (x.row(i) - mean.row(g));
      inplace_tri_mat_mult(z, rooti);
      out(g + i * n_gh) = other_terms - 0.5 * arma::dot(z, z);
    }
  }

  if (logd)
    return out;
  return exp(out);
}

// Poisson Density
arma::colvec dpois_cpp(arma::colvec y, arma::colvec mu_beta,
                       arma::colvec mu_Beta, arma::colvec mu_y_lv,
                       int const cores = 1) {

#if defined(_OPENMP)
  omp_set_num_threads(cores);
#endif

  // Poisson density
  int N = y.n_elem;
  int n_gh = mu_y_lv.n_elem;
  arma::colvec out(N * n_gh);

  double lambda;
  double yi, mu_y_z_i;
  double lgamma_yi_1;

  // Rcpp::Rcout << "...inside Poisson\n";
  // sleep(1);
#pragma omp parallel for schedule(static) private(lambda, yi, mu_y_z_i,        \
                                                      lgamma_yi_1) shared(y)
  for (int i = 0; i < N; i++) {
    mu_y_z_i = mu_beta(i) * mu_Beta(i);
    yi = y(i);
    lgamma_yi_1 = std::lgamma(yi + 1);
    for (int g = 0; g < n_gh; g++) {
      lambda = mu_y_z_i * mu_y_lv(g);
      out(g + i * n_gh) = yi * std::log(lambda) - lambda - lgamma_yi_1;
    }
  }
  // Rcpp::Rcout << "End of Poisson...\n";
  // sleep(1);
  return arma::exp(out);
}

// Negative Binomial Density
arma::colvec dnegbin_cpp(arma::colvec y, arma::colvec mu_y_z,
                         arma::colvec mu_y_lv, double alpha,
                         int const cores = 1) {

#if defined(_OPENMP)
  omp_set_num_threads(cores);
#endif

  // NegBin density (Hilbe, 2011, p.190)
  int N = y.n_elem;
  int n_gh = mu_y_lv.n_elem;
  double lambda;
  arma::colvec out(N * n_gh);
  double yi, mu_y_z_i;
  double lgamma_yi_1;
  double lgamma_yi_alpha;

  double inv_alpha = 1 / alpha;
  double lgamma_inv_alpha = std::lgamma(inv_alpha);

#pragma omp parallel for schedule(static) private(                             \
        lambda, yi, mu_y_z_i, lgamma_yi_1, lgamma_yi_alpha)                    \
    shared(y, inv_alpha, lgamma_inv_alpha)
  for (int i = 0; i < N; i++) {
    mu_y_z_i = mu_y_z(i);
    yi = y(i);
    lgamma_yi_1 = std::lgamma(yi + 1);
    lgamma_yi_alpha = std::lgamma(yi + inv_alpha);
    for (int g = 0; g < n_gh; g++) {
      lambda = mu_y_z_i * mu_y_lv(g);
      double alpha_x_lambda = alpha * lambda;
      out(g + i * n_gh) = yi * std::log(alpha_x_lambda / (1 + alpha_x_lambda)) -
                          inv_alpha * std::log(1 + alpha_x_lambda) +
                          lgamma_yi_alpha - lgamma_yi_1 - lgamma_inv_alpha;
    }
  }
  return arma::exp(out);
}
