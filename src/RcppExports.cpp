// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_IRC_groupcond_logl
double compute_IRC_groupcond_logl(arma::colvec y, arma::mat w, arma::mat z, int N, arma::colvec beta, arma::mat Beta, arma::colvec gamma, arma::mat Gamma, arma::mat Omega, double overdis, arma::colvec nu, arma::mat Lambda, arma::mat Theta, arma::colvec mu_eta, arma::mat Sigma_eta, arma::mat fixeta, arma::vec ghweight, arma::colvec mu_z, arma::mat Sigma_z, arma::mat Sigma_z_lv, String family, bool fixed_z, bool cfa, int const cores);
RcppExport SEXP _lavacreg_compute_IRC_groupcond_logl(SEXP ySEXP, SEXP wSEXP, SEXP zSEXP, SEXP NSEXP, SEXP betaSEXP, SEXP BetaSEXP, SEXP gammaSEXP, SEXP GammaSEXP, SEXP OmegaSEXP, SEXP overdisSEXP, SEXP nuSEXP, SEXP LambdaSEXP, SEXP ThetaSEXP, SEXP mu_etaSEXP, SEXP Sigma_etaSEXP, SEXP fixetaSEXP, SEXP ghweightSEXP, SEXP mu_zSEXP, SEXP Sigma_zSEXP, SEXP Sigma_z_lvSEXP, SEXP familySEXP, SEXP fixed_zSEXP, SEXP cfaSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< double >::type overdis(overdisSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type mu_eta(mu_etaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_eta(Sigma_etaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fixeta(fixetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ghweight(ghweightSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type mu_z(mu_zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_z(Sigma_zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_z_lv(Sigma_z_lvSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< bool >::type fixed_z(fixed_zSEXP);
    Rcpp::traits::input_parameter< bool >::type cfa(cfaSEXP);
    Rcpp::traits::input_parameter< int const >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_IRC_groupcond_logl(y, w, z, N, beta, Beta, gamma, Gamma, Omega, overdis, nu, Lambda, Theta, mu_eta, Sigma_eta, fixeta, ghweight, mu_z, Sigma_z, Sigma_z_lv, family, fixed_z, cfa, cores));
    return rcpp_result_gen;
END_RCPP
}
// compute_groupcond_logl
double compute_groupcond_logl(arma::colvec y, arma::mat w, arma::mat z, int N, arma::colvec beta, arma::mat Beta, arma::colvec gamma, arma::mat Gamma, arma::mat Omega, double overdis, arma::colvec nu, arma::mat Lambda, arma::mat Theta, arma::colvec mu_eta, arma::mat Sigma_eta, arma::mat fixeta, arma::vec ghweight, arma::colvec mu_z, arma::mat Sigma_z, arma::mat Sigma_z_lv, String family, bool fixed_z, bool cfa, int const cores);
RcppExport SEXP _lavacreg_compute_groupcond_logl(SEXP ySEXP, SEXP wSEXP, SEXP zSEXP, SEXP NSEXP, SEXP betaSEXP, SEXP BetaSEXP, SEXP gammaSEXP, SEXP GammaSEXP, SEXP OmegaSEXP, SEXP overdisSEXP, SEXP nuSEXP, SEXP LambdaSEXP, SEXP ThetaSEXP, SEXP mu_etaSEXP, SEXP Sigma_etaSEXP, SEXP fixetaSEXP, SEXP ghweightSEXP, SEXP mu_zSEXP, SEXP Sigma_zSEXP, SEXP Sigma_z_lvSEXP, SEXP familySEXP, SEXP fixed_zSEXP, SEXP cfaSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< double >::type overdis(overdisSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type mu_eta(mu_etaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_eta(Sigma_etaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type fixeta(fixetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ghweight(ghweightSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type mu_z(mu_zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_z(Sigma_zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_z_lv(Sigma_z_lvSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< bool >::type fixed_z(fixed_zSEXP);
    Rcpp::traits::input_parameter< bool >::type cfa(cfaSEXP);
    Rcpp::traits::input_parameter< int const >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_groupcond_logl(y, w, z, N, beta, Beta, gamma, Gamma, Omega, overdis, nu, Lambda, Theta, mu_eta, Sigma_eta, fixeta, ghweight, mu_z, Sigma_z, Sigma_z_lv, family, fixed_z, cfa, cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lavacreg_compute_IRC_groupcond_logl", (DL_FUNC) &_lavacreg_compute_IRC_groupcond_logl, 24},
    {"_lavacreg_compute_groupcond_logl", (DL_FUNC) &_lavacreg_compute_groupcond_logl, 24},
    {NULL, NULL, 0}
};

RcppExport void R_init_lavacreg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
