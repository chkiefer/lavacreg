#include <Rcpp.h>
using namespace Rcpp;

//' Compute the group-conditional log-likelihood function
//' 
//' @param x Model matrix containing the observed variables in order Y, W, Z
//' @param muy Matrix of conditional expectations of Y given Z and/or latent covariates
//' @param sigmayw Vector containing the size parameter (for overdispersion) and measurement error variances 
//' @param muwz Matrix of conditional expectations of W and Z given latent covariate (if existent) and given the group
//' @param sigmaz (Conditional) covariance matrix of observed covariates Z
//' @param ghweight Gauss-Hermite weights, if latent covariates included
//' @param detvarz Determinant of sigmaz (easier to compute in R...)
//' @param dims Integer vector with dimensional information, in order group size, integration points, number of Z, and number of W
//' @export
// [[Rcpp::export]]
double compute_groupcond_logl(NumericMatrix x, 
                              NumericMatrix muy, 
                              NumericVector sigmayw, 
                              NumericMatrix muwz, 
                              NumericMatrix sigmaz, 
                              NumericVector ghweight,
                              double detvarz,
                              IntegerVector dims){
  // If size (i.e., overdispersion parameter with size > 0) is provided
  // the indicator is set to FALSE (i.e., NOT Poisson), otherwise Poisson is used 
  bool poisson = sigmayw[0] ? 0 : 1;
  double size;
  double isize;
  double gamisize;
  if (!poisson){
    size = 1/sigmayw[0];
    isize = 1/size;
    gamisize = tgamma(isize);
  }
  
  
  
  // loop dimensions
  int N = dims[0];
  int ghpoints = dims[1];
  int zn = dims[2];
  int wn = dims[3];
  
  // constants
  double cz = pow(sqrt(2*PI), zn); 
  double cw = pow(sqrt(2*PI), wn); 
  
  // If latent variables present, precompute determinant for veps-matrix 
  // (note: measurement errors are independent, multiplicate diagonal elements)
  // else it is 1 and does not hurt in later multiplications
  double detveps = 1;
  if (wn){
    for (int v=0; v < wn; ++v){
      detveps *= sigmayw[v + 1];
    }
  }

  // memory allocation
  double isum;
  double expsumw;
  double expsumz;
  double prod;
  
  // Loops for computing densities and summing things up
  double out = 0;
  
  for (int i = 0; i < N; ++i){
    // Outer sum of log-likelihood (sums over individual observations)
    isum = 0;
    
    for (int g = 0; g < ghpoints; ++g){
      // Inner sum for Gauss-Hermite quadrature, if latent variables present
      // else: only 1 summand/point, weight = 1
      expsumz = 0;
      expsumw = 0;
      prod = 1;
      double y0 = x(i, 0);
      double mu0 = muy(g, i);

      // Product starts with GH-weight and choosen count density
      if (poisson){
        // Poisson density
        prod = ghweight[g] * pow(mu0, y0) * exp(-mu0) / tgamma(y0 + 1);
      } else {
        // Negative binomial density
        // Note: parts of density drawn to outer sum to avoid repeated computation
        double musize = 1/(mu0*size + 1);
        prod = ghweight[g] * tgamma(y0 + isize)/(tgamma(y0 + 1) * gamisize) * pow(musize, isize) * pow(1-musize, y0);
        //prod = ghweight[g] * tgamma(y0 + isize)/(tgamma(y0 + 1) * gamisize) * pow(musize, isize) * pow(1-musize, y0);
      }
      
      // Product for latent variable indicators W
      for (int v = 0; v < wn; ++v){
        //inner product over normally distributed indicator variables
        double w0 = x(i, v + 1);
        double mu0 = muwz(g, v);
        double muw = w0-mu0;
        double sig0 = sigmayw[v + 1];
        expsumw += (muw*muw)/(sig0);
      }
      
      // Product for manifest covariates Z
      for (int v = 0; v < zn; ++v){
        for (int w = 0; w < zn; ++ w){
          //inner product over normally distributed indicator variables
          double z0 = x(i, v+wn+1);
          double z1 = x(i, w+wn+1);
          double mu0 = muwz(g, v+wn);
          double mu1 = muwz(g, w+wn);
          double muz0 = z0-mu0;
          double muz1 = z1-mu1;
          double isig = sigmaz(v,w);
          expsumz += (muz0*muz1*isig);
        }
      }
      
      isum += prod * exp(-expsumw/2) * exp(-expsumz/2);
    }
    // if (!poisson){
    //   isum = isum/gamisize;
    // }
    out += log(isum);
  }
  
  return out - N*log(cz*sqrt(detveps)*cw*sqrt(detvarz));
}