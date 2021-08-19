#define STRICT_R_HEADERS
#include <Rcpp.h>
using namespace Rcpp;

//' Compute the group-conditional log-likelihood function
//' 
//' @param datalist list containing group-specific datasets
//' @param modellist list containing group-specific informations to compute likelihood
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector creg_group_logl_cpp(List datalist, List modellist){
  int Ng = datalist.size();
  NumericVector res(Ng);
  
  for (int group = 0; group < Ng; ++group){
    List modellist_g = as<List>(modellist[group]);
    
    
    // x Model matrix containing the observed variables in order Y, W, Z
    // muy Matrix of conditional expectations of Y given Z and/or latent covariates
    // sigmayw Vector containing the size parameter (for overdispersion) and measurement error variances 
    // muwz Matrix of conditional expectations of W and Z given latent covariate (if existent) and given the group
    // sigmaz (Conditional) covariance matrix of observed covariates Z
    // ghweight Gauss-Hermite weights, if latent covariates included
    // detvarz Determinant of sigmaz (easier to compute in R...)
    // dims Integer vector with dimensional information, in order group size, integration points, number of Z, and number of W
    
    NumericMatrix muy = as<NumericMatrix>(modellist_g["muy"]);
    NumericVector sigmayw = as<NumericVector>(modellist_g["sigmayw"]);
    NumericMatrix muwz = as<NumericMatrix>(modellist_g["muwz"]);
    NumericMatrix sigmaz = as<NumericMatrix>(modellist_g["sigmaz"]);
    NumericVector ghweight = as<NumericVector>(modellist_g["ghweight"]);
    double detvarz = as<double>(modellist_g["detvarz"]);
    IntegerVector dims = as<IntegerVector>(modellist_g["dims"]);
    
    NumericMatrix x = as<NumericMatrix>(datalist[group]);
        
        // If size (i.e., overdispersion parameter with size > 0) is provided
        // the indicator is set to FALSE (i.e., NOT Poisson), otherwise Poisson is used 
        bool poisson = sigmayw[0] ? 0 : 1;
        double size;
        double isize;
        double gamisize = 1.0;
        if (!poisson){
          size = 1/sigmayw[0];
          isize = 1/size;
          gamisize = std::tgamma(isize);
        }
        
        // loop dimensions
        int N = dims[0];
        int ghpoints = dims[1];
        int zn = dims[2];
        int wn = dims[3];
        
        // constants
        double sqrtpi = std::sqrt(2*M_PI);
        double cz = std::pow(sqrtpi, zn);
        double cw = std::pow(sqrtpi, wn);
        
        // If latent variables present, precompute determinant for veps-matrix 
        // (note: measurement errors are independent, multiplicate diagonal elements)
        // else it is 1 and does not hurt in later multiplications
        double detveps = 1.0;
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
        double out = 0.0;
        
        for (int i = 0; i < N; ++i){
          // Outer sum of log-likelihood (sums over individual observations)
          isum = 0.0;
          
          for (int g = 0; g < ghpoints; ++g){
            // Inner sum for Gauss-Hermite quadrature, if latent variables present
            // else: only 1 summand/point, weight = 1
            expsumz = 0.0;
            expsumw = 0.0;
            prod = 0.0;
            double y0 = x(i, 0);
            double mu0 = muy(g, i);
            
            // Product starts with GH-weight and choosen count density
            if (poisson){
              // Poisson density
              // prod = ghweight[g] * std::pow(mu0, y0) * exp(-mu0) / std::tgamma(y0 + 1);
              prod = std::log(ghweight[g]) + y0 * mu0 - std::exp(mu0) - std::lgamma(y0 + 1) ;
            } else {
              // Negative binomial density
              // Note: parts of density drawn to outer sum to avoid repeated computation
              
              // Version with differences between 32 and 64bit systems
              // double musize = 1/(mu0*size + 1);
              // prod = ghweight[g] * std::tgamma(y0 + isize)/std::tgamma(y0 + 1) * std::pow(musize, isize) * std::pow(1-musize, y0);
              
              // Test of numerical stable versions
              double musize = std::exp(mu0)*size + 1;
              prod = std::log(ghweight[g]) + y0 * std::log(1-1/musize) - isize * std::log(musize) + std::lgamma(y0 + isize) - std::lgamma(y0 + 1);
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
              for (int w = 0; w < zn; ++w){
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
            
            //isum +=  prod * std::exp(-expsumw/2) * std::exp(-expsumz/2);
            isum += std::exp(prod - expsumw/2 - expsumz/2);
          }
          
          out += std::log(isum);
        }
        
        res[group] = out - N*std::log(cz*std::sqrt(detveps)*cw*std::sqrt(detvarz)*gamisize);
        
  }
  return res;
}
