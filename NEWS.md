# lavacreg 0.2-2
* interaction terms between latent and / or manifest predictors can now be included in the model

# lavacreg 0.2-1
* added checks for OpenMP to avoid compilation errors on MacOS
* tests are skipped on CRAN, as minor numerical deviations cannot be avoided

# lavacreg 0.2-0

* refactoring of the whole codebase
* major change to the C++ files; likelihood estimation is now based on RcppArmadillo and matrix algebra
* major changes to the R files;  partable, modellist and so on were adapted to the matrix-based computations
* codebase now supports interaction terms, but this feature will be released with the next minor version (not tested enough right now)

# lavacreg 0.1-2

* minor changes to the C++ files due to upcoming changes in Rcpp package (STRICT_R_HEADERS)
* added some details in the documentation of the main function

# lavacreg 0.1-1

* skipping two tests failing on MacOS (i.e. applications in which only one latent covariate is considered)
* lowered numerical precision for test on Poisson regression with three manifest covariates (for Solaris)
* added a more detailed summary-function

# lavacreg 0.1-0

* Added a `NEWS.md` file to track changes to the package.

