# CountReg
CountReg is an R package for fitting count regression models (i.e., Poisson, negative binomial) with manifest as well as latent covariates and within multiple groups. It can be installed via GitHub.

## Installation
`CountReg` is currently not on CRAN. The development version of `CountReg` can be installed directly from this GitHub repository using the additional package `devtools`. Under Windows, please make sure Rtools (http://cran.r-project.org/bin/windows/Rtools) are installed and no odler version of `CountReg` is currently loaded:

```
install.packages("devtools")
library(devtools)

install_github("chkiefer/CountReg")
```

## Run CountReg
The main function of the package is `countreg()`. There is an article available here on GitHub to introduce you to its functionality.