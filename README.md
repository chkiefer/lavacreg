# lavacreg: Latent Variable Count Regression Models
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/chkiefer/lavacreg.svg?branch=master)](https://travis-ci.com/chkiefer/lavacreg)

lavacreg is an R package for fitting count regression models (i.e., Poisson, negative binomial) with manifest as well as latent covariates and within multiple groups. It can be installed via GitHub.

## Installation
`lavacreg` is currently on CRAN. The development version of `lavacreg` can be installed directly from this GitHub repository using the additional package `devtools`. Under Windows, please make sure Rtools (http://cran.r-project.org/bin/windows/Rtools) are installed and no older version of `lavacreg` is currently loaded:

```
install.packages("devtools")
library(devtools)

install_github("chkiefer/lavacreg")
```

## Run lavacreg
The main function of the package is `countreg()`. There is an article available here on GitHub (https://chkiefer.github.io/lavacreg/articles/intro.html) to introduce you to its functionality.
