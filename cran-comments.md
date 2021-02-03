## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)
- R-hub linux-x86_64-rocker-gcc-san (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE

  Maintainer: 'Christoph Kiefer <christoph.kiefer@uni-bielefeld.de>'
  
  Re-submission: I addressed your three requests:
  1. A reference was added to the description field of the DESCRIPTION
  2. documentation for functions is.count and summary was added. function creg_partable is no longer exported
  3. replaced dontrun with donttest for the examples, which run for 2-3 minutes
  
  The remaining note concerns that I submit the package to CRAN for the first time.


0 errors √ | 0 warnings √ | 1 note x