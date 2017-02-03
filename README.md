# msap
[![Build Status](https://travis-ci.org/anpefi/msap.png?branch=develop)](https://travis-ci.org/anpefi/msap)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/msap)](http://cran.r-project.org/package=msap)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/msap)](http://cran.r-project.org/package=msap)
  
**An R package for analyzing epigenetic and genetic diversity and differentiation in MSAP assays**

![](images/msap.logo.png)

This package provides a deep analysis of epigenetic variation starting from a binary data matrix indicating the presence or absence of EcoRI-HpaII and EcoRI-MspI fragments, typical of MSAP technique. After comparing the data from both enzyme combinations, the program determines if each fragment is susceptible to methylation (representative of epigenetic variation) or if there is no evidence of methylation (representative of genetic variation). Different analyses of the variation (genetic and epigenetic) among user-defined groups of samples are then performed, as well as the classification of the methylation occurrences in those groups. Statistical testing provide support to the analyses. A comprehensive report of the analyses and several useful plots could help researchers to asses the epigenetic variation in their experiments using MSAP technique.

The package is intended to be easy to use even for those people non-familiar to the R environment. Advanced users could take advantage of available source code to adapt msap for more complex analyses.
 

## Download and Installation
You can install msap automatically from within a R session:


* the latest released version from CRAN [Note that it could be a more recent stable version in GitHub] with

    ```R
    install.packages("msap")   
    ````

* the latest stable version from github with

    ```R
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("anpefi/msap")
 
    ```


## How to use msap
You can find the most recent version of the user's guide vignette in [this link](https://github.com/anpefi/msap/blob/master/pkg/msap/vignettes/msap.pdf).


## Contact

Please use the [issue tracker](https://github.com/anpefi/msap/issues) to report any bug or request any feature.

