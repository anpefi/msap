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


## How to use msap (QUICK GUIDE)

For a more detailed user guide of msap, please check this [vignette](vignettes/msap.pdf).

### Preparation of datafile

In order to use *msap* for analyzing your results from a MSAP experiment, you need to provide a CVS data file ([There is an example here](vignettes/example.csv)) with a binary matrix indicating the presence (coded as 1) or absence (coded as 0) of *Eco*RI-*Hpa*II and *Eco*RI-*Msp*I fragments in a bunch of samples of two or more populations/groups. 

The first row should include the markers name/references, these should be ordered by primer combination (if applicable) as analysis will be seperated for them. The **first column** should provide the label for the group where the sample is included, with the aim to make comparisons between different groups. **Second column** is reserved for an arbitrary label (i.e. to name the sample). **Third column** should identify the isoschizomer with ’HPA’ or ’MSP’ (**Note that both rows for each sample should be present**). The remaining columns are for the different MSAP markers. 

### Run msap

Once we are in the right working directory with an appropiate data file, we can run all analyses of *msap* with a single command (change “example.csv” by the name of your datafile, keeping the quotes, and change “Example” by a custom name, keeping quotes, to identify your data):
```r
library(msap)
myList <- msap("example.csv", name = "Example")
```
This will run all analyses, with default options, and will show an on screen text report with the results as well as returning a list called myList (or whatever you put before the <-) with useful data for further analyses. **Please refer to the [vignette](vignettes/msap.pdf) for further information about options and the different output files.**


## Contact

Please use the [issue tracker](https://github.com/anpefi/msap/issues) to report any bug or request any feature.

