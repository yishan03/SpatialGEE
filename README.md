# SpatialGEE

R package for spatial data analysis using generalized estimating equations (GEEs). This package enables robust statistical analysis of spatially correlated data and provides tools for hypothesis testing using Wald and generalized score tests.

## Reference

Wang, Y., Zang, C., Li, Z., Guo, C. C., Lai, D., & Wei, P. (2024). A comparative study of statistical methods for identifying differentially expressed genes in spatial transcriptomics. Unpublished manuscript.

## Installation

Install the SpatialGEE package via GitHub:

```r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("yishan03/SpatialGEE", build_vignettes = TRUE)
```

## To get started with the package

Load the package and open the package vignette:

```r
library("SpatialGEE")
vignette("SpatialGEE-Tutorial")
```

There will be examples for all the functions in SpatialGEE package in vignette. 
