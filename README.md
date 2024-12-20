# SpatialGEE

R package for spatial data analysis using generalized estimating equations (GEEs). This package enables robust statistical analysis of spatially correlated data and provides tools for hypothesis testing using Wald and generalized score tests.

## Reference

Wang, Y., Zang, C., Li, Z., Guo, C. C., Lai, D., & Wei, P. (2024). A comparative study of statistical methods for identifying differentially expressed genes in spatial transcriptomics. Unpublished manuscript.

## Installation

To install the SpatialGEE via GitHub, ensure all required dependencies are installed first. You can use the following commands in R:

```r
required_packages <- c("dplyr", "ggplot2", "geepack", "parallel", "rmarkdown", "knitr")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("yishan03/SpatialGEE")
```

## To Get Started with the Package

1. Load the package and explore its vignette:

   ```r
   library("SpatialGEE")
   vignette("SpatialGEE")
   ```

2. The vignette provides detailed examples for all the functions included in the `SpatialGEE` package. If you have any question, use help().

## Example Workflow

Here’s an example of how to use the package with sample data:

```r
library("SpatialGEE")

data(example_data)

results <- run_gee_gst(data = example_data, compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))
print(results)
```
