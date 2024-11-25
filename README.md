# SpatialGEE

R package for spatial data analysis using generalized estimating equations (GEEs). This package enables robust statistical analysis of spatially correlated data and provides tools for hypothesis testing using Wald and generalized score tests.

---

## Reference

Wang, Y., Zang, C., Li, Z., Guo, C. C., Lai, D., & Wei, P. (2024). A comparative study of statistical methods for identifying differentially expressed genes in spatial transcriptomics. Unpublished manuscript.

---

## Installation

Install the **SpatialGEE** package via GitHub:

```r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("yishan03/SpatialGEE", build_vignettes = TRUE)
```

---

## To Get Started with the Package

1. Load the package and explore its vignette:

   ```r
   library("SpatialGEE")
   vignette("SpatialGEE")
   ```

2. The vignette provides detailed examples for all the functions included in the `SpatialGEE` package. You can also access the function-specific documentation using `help()`.

---

## Example Workflow

Here’s an example of how to use the package with sample data:

```r
library("SpatialGEE")

data(example_data)

results <- gee_gst_test(data = example_data, compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))
print(results)
```