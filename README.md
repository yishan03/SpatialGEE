# SpatialGEE

R package for spatial data analysis using generalized estimating equations (GEE). This package enables robust statistical analysis of spatially correlated data and provides tools for hypothesis testing using Wald and generalized score tests.

---

## Reference

Yishan Wang and Peng Wei. **SpatialGEE: A generalized estimating equations framework for spatial data analysis.**  
(*Details of publication pending*).

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

# Load example data
data(example_data)

# Perform a Wald test using the GEE framework
results <- run_gee_wald(data = example_data, compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))
print(results)
```

---

### Feedback and Questions

If you have any questions or encounter any issues while using the package, please feel free to open an issue in this repository.

---

```
