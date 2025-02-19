# SpatialGEE

R package for spatial data analysis using generalized estimating equations (GEEs). This package enables robust statistical analysis of spatially correlated data and provides tools for hypothesis testing using Wald and generalized score tests. Additionally, SpatialGEE enables spatial co-profiling integration, supporting joint analysis of multi-modal spatial omics data with hypothesis testing and across-data-type FDR control for accurate differential expression in spatial transcriptomics and beyond.

## Reference

Wang, Y., Zang, C., Li, Z., Guo, C. C., Lai, D., & Wei, P. (2025). A comparative study of statistical methods for identifying differentially expressed genes in spatial transcriptomics. Unpublished manuscript.

Wang, Y., & Wei, P. (2025). A Mixture Model Approach for Integrating Spatial Transcriptomics and Epigenomics. Unpublished manuscript.

## Installation

To install the SpatialGEE via GitHub, ensure all required dependencies are installed first. You can use the following commands in R:

```r
required_packages <- c("dplyr", "ggplot2", "geepack", "IMIX", "parallel", "rmarkdown", "knitr")
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
   vignette("SpatialGEE Tutorial for Spatial Coprofiling Integration")
   ```

2. The vignette provides detailed examples for all the functions included in the `SpatialGEE` package. If you have any question, use help().

## Example Workflow

Here’s an example of how to use the package with sample data:

```r
library("SpatialGEE")

data(coprofile_example_data)
atac <- coprofile_example_data$ATAC
rna <- coprofile_example_data$RNA

atac_DE_res <- run_wilcoxon(
  atac, 
  compare_levels = c("non-Corpus callosum", "Corpus callosum"))

rna_DE_res <- run_gee_gst(
  rna, 
  compare_levels = c("non-Corpus callosum", "Corpus callosum"))

merged_pvalue <- merge(atac_DE_res, rna_DE_res, by = "gene") %>%
  rename_with(~ c("atac_pvalue", "rna_pvalue"), starts_with("p_value")) %>%
  tibble::column_to_rownames("gene") %>%
  na.omit()

integration_res <- IMIX::IMIX(data_input = merged_pvalue)
print(integration_res$significant_genes_with_FDRcontrol)
```
