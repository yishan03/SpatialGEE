# SpatialGEE

R package for spatial data analysis using generalized estimating equations (GEEs). This package enables robust statistical analysis of spatially correlated data and provides tools for hypothesis testing using Wald and generalized score tests. Additionally, SpatialGEE enables spatial co-profiling integration, supporting joint analysis of multi-modal spatial omics data with hypothesis testing and across-data-type FDR control for accurate differential expression in spatial transcriptomics and beyond.

## References

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
   vignette("SpatialGEE") # See Wang et al. (2025) for details on differential expression analysis for a single spatial data type
   vignette("SpatialGEEIntegration") # See Wang and Wei (2025) for details on integrating multiple spatial data types
   ```

2. The vignettes provide detailed examples for all the functions and pipeline included in the SpatialGEE package. If you have any question, use help().

## Example Workflow

1. Differential expression analysis for a single spatial data type (Wang et al., 2025)

This example demonstrates how to apply the generalized score test (GST) to identify differentially expressed genes in spatial transcriptomics data.

   ```r
   library("SpatialGEE")

   data(example_data)

   results_gst <- run_gee_gst(
     data = example_data, 
     compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))

   print(results_gst)
   ```

2. Integration of Multiple Spatial Data Types (Wang and Wei, 2025)

This example demonstrates how to integrate spatial transcriptomics and epigenomics data.

   ```r
   library("SpatialGEE")

   data(coprofile_example_data)
   atac <- coprofile_example_data$ATAC
   rna <- coprofile_example_data$RNA

   atac_DE_res <- run_wilcoxon(
     data = atac, 
     compare_levels = c("non-Corpus callosum", "Corpus callosum"))

   rna_DE_res <- run_gee_gst(
     data = rna, 
     compare_levels = c("non-Corpus callosum", "Corpus callosum"))

   merged_pvalue <- merge(atac_DE_res, rna_DE_res, by = "gene") %>%
     rename_with(~ c("atac_pvalue", "rna_pvalue"), starts_with("p_value")) %>%
     tibble::column_to_rownames("gene") %>%
     na.omit()

   integration_res <- IMIX::IMIX(data_input = merged_pvalue)
   print(integration_res$significant_genes_with_FDRcontrol)
   ```
