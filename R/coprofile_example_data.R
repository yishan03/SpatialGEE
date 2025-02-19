#' Example Data for Testing Spatial Epigenome-Transcriptome Co-Profiling Integration
#'
#' This dataset contains a subset of processed spatial epigenome-transcriptome mouse brain
#' co-profiling data (Zhang et al., 2023), including 100 genes, intended for testing and 
#' demonstrating the functionality of the package.
#'
#' @format A named list with two elements:
#' \describe{
#'   \item{ATAC}{A `data.frame` with metadata (first 4 columns) and normalized gene accessibility values.}
#'   \item{RNA}{A `data.frame` with metadata (first 4 columns) and gene expression counts.}
#' }
#' Each data frame has the following structure:
#' \describe{
#'   \item{Barcodes}{Unique spot barcodes for spatial transcriptomics (character).}
#'   \item{x}{x-coordinates of spatial locations (numeric).}
#'   \item{y}{y-coordinates of spatial locations (numeric).}
#'   \item{Pathology.Annotations}{Pathology labels (factor with levels `"non-Corpus callosum"` and `"Corpus callosum"`).}
#'   \item{Gabbr2, Pde7b, ..., Itga8}{Accessibility (ATAC) or expression (RNA) values for 100 selected genes (numeric).}
#' }
#'
#' @usage data(coprofile_example_data)
#' @examples
#' data(coprofile_example_data)
"coprofile_example_data"
