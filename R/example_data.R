#' Example Data for Testing Functions
#'
#' This dataset contains a subset of processed 10x Genomics breast cancer data, including 3 genes,
#' intended for demonstrating the functionality of the package.
#'
#' @format A data frame with 3168 rows and 7 columns:
#' \describe{
#'   \item{Barcodes}{Unique spot barcodes (character).}
#'   \item{x}{x-coordinates of spatial locations (numeric).}
#'   \item{y}{y-coordinates of spatial locations (numeric).}
#'   \item{Pathology.Annotations}{Pathology labels (factor with levels "Fibrous Tissue" and "Invasive Carcinoma").}
#'   \item{IGHG3}{Expression counts for IGHG3 (integer).}
#'   \item{IGKC}{Expression counts for IGKC (integer).}
#'   \item{IGLC2}{Expression counts for IGLC2 (integer).}
#' }
#' @usage data(example_data)
#' @examples
#' data(example_data)
"example_data"

