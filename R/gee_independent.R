#' Perform GEE analysis by treating each spatial location as an independent observation for spatial transcriptomics
#'
#' This function performs Generalized Estimating Equations (GEE) analysis for each gene in the dataset,
#' comparing two user-specified levels of `Pathology.Annotations` while treating each spatial location as independent.
#' The correlation structure is set to "independence," meaning no correlation is assumed between observations.
#'
#' @param data A `data.frame` containing metadata in the first 4 columns and gene expression data in the remaining columns.
#' The metadata must include the following columns with these exact names:
#'   - `Barcodes`: Spot barcodes for spatial transcriptomics.
#'   - `x`: Spatial x-coordinates of the spots.
#'   - `y`: Spatial y-coordinates of the spots.
#'   - `Pathology.Annotations`: Pathology annotations for each spot.
#' The remaining columns should contain gene expression data.
#' @param compare_levels A character vector of length 2 specifying the levels of
#' `Pathology.Annotations` to compare. The first element of `compare_levels`
#' will be set as the reference level, and the second element will be set as
#' the comparison level in the factor.
#' @param family A description of the error distribution and link function for the model (default: `poisson`).
#' @param cores Number of cores for parallel processing (default: 1).
#' @return A `data.frame` with the following columns:
#'   - `gene`: Gene names.
#'   - `p_value`: P-values from the independent GEE analysis for each gene.
#' @examples
#' data(example_data)
#' results <- run_indep_gee(data = example_data, compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))
#' print(results)
#' @export
run_indep_gee <- function(data, compare_levels, family = poisson, cores = 1) {
  library(dplyr)
  library(parallel)
  library(geepack)

  # Validate required metadata columns
  required_columns <- c("Barcodes", "x", "y", "Pathology.Annotations")
  first_four_columns <- colnames(data)[1:4]
  missing_columns <- setdiff(required_columns, first_four_columns)
  if (length(missing_columns) > 0) {
    stop(paste(
      "The following required metadata columns are missing from the first 4 columns:",
      paste(missing_columns, collapse = ", ")
    ))
  }

  # Validate input levels
  available_levels <- unique(data$Pathology.Annotations)
  if (length(compare_levels) != 2) {
    stop("Exactly two levels of Pathology.Annotations must be specified in `compare_levels`.")
  }

  if (!all(compare_levels %in% available_levels)) {
    stop(paste(
      "One or more specified levels are not present in Pathology.Annotations.",
      "Available levels are:", paste(available_levels, collapse = ", ")
    ))
  }

  # Filter data to include only the specified levels
  data <- data %>% filter(Pathology.Annotations %in% compare_levels) %>%
    mutate(Pathology.Annotations = factor(Pathology.Annotations, levels = compare_levels),
           Clusters = row_number()) %>%
    relocate(Clusters, .after = Pathology.Annotations) %>%
    arrange(Clusters)

  # Automatically detect gene expression columns
  gene_expression_columns <- 6:ncol(data)
  if (length(gene_expression_columns) == 0) {
    stop("No gene expression columns detected. Ensure that gene expression data follows the metadata.")
  }
  if (!all(sapply(data[, gene_expression_columns], is.numeric))) {
    stop("Non-numeric data detected in gene expression columns. Ensure that only numeric data is present.")
  }

  # Extract metadata and gene expression
  meta <- data %>% select(all_of(required_columns), Clusters)
  geneexp <- data[, gene_expression_columns, drop = FALSE]
  gene_names <- colnames(geneexp)

  # Collect error log
  error_log <- list()

  # Function to perform independent GEE analysis for a single gene
  gee_fn <- function(colname) {
    tryCatch({
      dat <- data.frame(meta, geneexp = geneexp[[colname]])

      # Extract groups based on Pathology.Annotations
      level1 <- dat %>% filter(Pathology.Annotations == compare_levels[1]) %>% pull(geneexp)
      level2 <- dat %>% filter(Pathology.Annotations == compare_levels[2]) %>% pull(geneexp)

      # Handle missing data
      level1 <- level1[!is.na(level1)]
      level2 <- level2[!is.na(level2)]

      # Validate that both groups have non-zero lengths
      if (length(level1) == 0 || length(level2) == 0) {
        stop("One or both levels have no data. Ensure valid compare_levels and data.")
      }

      # Fit GEE model with independent correlation structure
      mgee <- geeglm(
        geneexp ~ Pathology.Annotations,
        family = family,
        data = dat,
        id = dat$Clusters,
        corstr = "independence"
      )

      # Extract p-value for Pathology.Annotations
      p_value <- coef(summary(mgee))[2, 4]

      return(data.frame(gene = colname, p_value = p_value))
    }, error = function(e) {
      error_log[[length(error_log) + 1]] <<- list(gene = colname, message = e$message)
      return(data.frame(gene = colname, p_value = NA))
    })
  }

  # Run GEE in parallel or sequentially based on OS
  if (.Platform$OS.type == "windows") {
    results <- lapply(gene_names, gee_fn)
  } else {
    results <- mclapply(gene_names, gee_fn, mc.cores = cores)
  }

  results <- do.call(rbind, results)

  # Report errors if any
  if (length(error_log) > 0) {
    cat("\nThe following errors were encountered:\n\n")
    for (i in seq_along(error_log)) {
      cat(paste0("Error [", i, "]:\n"))
      cat(paste0("  Gene: ", error_log[[i]]$gene, "\n"))
      cat(paste0("  Message: ", error_log[[i]]$message, "\n\n"))
    }
  }

  return(results)
}
