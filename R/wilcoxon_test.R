#' Perform Wilcoxon rank-sum test for spatial transcriptomics
#'
#' This function performs the Wilcoxon rank-sum test for each gene in the dataset,
#' comparing two user-specified levels of `Pathology.Annotations`.
#'
#' @param data A `data.frame` containing metadata in the first 4 columns and gene expression data in the remaining columns.
#' The metadata must include the following columns with these exact names:
#'   - `Barcodes`: Spot barcodes for spatial transcriptomics.
#'   - `x`: Spatial x-coordinates of the spots.
#'   - `y`: Spatial y-coordinates of the spots.
#'   - `Pathology.Annotations`: Pathology annotations for each spot.
#' The remaining columns should contain gene expression data.
#' @param compare_levels A character vector of length 2 specifying the levels of
#' `Pathology.Annotations` to compare.
#' @param cores Number of cores for parallel processing (default: 1).
#' @return A `data.frame` with the following columns:
#'   - `gene`: Gene names.
#'   - `p_value`: P-values from the Wilcoxon rank-sum test for each gene.
#' @examples
#' data(example_data)
#' results <- run_wilcoxon(data = example_data, compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))
#' print(results)
#' @export
run_wilcoxon <- function(data, compare_levels, cores = 1) {
  library(dplyr)
  library(parallel)

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
  data <- data %>% filter(Pathology.Annotations %in% compare_levels)

  # Automatically detect gene expression columns
  gene_expression_columns <- 5:ncol(data)
  if (length(gene_expression_columns) == 0) {
    stop("No gene expression columns detected. Ensure that gene expression data follows the metadata.")
  }
  if (!all(sapply(data[, gene_expression_columns], is.numeric))) {
    stop("Non-numeric data detected in gene expression columns. Ensure that only numeric data is present.")
  }

  # Extract metadata and gene expression
  meta <- data[, required_columns]
  geneexp <- data[, gene_expression_columns, drop = FALSE]
  gene_names <- colnames(geneexp)

  # Collect error log
  error_log <- list()

  # Function to perform Wilcoxon test for a single gene
  wilcoxon_fn <- function(colname) {
    tryCatch({
      dat <- data.frame(geneexp = geneexp[[colname]], meta)
      group1 <- dat %>% filter(Pathology.Annotations == compare_levels[1]) %>% pull(geneexp)
      group2 <- dat %>% filter(Pathology.Annotations == compare_levels[2]) %>% pull(geneexp)

      # Handle missing data
      group1 <- group1[!is.na(group1)]
      group2 <- group2[!is.na(group2)]

      # Validate that both groups have non-zero lengths
      if (length(group1) == 0 || length(group2) == 0) {
        stop("One or both groups have no data. Ensure valid compare_levels and data.")
      }

      # Perform Wilcoxon test
      p_value <- wilcox.test(group1, group2)$p.value
      return(data.frame(gene = colname, p_value = p_value))
    }, error = function(e) {
      error_log[[length(error_log) + 1]] <<- list(gene = colname, message = e$message)
      return(data.frame(gene = colname, p_value = NA))
    })
  }

  # Run Wilcoxon test in parallel or sequentially based on OS
  if (.Platform$OS.type == "windows") {
    results <- lapply(gene_names, wilcoxon_fn)
  } else {
    results <- mclapply(gene_names, wilcoxon_fn, mc.cores = cores)
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
