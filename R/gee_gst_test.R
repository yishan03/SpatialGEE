#' Perform GEE analysis with the GST for spatial transcriptomics
#'
#' This function performs Generalized Estimating Equations (GEE) analysis for each gene in the dataset,
#' comparing two user-specified levels of `Pathology.Annotations` using a generalized score test (GST).
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
#' @param k Number of clusters for k-means clustering (default: 100).
#' @param family A description of the error distribution and link function for the model (default: `poisson`).
#' @param corstr A character string specifying the correlation structure (default: `"independence"`).
#' @param cores Number of cores for parallel processing (default: 1).
#' @return A `data.frame` with the following columns:
#'   - `gene`: Gene names.
#'   - `p_value`: P-values from the generalized score test for each gene
#' data(example_data)
#' results <- run_gee_gst(data = example_data, compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))
#' print(results)
#' @export
run_gee_gst <- function(data, compare_levels, k = 100, family = poisson, corstr = "independence", cores = 1) {
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
  data <- data %>%
    filter(Pathology.Annotations %in% compare_levels) %>%
    mutate(Pathology.Annotations = factor(Pathology.Annotations, levels = compare_levels))

  # Add k-means clusters
  if (k < 50) {
    warning("The number of clusters (k) is less than 50. The number of clusters might be too small for asymptotic behavior in GEEs. It is recommended to use approximately 100 clusters.")
  }
  set.seed(123)
  kmeans_result <- kmeans(data[, c("x", "y")], centers = k)
  data$Clusters <- kmeans_result$cluster
  data <- data %>%
    relocate(Clusters, .after = Pathology.Annotations)
  data <- data %>%
    group_by(Clusters) %>%
    filter(n() > 1) %>%  # Remove clusters with only 1 observation
    ungroup() %>%
    arrange(Clusters) %>%
    select(1:5, which(colSums(data[, 6:ncol(data)] != 0) > 0) + 5)

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
  geneexp <- data.frame(lapply(geneexp, as.numeric))
  gene_names <- colnames(geneexp)

  # Collect error log1
  error_log1 <- list()

  # Function to fit null model for a single gene
  fit_null_model  <- function(colname) {
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

      # Fit GEE model
      mgee <- geeglm(
        geneexp ~ 1,
        family = family,
        data = dat,
        id = dat$Clusters,
        corstr = corstr
      )

      return(data.frame(gene = colname, est_inter = coef(summary(mgee))[1,1]))
    }, error = function(e) {
      error_log1[[length(error_log1) + 1]] <<- list(gene = colname, message = e$message)
      return(data.frame(gene = colname, est_inter = NA))
    })
  }

  # Run GEE in parallel or sequentially based on OS
  if (.Platform$OS.type == "windows") {
    null_estimates <- lapply(gene_names, fit_null_model)
  } else {
    null_estimates <- mclapply(gene_names, fit_null_model, mc.cores = cores)
  }

  null_estimates <- do.call(rbind, null_estimates)

  # Process null estimates
  gee_est_inter <- null_estimates %>%
    dplyr::select(gene, est_inter)

  # Prepare beta1 estimates as zero (under null hypothesis)
  gee_est_beta1 <- gee_est_inter %>%
    dplyr::select(gene) %>%
    mutate(est_beta1 = 0)

  # Convert null estimates into named lists for easy lookup
  est_inter <- setNames(as.list(gee_est_inter$est_inter), gee_est_inter$gene)
  est_beta1 <- setNames(as.list(gee_est_beta1$est_beta1), gee_est_beta1$gene)

  # Generalized Score Test (GST) function
  gst <- function(gene_expression_column, interEst, beta1Est, data) {

    metadata <- data[, 1:5]
    tempdata <- cbind(metadata, gene_expression_column)

    m <- unique(data$Clusters)

    a_beta <- matrix(0, 2, 2)
    b_beta <- matrix(0, 2, 2)
    score <- matrix(0, 2, 1)

    for(i in m) {

      X <- model.matrix(~Pathology.Annotations, data = subset(data, Clusters == i))

      beta <- t(as.matrix(c(as.numeric(interEst), as.numeric(beta1Est))))

      mu <- exp(beta %*% t(X))

      D <- sweep(X, MARGIN = 1, mu, `*`)

      c_beta <- diag(as.numeric(sqrt(mu)))
      r_alpha <- diag(dim(c_beta)[1])
      V <- c_beta %*% r_alpha %*% c_beta

      a_beta_i <- t(D) %*% solve(V) %*% D
      a_beta <- a_beta + a_beta_i

      r <- t(as.matrix(subset(tempdata, Clusters == i)$gene_expression_column - mu))

      b_beta_i <- t(D) %*% solve(V) %*% r %*% t(r) %*% solve(V) %*% D
      b_beta <- b_beta + b_beta_i

      score_i <- t(D) %*% solve(V) %*% r
      score <- score + score_i
    }

    var_beta <- solve(a_beta) %*% b_beta %*% solve(t(a_beta))

    t_gs <- t(score) %*% solve(t(a_beta)) %*% solve(var_beta) %*% solve(a_beta) %*% score

    p_value <- pchisq(t_gs, 1, lower.tail = FALSE) # double check

    return(p_value)
  }

  # Run GST for each gene
  gene_results <- list()
  # Run GST for each gene based on the operating system
  if (.Platform$OS.type == "windows") {
    result <- lapply(names(geneexp), function(colname) {
      tryCatch({
        p_value <- gst(geneexp[[colname]], est_inter[[colname]], est_beta1[[colname]], data)
        gene_results[[colname]] <- data.frame(gene = colname, p_value = p_value)
      }, error = function(e) {
        gene_results[[colname]] <- data.frame(gene = colname, p_value = NA)  # Return NA if an error occurs
      })
    })
  } else {
    result <- mclapply(names(geneexp), function(colname) {
      tryCatch({
        p_value <- gst(geneexp[[colname]], est_inter[[colname]], est_beta1[[colname]], data)
        gene_results[[colname]] <- data.frame(gene = colname, p_value = p_value)
      }, error = function(e) {
        gene_results[[colname]] <- data.frame(gene = colname, p_value = NA)  # Return NA if an error occurs
      })
    }, mc.cores = cores)
  }

  final_results <- do.call(rbind, result)

  return(final_results)
}
