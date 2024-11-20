#' Visualize spatial transcriptomics data with spatial clusters
#'
#' This function visualizes the spatial coordinates of the transcriptomics data,
#' colored by the k-means clusters.
#'
#' @param data A `data.frame` containing metadata in the first 4 columns and optionally gene expression data.
#' The metadata must include the following columns with these exact names:
#'   - `Barcodes`: Spot barcodes for spatial transcriptomics.
#'   - `x`: Spatial x-coordinates of the spots.
#'   - `y`: Spatial y-coordinates of the spots.
#'   - `Pathology.Annotations`: Pathology annotations for each spot.
#' @param k Number of clusters for k-means clustering (default: 100).
#' @return A ggplot object showing the spatial distribution of clusters.
#' @examples
#' data(example_data)
#' plot <- plot_spatial_clusters(example_data)
#' print(plot)
#' @export
plot_spatial_clusters <- function(data, k = 100) {
  library(dplyr)
  library(ggplot2)

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

  # Add k-means clusters
  set.seed(123)
  kmeans_result <- kmeans(data[, c("x", "y")], centers = k)
  data$Clusters <- as.factor(kmeans_result$cluster) # Convert clusters to a factor
  data <- data %>% relocate(Clusters, .after = Pathology.Annotations) # Move Clusters column

  # Subset relevant columns for clustering and visualization
  xyz <- data %>%
    select(x, y, Clusters)

  # Function to calculate convex hull for each cluster
  find_hull <- function(data) {
    data %>%
      group_by(Clusters) %>%
      slice(chull(x, y))
  }

  # Get convex hull points for cluster boundaries
  hull_points <- find_hull(xyz)

  # Generate the plot
  plot <- ggplot() +
    geom_raster(data = xyz, aes(x = x, y = y, fill = factor(Clusters))) +  # Raster plot of clusters
    geom_polygon(
      data = hull_points, aes(x = x, y = y, group = Clusters),
      color = "black", fill = NA, linewidth = 0.8
    ) +  # Cluster boundaries
    labs(
      title = "Tissue with Spatial Clusters",
      fill = "Cluster ID"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )

  return(plot)
}
