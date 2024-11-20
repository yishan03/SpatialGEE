test_that("run_wilcoxon works with sample data", {
  data(example_data) # Load the sample data
  results <- run_wilcoxon(data = sample_data, compare_levels = c("Fibrous Tissue", "Invasive Carcinoma"))

  # Perform checks
  expect_s3_class(results, "data.frame")
  expect_named(results, c("gene", "p_value"))
  expect_true(all(results$p_value >= 0 & results$p_value <= 1, na.rm = TRUE))
  expect_equal(nrow(results), 3) # Check if all genes were tested
})
