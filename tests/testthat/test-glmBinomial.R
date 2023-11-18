context("Testing glmBinomial function")

test_that("glmBinomial estimates and evaluates binomial regression model correctly", {
  # Load necessary libraries
  library(redo.glm)
  library(testthat)

  # Create example data
  Cell_Differentiation <- read.csv(system.file("data", "Cell_Differentiation.csv", package = "redo.glm"))

  # Setup model parameters using glmSetup (assuming you have already tested glmSetup)
  glm_full <- glmSetup(res_var = "Count", base_var = "TNF", input_data = Cell_Differentiation, tolerance = 0.00001, log_offset = 0.1, max_iterations = 25)

  # Test glmBinomial function
  glm_binomial <- glmBinomial(glm_full)

  # Check the class of the returned object
  expect_true(is(glm_binomial, "list"))

  # Check the presence of expected elements in the list
  expect_true("variance_covariance_matrix" %in% names(glm_binomial))
  expect_true("estimates" %in% names(glm_binomial))
  expect_true("std" %in% names(glm_binomial))
  expect_true("confidence_interval" %in% names(glm_binomial))
  expect_true("residuals" %in% names(glm_binomial))
  expect_true("leverage" %in% names(glm_binomial))
  expect_true("standardized_residuals" %in% names(glm_binomial))
  expect_true("cooks_distance" %in% names(glm_binomial))
  expect_true("goodness_of_fit" %in% names(glm_binomial))

  # Cleanup, detach libraries if necessary
  detach("package:redo.glm", unload = TRUE)
})
