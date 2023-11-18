# testthat tests for glmPoisson

context("Testing glmPoisson function")

test_that("glmPoisson estimates Poisson regression model correctly", {
  # Load necessary libraries
  library(redo.glm)

  # Create example data
  Cell_Differentiation <- read.csv(system.file("data", "Cell_Differentiation.csv", package = "redo.glm"))

  # Set up model parameters
  glm_full <- glmSetup(res_var = "Count", base_var = "TNF", input_data = Cell_Differentiation, tolerance = 0.00001, log_offset = 0.1, max_iterations = 25)

  # Fit Poisson model
  glm_poisson <- glmPoisson(glm_full)

  # Check the class of the returned object
  expect_true(is(glm_poisson, "list"))

  # Check the presence of expected elements in the list
  expect_true("variance_covariance_matrix" %in% names(glm_poisson))
  expect_true("estimates" %in% names(glm_poisson))
  expect_true("std" %in% names(glm_poisson))
  expect_true("confidence_interval" %in% names(glm_poisson))
  expect_true("residuals" %in% names(glm_poisson))
  expect_true("leverage" %in% names(glm_poisson))
  expect_true("standardized_residuals" %in% names(glm_poisson))
  expect_true("cooks_distance" %in% names(glm_poisson))
  expect_true("goodness_of_fit" %in% names(glm_poisson))

  # Cleanup, detach libraries if necessary
  detach("package:redo.glm", unload = TRUE)
})
