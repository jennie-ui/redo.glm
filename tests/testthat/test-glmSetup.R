context("Testing glmSetup function")

test_that("glmSetup sets up model parameters correctly", {
  # Load necessary libraries
  library(redo.glm)

  # Create example data
  Cell_Differentiation <- read.csv(system.file("data", "Cell_Differentiation.csv", package = "redo.glm"))

  # Test glmSetup function
  glm_full <- glmSetup(res_var = "Count", base_var = "TNF", input_data = Cell_Differentiation, tolerance = 0.00001, log_offset = 0.1, max_iterations = 25)

  # Check the class of the returned object
  expect_true(is(glm_full, "list"))

  # Check the presence of expected elements in the list
  expect_true("histogram" %in% names(glm_full))
  expect_true("transformed_response" %in% names(glm_full))
  expect_true("design_matrix" %in% names(glm_full))
  expect_true("n" %in% names(glm_full))
  expect_true("num_predictors" %in% names(glm_full))

  # Cleanup, detach libraries if necessary
  detach("package:redo.glm", unload = TRUE)
})
