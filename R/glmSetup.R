#' Set up for IRWLS generalized linear model
#'
#' This function prepares the data and parameters for fitting a generalized
#' linear model using Iteratively Reweighted Least Squares (IRWLS). Note that
#' this model setup currently does not handle interaction terms.
#'
#' @param res_var Response variable name
#' @param base_var Baseline variable name
#' @param input_data Input data frame
#' @param tolerance Tolerance for convergence
#' @param log_offset Offset for log transformation
#' @param max_iterations Maximum number of iterations
#' @return List containing setup information
#'
#' @seealso
#' Use glmPoisson() or glmBinomial() to perform model diagnostics and obtain detailed results.
#'
#' @export
#' @examples
#' # Example usage:
#' Cell_Differentiation <- read.csv(system.file("data", "Cell_Differentiation.csv", package = "redo.glm"))
#' glm_full <- glmSetup(res_var = "Count", base_var = "TNF", input_data = Cell_Differentiation, tolerance = 0.00001, log_offset = 0.1, max_iterations = 25)
#'
glmSetup <- function(res_var, base_var, input_data, tolerance, log_offset, max_iterations) {
  # Input validation
  if (!is.character(res_var) || !is.character(base_var)) {
    stop("res_var and base_var should be character strings.")
  }

  if (!is.data.frame(input_data)) {
    stop("input_data should be a data frame.")
  }

  if (!(res_var %in% colnames(input_data) && base_var %in% colnames(input_data))) {
    stop("res_var and base_var should exist in input_data.")
  }

  # Extract relevant columns
  df <- data.frame(input_data)
  n <- nrow(df)
  baseline_position <- which(colnames(input_data) == base_var)
  y_position <- which(colnames(input_data) == res_var)
  baseLine <- input_data[, baseline_position]
  Y <- input_data[, y_position]
  log_Y <- log(Y + log_offset)
  X <- as.matrix(cbind(1, df[-c(y_position)]))
  num_predictors <- ncol(X)

  # Set up convergence parameters
  tol <- tolerance
  epsilon <- 99
  ite_max <- max_iterations
  ite <- 0

  # Plotting
  distribution_plot <- hist(Y,
                            main = paste0("Histogram of ", colnames(df)[1]),
                            xlab = colnames(df)[1])

  # Create and return a list
  my_object <- list(
    "histogram" = distribution_plot,
    "response" = Y,
    "transformed_response" = log_Y,
    "baseline_variable" = baseLine,
    "design_matrix" = X,
    "n" = n,
    "num_predictors" = num_predictors,
    "tol" = tol,
    "epsilon" = epsilon,
    "ite_max" = ite_max,
    "ite" = ite
  )

  return(my_object)
}
