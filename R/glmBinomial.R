#' Perform diagnostics and obtain results for a fitted binomial regression model.
#'
#' This function takes the output from glmSetup() and performs diagnostics for a binomial regression model.
#'
#' @param myObject Output of glmSetup function
#'
#' @return A list containing various diagnostic results and information about the fitted binomial regression model.
#'
#' @examples
#' # Example usage:
#' Cell_Differentiation <- read.csv(system.file("data", "Cell_Differentiation.csv", package = "redo.glm"))
#' glm_full <- glmSetup(res_var = "Count", base_var = "TNF", input_data = Cell_Differentiation, tolerance = 0.00001, log_offset = 0.1, max_iterations = 25)
#' glmBinomial(glm_full)
#'
#' @export
#'
glmBinomial <- function(myObject) {
  # Input validation
  if (!is.list(myObject) || !all(c("response", "baseline_variable", "design_matrix", "n", "num_predictors", "tol", "epsilon", "ite_max", "ite") %in% names(myObject))) {
    stop("Invalid myObject. Please provide a valid output from glmSetup.")
  }

  Y <- myObject$response
  baseLine <- myObject$baseline_variable
  Y_star <- Y < baseLine
  X <- myObject$design_matrix
  n <- myObject$n
  num_predictors <- myObject$num_predictors

  tol <- myObject$tol
  epsilon <- myObject$epsilon
  ite_max <- myObject$ite_max
  ite <- myObject$ite

  beta <- rep(0, num_predictors)

  while (epsilon > tol & ite <= ite_max) {
    eta <- X %*% beta
    mu <- exp(eta) / (1 + exp(eta))
    nu <- mu * (1 - mu)
    V <- diag(x = as.vector(nu))
    Z <- eta + solve(V) %*% (Y_star - mu)
    beta_new <- solve(t(X) %*% V %*% X) %*% t(X) %*% V %*% Z
    epsilon = sqrt(sum((beta_new - beta)^2))
    beta <- beta_new
    ite <- ite + 1
  }

  Var_beta <- solve(t(X) %*% V %*% X)
  rownames(Var_beta)[1] <- c("Intercept")
  colnames(Var_beta)[1] <- c("Intercept")
  Var_dis <- format(Var_beta, digits = 4, scientific = 4)

  Std_Deviation <- sqrt(diag(Var_beta))

  Wald_CI<-data.frame(CI_Lower=beta-1.96*Std_Deviation, CI_Upper=beta+1.96*Std_Deviation)
  Wald_CI = format(Wald_CI, digits = 7)

  r_p <- (Y - mu) / sqrt(mu)

  H <- sqrt(V) %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% sqrt(V)
  Leverage <- diag(H)

  r_PS <- r_p / sqrt(1 - Leverage)

  CookD <- Leverage / (num_predictors * (1 - Leverage)) * (r_PS)^2

  chiP <- sum(r_p^2)

  res <- data.frame(
    chiP = chiP,
    chi = qchisq(0.95, (n - num_predictors)) - chiP
  )
  res <- format(res, digits = 3)
  colnames(res) <- c("chi^2_P", "chi^2_{n-q(num_predictors),0.95}")

  newList <- list(
    "variance_covariance_matrix" = Var_dis,
    "estimates" = beta,
    "std" = Std_Deviation,
    "confidence_interval" = Wald_CI,
    "residuals" = list(pearson = r_p),
    "leverage" = Leverage,
    "standardized_residuals" = list(pearson = r_PS),
    "cooks_distance" = CookD,
    "goodness_of_fit" = res
  )

  return(newList)
}
