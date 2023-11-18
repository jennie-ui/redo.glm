#' Perform diagnostics and obtain results for a fitted Poisson regression model.
#'
#'This function takes the output from glmSetup() and performs diagnostics for a Poisson regression model.
#'
#' @param myObject Output of glmSetup function
#' @return A list containing various diagnostic results and information about the fitted Poisson regression model.
#'
#' @examples
#' # Example usage:
#' glm_full <- glmSetup(res_var = "Count", base_var = "TNF", input_data = Cell_Differentiation, tolerance = 0.00001, log_offset = 0.1, max_iterations = 25)
#' poisson_diagnostics <- glmPoisson(glm_full)
#'
#' @export
#'
glmPoisson <- function(myObject) {
  # Input validation
  if (!is.list(myObject) || !all(c("response", "transformed_response", "design_matrix", "n", "num_predictors", "tol", "epsilon", "ite_max", "ite") %in% names(myObject))) {
    stop("Invalid myObject. Please provide a valid output from glmSetup.")
  }

  Y <- myObject$response
  log_Y <- myObject$transformed_response
  X <- myObject$design_matrix
  n <- myObject$n
  num_predictors <- myObject$num_predictors

  tol <- myObject$tol
  epsilon <- myObject$epsilon
  ite_max <- myObject$ite_max
  ite <- myObject$ite

  beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y

  while (epsilon > tol & ite <= ite_max){
    eta = X %*% beta
    mu = exp(eta)
    nu = exp(eta)
    V = diag(x = as.vector(nu))
    Z = eta + solve(V) %*% (Y-mu)
    beta_new = solve(t(X) %*% V %*% X) %*% t(X) %*% V %*% Z
    epsilon = sqrt(t(beta_new-beta)%*%(beta_new-beta))
    ite = ite + 1
    beta_t = t(beta)
  }

  Var_beta = solve(t(X) %*% V %*% X)
  rownames(Var_beta)[1] = c("Intercept")
  colnames(Var_beta)[1] = c("Intercept")
  Var_dis<-format(Var_beta, digits = 4, scientific = 4)

  Wald_CI<-data.frame(Estimate=beta, Std_Deviation=sqrt(diag(Var_beta)))
  rownames(Wald_CI)[1] = c("Intercept")
  Wald_CI$CI_Lower = Wald_CI$Estimate-1.96*Wald_CI$Std_Deviation
  Wald_CI$CI_Upper = Wald_CI$Estimate+1.96*Wald_CI$Std_Deviation
  Wald_CI = format(Wald_CI, digits= 7)

  r_p = (Y-mu)/sqrt(mu)

  D_i = 2*(Y*(log(Y/mu)-(Y-mu)))
  r_D = sign(Y-mu)*sqrt(abs(D_i))

  H=sqrt(V) %*% X %*% solve(t(X) %*% V %*% X) %*% t(X) %*% sqrt(V)
  Leverage=diag(H)

  r_PS = r_p/sqrt(1-Leverage)
  r_DS = r_D/sqrt(1-Leverage)

  CookD = myObject$Leverage / (num_predictors * (1 - myObject$Leverage)) * (r_PS)^2

  D = sum(D_i)
  chiP = sum(r_p^2)

  res = data.frame(D = D, chiP = chiP, chi = qchisq(0.95, (n - num_predictors)))

  res = format(res, digits= 3)
  colnames(res) = c("Deviance", "chi^2_P", "chi^2_{n-q(num_predictors),0.95}")

  # Return a more structured result
  return(
    variance_covariance_matrix = Var_dis,
    estimates_and_intervals = Wald_CI,
    residuals = list(pearson = r_p, deviance = r_D),
    leverage = Leverage,
    standardized_residuals = list(pearson = r_PS, deviance = r_DS),
    cooks_distance = CookD,
    goodness_of_fit = res
  )
}
