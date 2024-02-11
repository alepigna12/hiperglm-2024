calc_loglik <- function(model, reg_coef, ...) {
  UseMethod("calc_loglik")
}

calc_loglink_deriv <- function(model, reg_coef, order, ...) {
  UseMethod("calc_loglink_deriv")
}

calc_loglik.linear_model <- function(model, reg_coef, noise_var = 1) {
  outcome <- model$outcome
  design <- model$design
  predicted_val <- design %*% reg_coef
  loglik <- - 0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_loglink_deriv.linear_model <- function(model, reg_coef, order, noise_var = 1) {
  if (order > 1) {
    stop("2nd+ order derivative calculations are not supported for linear models")
  }
  predicted_val <- model$design %*% reg_coef
  deriv <- (model$outcome - predicted_val) / noise_var
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_loglik.logit_model <- function(model, reg_coef) {
  outcome <- model$outcome
  design <- model$design
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }
  logit_prob <- design %*% reg_coef
  loglik <- sum(n_success * logit_prob - n_trial * log(1 + exp(logit_prob)))
    # TODO: improve numerical stability for logit_prob >> 1
  return(loglik)
}

calc_loglink_deriv.logit_model <- function(model, reg_coef, order) {
  outcome <- model$outcome
  design <- model$design
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }
  logit_prob <- as.vector(design %*% reg_coef)
  predicted_prob <- 1 / (1 + exp(-logit_prob))
  if (order == 1) {
    deriv <- n_success - n_trial * predicted_prob
  } else if (order == 2) {
    deriv <- n_trial * predicted_prob * (1 - predicted_prob)
  } else {
    stop("3rd+ order derivative calculations are not supported for logistic models")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_grad <- function(model, reg_coef, ...) {
  design <- model$design
  loglink_grad <- calc_loglink_deriv(model, reg_coef, order = 1, ...)
  grad <- t(design) %*% loglink_grad
  grad <- as.vector(grad)
  return(grad)
}

calc_hessian <- function(model, reg_coef) {
  design <- model$design
  weight <- calc_loglink_deriv(model, reg_coef, order = 2)
  hess <- - t(design) %*% (outer(weight, rep(1, ncol(design))) * design)
  return(hess)
}

calc_hessian_inverse <- function(model, reg_coef) {
  design <- model$design
  weight <- calc_loglink_deriv(model, reg_coef, order = 2)
  sqrt_weighted_design <- outer(sqrt(weight), rep(1, ncol(design))) * design
  R <- qr_wrapper(sqrt_weighted_design)$R
  inverse <- - invert_gram_mat_from_qr(R)
  return(inverse)
}
