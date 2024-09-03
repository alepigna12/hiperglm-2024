calc_loglik <- function(model, reg_coef, ...) {
  UseMethod("calc_loglik")
}

calc_loglink_deriv <- function(model, reg_coef, order, subset_ind = NULL, via_transp = TRUE, use_rcpp = TRUE, ...) {
  UseMethod("calc_loglink_deriv")
}

calc_loglik.linear_model <- function(model, reg_coef, noise_var = 1) {
  outcome <- model$outcome
  predicted_val <- matvec_by_design(model, reg_coef)
  loglik <- - 0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_loglink_deriv.linear_model <- function(model, reg_coef, order, noise_var = 1, subset_ind = NULL, via_transp = TRUE, use_rcpp = TRUE) {
  outcome <- get_outcome(model, subset_ind)
  if (order > 1) {
    stop("2nd+ order derivative calculations are not supported for linear models")
  }
  predicted_val <- matvec_by_design(model, reg_coef, subset_ind, via_transp, use_rcpp)
  deriv <- (outcome - predicted_val) / noise_var
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_loglik.logit_model <- function(model, reg_coef) {
  outcome_pair <- get_logit_outcome_pair(model)
  n_success <- outcome_pair$n_success
  n_trial <- outcome_pair$n_trial
  logit_prob <- matvec_by_design(model, reg_coef)
  loglik <- sum(n_success * logit_prob - n_trial * log(1 + exp(logit_prob)))
    # TODO: improve numerical stability for logit_prob >> 1
  return(loglik)
}

calc_loglink_deriv.logit_model <- function(model, reg_coef, order, subset_ind = NULL, via_transp = TRUE, use_rcpp = TRUE) {
  outcome_pair <- get_logit_outcome_pair(model, subset_ind)
  n_success <- outcome_pair$n_success
  n_trial <- outcome_pair$n_trial
  logit_prob <- matvec_by_design(model, reg_coef, subset_ind, via_transp, use_rcpp)
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

calc_grad <- function(model, reg_coef, subset_ind = NULL, via_transp = TRUE, use_rcpp = TRUE, ...) {
  if (is.null(via_transp)) {
    via_transp = TRUE
  }
  if (is.null(use_rcpp)) {
    use_rcpp = TRUE
  }
  loglink_grad <- calc_loglink_deriv(model, reg_coef, order = 1, subset_ind = subset_ind, via_transp = via_transp, use_rcpp = use_rcpp, ...)
  grad <- matvec_by_design_transp(model, loglink_grad, subset_ind = subset_ind, via_transp = via_transp, use_rcpp = use_rcpp)
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
