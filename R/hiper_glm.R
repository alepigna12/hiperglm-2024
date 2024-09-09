#' @export
hiper_glm <- function(design, outcome, model_name = "linear", option = list()) {
  supported_model <- c("linear", "logit")
  if (!(model_name %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model_name))
  }
  option <- initialize_option(option)
  model <- new_regression_model(design, outcome, model_name, option$use_matvec_via_transp)
  hglm_out <- find_mle(model, option)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}

find_mle <- function(model, option) {
  if (option$mle_solver == "newton") {
    if (model$name == 'linear') {
      result <- solve_via_least_sq(model$design, model$outcome)
    } else {
      result <- solve_via_newton(
        model, option$n_max_iter, option$rel_tol, option$abs_tol
      )
    }
  } else if (option$mle_solver %in% c("BFGS", "CG", "L-BFGS-B")) {
    result <- solve_via_optim(model, option$mle_solver)
  } else if (option$mle_solver == "SGD") {
    result <- solve_via_SGD(model, option$n_batch, option$stepsize, option$n_epoch, option$subsample_with_replacement, option$use_matvec_via_transp, option$use_cpp_matvec)
  } else {
    stop("Unsupported MLE solver type.")
  }
  return(result)
}

solve_via_least_sq <- function(design, outcome) {
  ls_result <- solve_least_sq_via_qr(design, outcome)
  mle_coef <- ls_result$solution
  noise_var <- mean((outcome - design %*% mle_coef)^2)
  n_obs <- nrow(design); n_pred <- ncol(design)
  noise_var <- noise_var / (1 - n_pred / n_obs) 
    # Use the same nearly-unbiased estimator as in `stats::lm`
  cov_est <- noise_var * invert_gram_mat_from_qr(ls_result$R)
  return(list(coef = mle_coef, cov_est = cov_est))
}

solve_via_newton <- function(model, n_max_iter, rel_tol, abs_tol) {
  if (is.null(n_max_iter)) { n_max_iter <- 25L }
  if (is.null(rel_tol)) { rel_tol <- 1e-6 }
  if (is.null(abs_tol)) { abs_tol <- 1e-6 }
  coef_est <- rep(0, ncol(model$design))
  n_iter <- 0L
  max_iter_reached <- FALSE
  converged <- FALSE
  curr_loglik <- calc_loglik(model, coef_est)
  while (!(converged || max_iter_reached)) {
    prev_loglik <- curr_loglik
    coef_est <- take_one_newton_step(model, coef_est)
    curr_loglik <- calc_loglik(model, coef_est)
    converged <- (
      2 * abs(curr_loglik - prev_loglik) < (abs_tol + rel_tol * abs(curr_loglik))
    )
    n_iter <- n_iter + 1L
    max_iter_reached <- (n_iter == n_max_iter)
  }
  if (max_iter_reached && !converged) {
    warning("Newton's method did not converge. The estimates may be meaningless.")
  }
  cov_est <- - calc_hessian_inverse(model, coef_est)
  return(list(
    coef = coef_est, cov = cov_est, 
    converged = converged, n_iter = n_iter
  ))
}

take_one_newton_step <- function(model, coef_est, solver = "weighted-least-sq") {
  if (solver == "weighted-least-sq") {
    loglink_grad <- calc_loglink_deriv(model, coef_est, order = 1)
    weight <- calc_loglink_deriv(model, coef_est, order = 2)
    if (any(weight == 0)) {
      stop("Exact 0 or 1 found in predicted probability while solving for MLE.")
        # TODO: pursue alternative path forward in this case. Maybe just fall 
        # back on a Newton step with explicit computation of weighted Hessian. 
    }
    ls_target_vec <- loglink_grad / weight
    coef_update <- solve_least_sq_via_qr(model$design, ls_target_vec, weight)$solution
  } else if (solver == "normal-eq") {
    grad <- calc_grad(model, coef_est)
    hess <- calc_hessian(model, coef_est)
    coef_update <- - solve(hess, grad)
  } else {
    stop("Unsupported solver type.")
  }
  coef_est <- coef_est + coef_update
  return(coef_est)
}

solve_via_optim <- function(model, method) {
  init_coef <- rep(0, ncol(model$design))
  obj_fn <- function (coef) {
    calc_loglik(model, coef) 
  }
  obj_grad <- function (coef) {
    calc_grad(model, coef)
  }
  optim_result <- stats::optim(
    init_coef, obj_fn, obj_grad, method = method,
    control = list(fnscale = -1) # Maximize the function
  )
  optim_converged <- (optim_result$convergence == 0L)
  if (!optim_converged) {
    warning("Optimization did not converge. The estimates may be meaningless.")
  }
  return(list(coef = optim_result$par))
}

solve_via_SGD <- function(model, n_batch, stepsize, n_epoch, replacement=FALSE, via_transp, use_rcpp) {
  n_obs <- model$n_obs
  n_pred <- model$n_pred
  coef_est <- rep(0, n_pred)
  for (epoch in 1:n_epoch) {
    rand_indices <- sample.int(n_obs, size = n_obs, replace = replacement)
    batch_cum_size <- floor(seq(0, n_obs, length.out = n_batch + 1L))
    subset_ind_per_batch <- lapply(
      1:n_batch,
      function(i) rand_indices[(batch_cum_size[i] + 1):batch_cum_size[i + 1]]
    )
    for(batch in 1:n_batch){
      subset_index <- subset_ind_per_batch[[batch]]
      n_sub <- length(subset_index)
      coef_est <- coef_est + stepsize / n_sub * calc_grad(model, coef_est, subset_index, via_transp, use_rcpp)
    }
  }
  return(list(coef = coef_est))
}
