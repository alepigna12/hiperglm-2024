compare_analytical_and_numerical_grad <- function(
  model_name, n_obs = 32, n_pred = 4, n_test = 10, data_seed = 1918, loc_seed = 615
) {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model_name, seed = data_seed)
  design <- data$design; outcome <- data$outcome
  model <- new_regression_model(design, outcome, model_name)
  loglik_func <- function (coef) { calc_loglik(model, coef) }
  set.seed(loc_seed)
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    regcoef <- rnorm(n_pred)
    analytical_grad <- calc_grad(model, regcoef)
    numerical_grad <- approx_grad_via_finite_diff(loglik_func, regcoef)
    grads_are_close <- are_all_close(
      analytical_grad, numerical_grad, abs_tol = Inf, rel_tol = 1e-3
    )
  }
  return(grads_are_close)
}

test_that("linear model's analytical gradient is close to numerical one", {
  expect_true(
    compare_analytical_and_numerical_grad("linear")
  )
})

test_that("logit model's analytical gradient is close to numerical one", {
  expect_true(
    compare_analytical_and_numerical_grad("logit")
  )
})

test_that("subset grad of full data coincide with full grad of subset data", {
  
  n_obs <- 16; subset_size <- 5; n_pred <- 4
  set.seed(615)
  coef <- rnorm(n_pred)
  subset_ind <- sample.int(n_obs, size = subset_size)
  
  for (model_name in c("linear", "logit")) {
    data <- simulate_data(n_obs, n_pred, model_name, seed = 1918)
    design <- data$design; outcome <- data$outcome
    
    model <- new_regression_model(design, outcome, model_name)
    sub_model <- new_regression_model(
      design[subset_ind, ], outcome[subset_ind], model_name
    )
    expect_true(are_all_close(
      calc_grad(model, coef, subset_ind), 
      calc_grad(sub_model, coef, NULL)
    ))
  }
})
