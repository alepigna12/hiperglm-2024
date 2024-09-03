are_default_and_specified_execution_times_close <- function(
    model_name, stepsize, n_obs = 10000, n_pred = 100, data_seed = 1918
) {
  data <- simulate_data(n_obs, n_pred, model_name, seed = data_seed)
  design <- data$design; outcome <- data$outcome
  via_default_sgd <- system.time(hiper_glm(
    design, outcome, model_name, option = list(
    mle_solver = "SGD",
    n_batch = 100L,
    stepsize = stepsize,
    n_epoch = 50,
    subsample_with_replacement = FALSE
  )))
  time_default <- via_default_sgd["user.self"]
  via_specified_sgd <- system.time(hiper_glm(
    design, outcome, model_name, option = list(
      mle_solver = "SGD",
      n_batch = 100L,
      stepsize = stepsize,
      n_epoch = 50,
      subsample_with_replacement = FALSE,
      use_matvec_via_transp = TRUE,
      use_cpp_matvec = TRUE
    )))
  time_specified <- via_specified_sgd["user.self"]
  return(time_default/time_specified > 0.6 & time_default/time_default < 5/3)
}

test_that("SGD defaults to via_transp and use_rcpp options, linear model", {
  expect_true(
    are_default_and_specified_execution_times_close("linear", 0.01)
  )
})

test_that("SGD defaults to via_transp and use_rcpp options, logit model", {
  expect_true(
    are_default_and_specified_execution_times_close("logit", 0.2)  
  )
})