compare_grad_between_options <- function(
    model_name, option_selection_1 = list(), option_selection_2 = list(), n_obs = 32, n_pred = 4, n_test = 10, data_seed = 1918, loc_seed = 615
) {
  data <- simulate_data(n_obs, n_pred, model_name, seed = data_seed)
  design <- data$design; outcome <- data$outcome
  model <- new_regression_model(design, outcome, model_name)
  set.seed(loc_seed)
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    regcoef <- rnorm(n_pred)
    subset_ind <- sample.int(n_obs, 0.5*n_obs)
    first_option_grad <- calc_grad(model, regcoef, subset_ind, option_selection_1$via_transp, option_selection_1$use_rcpp)
    second_option_grad <- calc_grad(model, regcoef, subset_ind, option_selection_2$via_transp, option_selection_2$use_rcpp)
    grads_are_close <- are_all_close(
      first_option_grad, second_option_grad, abs_tol = Inf, rel_tol = 1e-3
    )
  }
  return(grads_are_close)
}

test_that("Gradients with use_matvec_via_transp=TRUE and use_cpp_matvec=TRUE coincide with radients with use_matvec_via_transp=TRUE and use_cpp_matvec=FALSE for linear model", {
  expect_true(
    compare_grad_between_options("linear", list(via_transp = TRUE, use_rcpp = TRUE), list(via_transp = TRUE, use_rcpp = FALSE))
  )
})

test_that("Gradients with use_matvec_via_transp=TRUE and use_cpp_matvec=TRUE coincide with radients with use_matvec_via_transp=FALSE and use_cpp_matvec=TRUE for linear model", {
  expect_true(
    compare_grad_between_options("linear", list(via_transp = TRUE, use_rcpp = TRUE), list(via_transp = FALSE, use_rcpp = TRUE))
  )
})

test_that("Gradients with use_matvec_via_transp=TRUE and use_cpp_matvec=TRUE coincide with radients with use_matvec_via_transp=FALSE and use_cpp_matvec=FALSE for linear model", {
  expect_true(
    compare_grad_between_options("linear", list(via_transp = TRUE, use_rcpp = TRUE), list(via_transp = FALSE, use_rcpp = FALSE))
  )
})

test_that("Gradients with use_matvec_via_transp=TRUE and use_cpp_matvec=TRUE coincide with radients with use_matvec_via_transp=TRUE and use_cpp_matvec=FALSE for logit model", {
  expect_true(
    compare_grad_between_options("logit", list(via_transp = TRUE, use_rcpp = TRUE), list(via_transp = TRUE, use_rcpp = FALSE))
  )
})

test_that("Gradients with use_matvec_via_transp=TRUE and use_cpp_matvec=TRUE coincide with radients with use_matvec_via_transp=FALSE and use_cpp_matvec=TRUE for logit model", {
  expect_true(
    compare_grad_between_options("logit", list(via_transp = TRUE, use_rcpp = TRUE), list(via_transp = FALSE, use_rcpp = TRUE))
  )
})

test_that("Gradients with use_matvec_via_transp=TRUE and use_cpp_matvec=TRUE coincide with radients with use_matvec_via_transp=FALSE and use_cpp_matvec=FALSE for logit model", {
  expect_true(
    compare_grad_between_options("logit", list(via_transp = TRUE, use_rcpp = TRUE), list(via_transp = FALSE, use_rcpp = FALSE))
  )
})