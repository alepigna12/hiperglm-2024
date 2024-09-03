is_default_equal_to_option_model <- function(model_name, n_obs = 32, n_pred = 4, data_seed = 1918) {
  data <- simulate_data(n_obs, n_pred, model_name, seed = data_seed)
  design <- data$design; outcome <- data$outcome
  default_model <- new_regression_model(design, outcome, model_name, SGD_solver = TRUE)
  option_model <- new_regression_model(design, outcome, model_name, SGD_solver = TRUE, via_transp = TRUE)
  return(identical(default_model, option_model))
}

test_that("SGD defaults to via_transp=TRUE option, linear model", {
  expect_true(
    is_default_equal_to_option_model("linear")
  )
})

test_that("SGD defaults to via_transp=TRUE options, logit model", {
  expect_true(
    is_default_equal_to_option_model("logit") 
  )
})