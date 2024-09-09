test_that("SGD defaults to option$use_matvec_via_transp=TRUE and option$use_cpp_matvec=TRUE", {
  option_default = list(
    mle_solver = "SGD",
    n_batch = 100L,
    stepsize = 0.01,
    n_epoch = 10,
    subsample_with_replacement = FALSE
  )
  option_default = initialize_option(option_default)
  
  expect_true(
    option_default$use_matvec_via_transp && option_default$use_cpp_matvec
  )
})
