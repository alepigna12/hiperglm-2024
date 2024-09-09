test_that("Rcpp implementation of row subset matvec coincides with base R result", {
  set.seed(111)
  
  n <- 10^4; m <- 10; subset_frac = 0.1
  A <- matrix(rnorm(n * m), n, m)
  b <- rnorm(m)
  row_index <- sample.int(n, subset_frac * n)
  
  expect_true(are_all_close(
    A[row_index, ] %*% b,
    row_subset_matvec(A, b, row_index)
  ))
})


test_that("Effficient Rcpp implementation of row subset matvec via transpose coincides with base R result", {
  set.seed(111)
  
  n <- 10^4; m <- 10; subset_frac = 0.1
  A <- matrix(rnorm(n * m), n, m)
  b <- rnorm(m)
  row_index <- sample.int(n, subset_frac * n)
  
  expect_true(are_all_close(
    A[row_index, ] %*% b,
    row_subset_matvec_via_transpose(t(A), b, row_index)
  ))
})


test_that("Effficient Rcpp implementation of col subset matvec coincides with base R result", {
  set.seed(111)
  
  n <- 10^2; m <- 10^3; subset_frac = 0.1
  A <- matrix(rnorm(n * m), n, m)
  b <- rnorm(subset_frac * m)
  col_index <- sample.int(m, subset_frac * m)
  
  expect_true(are_all_close(
    A[, col_index] %*% b,
    col_subset_matvec(A, b, col_index)
  ))
})