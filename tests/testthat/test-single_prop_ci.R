test_that("Error message should be given if x < 0", {
  expect_error(single_prop_ci(-1, 100))
})

test_that("Error message should be given if x > n", {
  expect_error(single_prop_ci(101, 100))
})
