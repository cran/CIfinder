test_that("Error message should be given for gart and name method if the input x1=0 and/or x0=0 and continuity correction is FALSE", {
  expect_error(ppv_npv_ci(0, 100, 99, 100, prevalence = 0.1, method = "gart and nam", continuity.correction = FALSE))
})

test_that("Error message should be given if bias_correction is set except gart and nam method", {
  expect_error(ppv_npv_ci(88, 100, 99, 100, prevalence = 0.1, method = "mover-j", bias_correction = TRUE))
})

test_that("Warning message should be given for walter method if continuit correction is set to TRUE", {
  expect_warning(ppv_npv_ci(88, 100, 99, 100, prevalence = 0.1, method = "walter", continuity.correction = TRUE))
})
