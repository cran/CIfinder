bca_ci <- function(theta, alternative = "two.sided", conf.level = 0.95) {
  if (alternative == "two.sided") {
    alpha <- 1- conf.level
  } else {
    alpha <- (1 - conf.level) * 2
  }

  low <- alpha/2
  high <- 1- alpha/2
  sims <- length(theta)

  z.inv <- length(theta[theta < mean(theta)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(theta) - theta)
  top <- sum(U^3)
  under <- 6 * (sum(U^2))^{3/2}
  a <- top/under
  lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower <- quantile(theta, lower.inv)
  upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper <- quantile(theta, upper.inv)
  return(c(lower, upper))
}
