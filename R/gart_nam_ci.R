gart_nam_ci <- function(x1, n1, x0, n0, skew = TRUE, level = 0.95, alternative = "two.sided") {
  p1hat <- x1 / n1
  p0hat <- x0 / n0

  if (alternative == "two.sided") {
    alpha <- 1 - level
  } else {
    alpha <- (1 - level) * 2
  }

  myfun <- function(phi, skewswitch = skew, lev = level) {
    scorephi(phi = phi, x1 = x1, x0 = x0, n1 = n1, n0 = n0, conf.level = level, alternative = alternative)$score
  }
  mydsct <- function(phi, skewswitch = skew, lev = level) {
    scorephi(phi = phi, x1 = x1, x0 = x0, n1 = n1, n0 = n0, conf.level = level, alternative = alternative)$dsct
  }

  point_FE <- bisect(
    ftn = function(phi) {
      myfun(phi, lev = 0) - 0
    },
    precis = 6,
    uplow = "low"
  )

  point_FE[myfun(phi = 1, lev = 0) == 0] <- 1
  point <- point_FE
  point[myfun(phi = 1, lev = 0) == 0] <- 1


  point[x1 == 0 & x0 == 0] <- NA
  point[(x1 > 0 & x0 == 0) & skew == FALSE] <- Inf
  point[x1 == n1 & x0 == n0] <- 1

  at_MLE <- scorephi(
    phi = ifelse(is.na(point), 1, point),
    x1, n1, x0, n0, skew = TRUE, conf.level = 0.95, alternative = alternative
  )

  p1t_MLE <- at_MLE$p1t
  p0t_MLE <- at_MLE$p0t
  wt_MLE <- at_MLE$wt
  p1hat_w <- p1hat
  p0hat_w <- p0hat
  p1t_w <- p1t_MLE
  p0t_w <- p0t_MLE
  wt_MLE <- NULL

  qtnorm <- qnorm(1 - alpha / 2)

  lower <- bisect(ftn = function(phi) {
    myfun(phi) - qtnorm
  }, precis = 6, uplow = "low")

  upper <- bisect(ftn = function(phi) {
    myfun(phi) + qtnorm
  }, precis = 6, uplow = "up")

  lower[x1 == 0] <- 0
  upper[x0 == 0] <- Inf

  inputs <- cbind(x1 = x1, n1 = n1, x0 = x0, n0 = n0)
  estimates <- list(
    phi_mle = point,
    Lower = lower,
    Upper = upper,
    level = level,
    alternative = alternative,
    inputs,
    p1hat = p1hat_w,
    p0hat = p0hat_w,
    p1mle = p1t_w,
    p0mle = p0t_w,
    V = at_MLE$V
  )
  return(estimates)
}
