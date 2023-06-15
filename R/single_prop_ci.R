#' Compute the confidence interval for a single proportion based on different methods
#'
#' @param x number of successes
#' @param n number of trials
#' @param method one of these options "all", "clopper.pearson", "wald", "wislon", "wislon.correct", "agresti", or "beta"
#' @param alternative indicates "two.sided", "one.sided"
#' @param conf.level confidence level
#' @param prior the prior values for "beta" method
#'
#' @return Estimated confidence intervals for the probability of success
#' @importFrom stats qbeta binom.test qnorm prop.test
#' @export
#'
#' @examples
#' single_prop_ci(53, 57, method = "all")
single_prop_ci <- function(x, n, method = "all", alternative = "two.sided", conf.level=0.95, prior = c(1, 1)) {
  if (method == "all"){
    return(list(clopper_ci = clopper_ci(x, n, test = alternative, conf.level = conf.level),
                wilson = wilson_ci(x, n, test = alternative, conf.level = conf.level),
                wilson.correct =wilson.correct_ci(x, n, test = alternative, conf.level = conf.level),
                wald = wald_ci(x, n, test = alternative, conf.level = conf.level),
                agresti = agresti_ci(x, n, test = alternative, conf.level = conf.level),
                beta_binomial = beta_binomial(x, n, test = alternative, conf.level = conf.level, prior = prior)))
  }
  outrange <- x > n | x < 0
  if (any(outrange)) {
    stop("x must be non-nagative and less or equal to n;")
  }
  if (method == "clopper.pearson"){
    return(clopper_ci = clopper_ci(x, n, test = alternative, conf.level = conf.level))
  }
  if (method == "wilson"){
    return(wilson = wilson_ci(x, n, test = alternative, conf.level = conf.level))
  }
  if (method == "wilson.correct"){
    return(wilson.correct = wilson.correct_ci(x, n, test = alternative, conf.level = conf.level))
  }
  if (method == "wald"){
    return(wald = wald_ci(x, n, test = alternative, conf.level = conf.level))
  }
  if (method == "agresti"){
    return(agresti = agresti_ci(x, n, test = alternative, conf.level = conf.level))
  }
  if (method == "beta"){
    return(beta_binomial = beta_binomial(x, n, test = alternative, conf.level = conf.level, prior = prior))
  }
}


wilson_ci <- function(x, n, test = "two.sided", conf.level = 0.95) {
  limits <- prop.test(x, n, alternative = test, conf.level = conf.level, correct = FALSE)$conf.int
  limits <- c(limits[1], limits[2])
  names(limits) <- c('wilson_lower','wilson_upper')
  limits
}

wilson.correct_ci <- function(x, n, test = "two.sided", conf.level = 0.95) {
  limits <- prop.test(x, n, alternative = test, conf.level = conf.level, correct = TRUE)$conf.int
  limits <- c(limits[1], limits[2])
  names(limits) <- c('wilson.correct_lower','wilson.correct_upper')
  limits
}


wald_ci <- function(x, n, test = "two.sided", conf.level = 0.95) {
  p <- x/n
  se <- sqrt(p * (1 - p)/n)
  if (test == "two.sided") {
    conf.level <- conf.level
    z <- c(qnorm((1 - conf.level)/2), qnorm(1 - (1 - conf.level)/2))
    limits <- p + z * se
    names(limits) <- c('wald_lower','wald_upper')
    limits
  } else {
    conf.level <- 1-(1-conf.level)*2
    z <- c(qnorm((1 - conf.level)/2), qnorm(1 - (1 - conf.level)/2))
    limits <- p + z * se
    if (test == "less"){
      ci <- c(0, max(limits))
      names(ci) <- c('wald_lower','wald_upper')
      return(ci)
    }
    if (test == "greater"){
      ci <- c(min(limits), 1)
      names(ci) <- c('wald_lower','wald_upper')
      return(ci)
    }
  }
}

agresti_ci <- function(x, n, test = "two.sided", conf.level = 0.95) {
  if (test == "two.sided") {
    conf.level <- conf.level
    z <- c(qnorm((1 - conf.level)/2), qnorm(1 - (1 - conf.level)/2))
    n_tilde <- n + (z[1])^2
    p_tilde <- (1/n_tilde)*(x + (z[1])^2/2)
    se <- sqrt((p_tilde/n_tilde)*(1-p_tilde))
    limits <- p_tilde + z * se
    names(limits) <- c('agresti_lower','agresti_upper')
    limits
  } else {
    conf.level <- 1-(1-conf.level)*2
    z <- c(qnorm((1 - conf.level)/2), qnorm(1 - (1 - conf.level)/2))
    n_tilde <- n + (z[1])^2
    p_tilde <- (1/n_tilde)*(x + (z[1])^2/2)
    se <- sqrt((p_tilde/n_tilde)*(1-p_tilde))
    limits <- p_tilde + z * se
    if (test == "less"){
      ci <- c(0, max(limits))
      names(ci) <- c('agresti_lower','agresti_upper')
      return(ci)
    }
    if (test == "greater"){
      ci <- c(min(limits), 1)
      names(ci) <- c('agresti_lower','agresti_upper')
      return(ci)
    }
  }
}

clopper_ci <- function(x, n, test = "two.sided", conf.level=0.95) {
  limits <- as.numeric(binom.test(x, n, alternative = test, conf.level=conf.level)$conf.int)
  names(limits) <- c('cp_lower','cp_upper')
  limits
}


beta_binomial <- function(x, n, test = "two.sided", conf.level = 0.95, prior = c(1,1)){
  if (test == "two.sided") {
    conf.level <- conf.level
    limit <- (1-conf.level)/2
    limits <- qbeta(c(limit, 1-limit), x+prior[1],n-x+prior[2])
    names(limits) <- c('beta_lower','beta_upper')
    limits
  } else {
    conf.level <- 1-(1-conf.level)*2
    limit <- (1-conf.level)/2
    limits <- qbeta(c(limit, 1-limit), x+prior[1],n-x+prior[2])
    if (test == "less"){
      ci <- c(0, max(limits))
      names(ci) <- c('beta_lower','beta_upper')
      return(ci)
    }
    if (test == "greater"){
      ci <- c(min(limits), 1)
      names(ci) <- c('beta_lower','beta_upper')
      return(ci)
    }
  }
}

