#' Estimate the confidence intervals for positive predictive value (PPV) and negative predictive value (NPV) based on different methods
#'
#' @param x1 number of positives for both reference(true) marker and testing marker.
#' @param n1 number of positives for reference (true) marker.
#' @param x0 number of negatives for both reference(true) marker and testing marker.
#' @param n0 number of negatives for reference (true) marker.
#' @param prevalence disease prevalence.
#' @param method current support "gart and nam", "walter", "mover-j", "pepe", "zhou", "delta", "fieller", "bootstrap"; Default is "gart and nam"; Check the \code{Details} for additional information for each method.
#' @param conf.level confidence level. default 0.95.
#' @param bias_correction Logical, indicating whether to apply bias correction in the score denominator. default FALSE. This argument can be used only for `gart and nam` method.
#' @param continuity.correction logical. default FALSE. 0.5 will be applied if TRUE except the \code{zhou}'s method where \eqn{\frac{z_{\alpha/2}^2}{2}} is used.
#' @param n_bootstraps number of bootstraps for bootstrap method, default 1000
#' @param boot.type A character string representing the type of intervals for bootstrap. The value should be any of the values c("bca", "norm", "perc"), default bca.
#' @param ... Other arguments passed on to method (e.g., defining the `Beta(ai,bi)` prior distributions in mover-j method for each group (default `ai = bi = 0.5` for Jeffreys method))
#' @details
#' Eight methods are supported in current version: "gart and nam", "walter", "mover-j", "pepe", "zhou", "delta", "fieller", and "bootstrap".
#'
#' Among those, \strong{gart and nam}, \strong{mover-j}, \strong{walter}, and \strong{fieller} construct the confidence intervals for PPV and NPV by converting the confidence intervals for the ratio of two binomial proportions (\eqn{\phi=\frac{p_1}{p_0}}) where
#'\eqn{\phi_{PPV}=\frac{(1-specificity)}{sensitivity}} and \eqn{\phi_{NPV}=\frac{(1-sensitivity)}{specificity}}.
#'The sensitivity is estimated by \eqn{sensitivity=\frac{x_1}{n_1}} and the specificity by \eqn{specificity=\frac{x_0}{n_0}}. The confidence intervals for \eqn{\phi_{PPV}} and \eqn{\phi_{NPV}}
#'are converted to the corresponding confidence intervals for PPV and NPV using the following equations:
#'\itemize{
#'\item \eqn{PPV=\frac{\rho}{\rho+(1-\rho)*\phi_{PPV}}}
#'
#'\item \eqn{NPV=\frac{1-\rho}{(1-\rho)+\rho*\phi_{NPV}}}
#'
#'where \eqn{\rho} denotes the \code{prevalence}.
#'
#'}
#'\enumerate{
#' \item The \strong{gart and nam} method constructs the confidence interval for \eqn{\phi} based on score method with skewness correction. See the details in the paper listed in the \code{Reference} section.
#' This method can be applied to special situations where `x1=n1` or `x0=n0` but not for `x1=0` and/or `x0=0`. `continuity.correction` can be considered where `x1=0` and/or `x0=0`.
#'
#' \item The \strong{mover-j} method constructs the confidence interval for \eqn{\phi} from separate intervals for the individual group rates (i.e., \eqn{p_1} and \eqn{p_0}).
#' By applying the equal-tailed Jeffreys method (default 0.5 to each group), it may achieve a skewness-corrected interval for \eqn{\phi}. \cr
#'
#' \item The \strong{walter} method constructs the confidence interval for \eqn{\phi} based on \eqn{log(\phi)}. 0.5 is added to `x1`, `x0`, `n1`, and `n0`. Thus, no continuity correct should be applied additionally. This method has shown skewness concerns for
#' small ratios and sample sizes.
#'
#' \item The \strong{fieller} method constructs the confidence interval for \eqn{\phi} by solving the roots of a quadratic equation, see the reference for details. \cr
#' \cr
#' The \strong{pepe}, \strong{zhou}, and \strong{delta} are three direct confidence interval methods for PPV and NPV.
#' \cr
#' \item The \strong{pepe} method finds the confidence intervals for PPV and NPV via diagnostic likelihood ratios (DLR) and associated logit(PPV) and logit(NPV). This method is not applicable for special
#' cases when `x1=0`, `x0=0`, `x1=n1` or `x0=n0`. Continuity correction should be considered for those special cases.
#'
#' \item The \strong{zhou} method can return confidence intervals from the four methods described in the paper. Without continuity correction, it will return the confidence intervals
#' for PPV and NPV based on standard delta method and based on logit transformed method. If \code{continuity.correction=TRUE}, \eqn{\frac{z_{\alpha/2}^2}{2}} will be added to `x1`, `x0`, `n1`, and `n0`, and the function will return the adjusted and adjusted logit confidence intervals
#' as described in the paper.
#'
#' \item The \strong{delta} method constructs the confidence intervals based on the Wald-type formulation. The estimates and variances of PPV and NPV are calculated based equations described in the Zhou's paper listed in the \code{Reference} section.
#'
#' \item The \strong{boot} method constructs the confidence intervals based on bootstrap method.
#'}
#' @return A list object contains the method and the estimates of sensitivity, specificity, PPV, NPV and their confidence intervals.
#' @references
#' \enumerate{
#' \item \insertRef{gn_1988}{CIfinder}
#'
#' \item \insertRef{zhou_2007}{CIfinder}
#'
#' \item \insertRef{laud_2017}{CIfinder}
#'
#' \item \insertRef{pepe_2003}{CIfinder}
#'
#' \item \insertRef{fieller_1954}{CIfinder}
#'}
#' @importFrom ratesci moverci
#' @importFrom stats pnorm quantile
#' @importFrom boot boot boot.ci
#' @import Rdpack kableExtra
#' @export
#' @examples
#' ppv_npv_ci(60, 65, 113, 113, prevalence = 0.02)
ppv_npv_ci <- function (x1, n1, x0, n0, prevalence, method = "gart and nam", conf.level = 0.95, bias_correction = FALSE, continuity.correction = FALSE, n_bootstraps = 1000, boot.type = "bca", ...) {

  if(any(c(x1, n1, x0, n0) < 0)) {
    stop("All input numbers must be non-negative")
  }

  if(n1 < x1 | n0 < x0) {
    stop("n1 or n0 cann't be less than x1 or x0")
  }

  if(prevalence < 0 | prevalence > 1) {
    stop("prevalence must be in [0, 1]")
  }

  if (method == "walter" & continuity.correction) {
    warning("0.5 has been applied as default in Walter method. continuity.correction is not used")
  }

  if (method != "gart and nam" & bias_correction) {
    stop("bias correction is only applicable to gart and nam method")
  }

  prev <- prevalence

  # confidence interval for phi_ppv
  x11 <- x1
  x10 <- n0-x0
  x01 <- n1-x1
  x00 <- x0

  if (method == "gart and nam") {
    out <- get_gn_ci(x10 = x10, x11 = x11, x01 = x01, x00=x00, n1 = n1, n0 = n0, prev = prev, conf.level = conf.level, bias_correction = bias_correction, continuity.correction = continuity.correction)
  }

  if (method == "walter") {
      out <- get_walter_ci(x10 = x10, x11 = x11, x01 = x01, x00 = x00, n1 = n1, n0 = n0, prev = prev, conf.level = conf.level)
  }

  if (method == "pepe") {
      out <- get_pepe_ci(x11 = x11, x10 = x10, x01 = x01, x00 = x00, n1 = n1, n0 = n0, prev = prev, conf.level = conf.level, continuity.correction = continuity.correction)
  }

  if (method == "mover-j") {
    out <- get_moverj_ci(x10 = x10, x11 = x11, x01 = x01, x00 = x00, n1 = n1, n0 = n0, prev =prev, conf.level = conf.level, continuity.correction = continuity.correction, ...)
  }

  if (method == "zhou") {
    out <- get_zhou_ci(x10 = x10, x11 = x11, x01 = x01, x00 = x00, n1 = n1, n0 = n0, prev = prev, conf.level = conf.level, continuity.correction = continuity.correction)
  }

  if (method == "delta") {
    out <- get_delta_ci(x10 = x10, x11 = x11, x01 = x01, x00 = x00, n1 = n1, n0 = n0, conf.level = conf.level, continuity.correction = continuity.correction, prev = prev)
  }

  if (method == "fieller") {
    out <- get_fieller_ci(x10 = x10, x11 = x11, x01 = x01, x00 = x00, n1 = n1, n0 = n0, conf.level = conf.level, continuity.correction = continuity.correction, prev = prev)
  }

  if (method == "boot") {
    out <- get_boot_ci(x10 = x10, x11 = x11, x01 = x01, x00 = x00, n1 = n1, n0 = n0, conf.level = conf.level, continuity.correction = continuity.correction, prev = prev, n_bootstraps = n_bootstraps, type = boot.type)
  }
  return(out)
}

get_gn_ci <- function (x10, x11, x01, x00, n1, n0, conf.level, prev, bias_correction, continuity.correction) {

  if (x11==0 | x00==0) {
    continuity.correction = TRUE
    warning("x1 or x0 equals zero, continuity correction is applied")
  }

  # sen, spe, phi, etc
  sen <- x11/n1
  spe <- x00/n0
  phi_ppv <- (1-spe)/sen
  phi_npv <- (1-sen)/spe
  ppv_est <- prev/(prev + (1 - prev) * phi_ppv)
  npv_est <- (1 - prev)/(1 - prev + prev * phi_npv)
  ##calculate the confidence interval for phi_ppv = (1-spe)/sen
  fit_phi_ppv <- ratesci::scoreci(x1 = x10, n1 = n0, x2 = x11, n2 = n1, distrib = "bin",
                                  contrast = "RR", level = conf.level, bcf = bias_correction, cc = continuity.correction)
  phi_ppv_mle <- fit_phi_ppv$estimates[2]
  phi_ppv_l <- fit_phi_ppv$estimates[1]
  phi_ppv_u <- fit_phi_ppv$estimates[3]
  ppv_mle <- prev/(prev + (1 - prev) * fit_phi_ppv$estimates[2])
  ppv_u <- prev/(prev + (1 - prev) * fit_phi_ppv$estimates[1])
  ppv_l <- prev/(prev + (1 - prev) * fit_phi_ppv$estimates[3])

  ##calculate the confidence interval for phi_npv = (1-sen)/spe
  fit_phi_npv <- ratesci::scoreci(x1 = x01, n1 = n1, x2 = x00, n2 = n0,
                                  distrib = "bin", contrast = "RR", level = conf.level, bcf = bias_correction, cc = continuity.correction)
  phi_npv_mle <- fit_phi_npv$estimates[2]
  phi_npv_l <- fit_phi_npv$estimates[1]
  phi_npv_u <- fit_phi_npv$estimates[3]
  npv_mle <- (1 - prev)/(1 - prev + prev * fit_phi_npv$estimates[2])
  npv_l <- (1 - prev)/(1 - prev + prev * fit_phi_npv$estimates[3])
  npv_u <- (1 - prev)/(1 - prev + prev * fit_phi_npv$estimates[1])

  return(list(method = "gart and nam",
              sensitivity = sen,
              specificity = spe,
              phi_ppv = c(phi_ppv_est = phi_ppv,
                          phi_ppv_l = phi_ppv_l,
                          phi_ppv_u = phi_ppv_u,
                          phi_ppv_mle = phi_ppv_mle),
              ppv=c(ppv_est = ppv_est,
                    ppv_l = ppv_l,
                    ppv_u = ppv_u,
                    ppv_mle=ppv_mle),
              phi_npv = c(phi_npv_est = phi_npv,
                          phi_npv_l = phi_npv_l,
                          phi_npv_u = phi_npv_u,
                          phi_npv_mle = phi_npv_mle),
              npv=c(npv_est = npv_est,
                    npv_l = npv_l,
                    npv_u = npv_u,
                    npv_mle=npv_mle)))
}

get_pepe_ci <- function(x11, x10, x01, x00, n1, n0, prev, conf.level, continuity.correction) {
  z <- qnorm(1-(1-conf.level)/2)
  if (continuity.correction) {
    x11 <- x11 + 0.5
    x10 <- x10 + 0.5
    x01 <- x01 + 0.5
    x00 <- x00 + 0.5
    n1 <- x11 + x01
    n0 <- x10 + x00
  }

  sen <- x11 / n1
  spe <- x00 / n0
  phi_ppv <- (1-spe)/sen
  phi_npv <- (1-sen)/spe
  ppv_est <- prev/(prev + (1 - prev) * phi_ppv)
  npv_est <- (1 - prev)/(1 - prev + prev * phi_npv)

  dlr1 <- sen / (1 - spe)
  dlr0 <- (1 - sen) / spe

  var_log_dlr1 <- function(sen, spe, n1, n0) {
    (1 - sen) / (n1 * sen) + spe / (n0 * (1 - spe))
  }

  var_log_dlr0 <- function(sen, spe, n1, n0) {
    sen / (n1 * (1 - sen)) + (1 - spe) / (n0 * spe)
  }

  log_dlr1_u <- log(dlr1) + z * sqrt(var_log_dlr1(sen, spe, n1, n0))
  log_dlr1_l <- log(dlr1) - z * sqrt(var_log_dlr1(sen, spe, n1, n0))

  dlr1_ci <- exp(c(log_dlr1_l, log_dlr1_u))
  log_dlr1_ci <- c(log_dlr1_l, log_dlr1_u)
  logit_ppv_ci <- log(prev / (1 - prev)) + log_dlr1_ci
  ppv_ci <- exp(logit_ppv_ci) / (1 + exp(logit_ppv_ci))

  log_dlr0_u <- log(dlr0) + z * sqrt(var_log_dlr0(sen, spe, n1, n0))
  log_dlr0_l <- log(dlr0) - z * sqrt(var_log_dlr0(sen, spe, n1, n0))
  dlr0_ci <- exp(c(log_dlr0_l, log_dlr0_u))
  log_dlr0_ci <- c(log_dlr0_l, log_dlr0_u)
  logit_npv_ci <- (-log(prev / (1 - prev))) - log_dlr0_ci
  npv_ci <- 1/ (1 + exp(-logit_npv_ci))

  return(list(method = "pepe",
              sensitivity = sen,
              specificity = spe,
              log_dlr_pos = c(log_dlr_pos = log(dlr1),
                              log_dlr_pos_l = log_dlr1_l,
                              log_dlr_pos_u = log_dlr1_u),
              dlr_pos = c(dlr_pos = dlr1,
                          dlr_pos_l = dlr1_ci[1],
                          dlr_pos_u = dlr1_ci[2]),
              ppv=c(ppv_est = ppv_est,
                    ppv_l = ppv_ci[1],
                    ppv_u = ppv_ci[2]),
              log_dlr_neg = c(log_dlr_neg = log(dlr0),
                              log_dlr_neg_l = log_dlr0_l,
                              log_dlr_neg_u = log_dlr0_u),
              dlr_neg = c(dlr_neg = dlr0,
                          dlr_neg_l = dlr0_ci[1],
                          dlr_neg_u = dlr0_ci[2]),
              npv=c(npv_est = npv_est,
                    npv_l = npv_ci[2],
                    npv_u = npv_ci[1])))
}

get_moverj_ci <- function(x10, x11, x01, x00, n1, n0, prev, conf.level, continuity.correction, ...) {
  # sen, spe, phi, etc
  sen <- x11/n1
  spe <- x00/n0
  phi_ppv <- (1-spe)/sen
  phi_npv <- (1-sen)/spe
  ppv_est <- prev/(prev + (1 - prev) * phi_ppv)
  npv_est <- (1 - prev)/(1 - prev + prev * phi_npv)

  phi_ppv_ests <- moverci(x10, n0, x11, n1, level = conf.level, contrast = "RR", distrib = "bin", type = "jeff", cc = continuity.correction, ...)
  phi_ppv_l <- phi_ppv_ests[1]
  phi_ppv_u <- phi_ppv_ests[3]
  ppv_u <- prev/(prev + (1 - prev) * phi_ppv_l)
  ppv_l <- prev/(prev + (1 - prev) * phi_ppv_u)

  phi_npv_ests <- moverci(x01, n1, x00, n0, level = conf.level, contrast = "RR", distrib = "bin", type = "jeff", cc = continuity.correction, ...)
  phi_npv_l <- phi_npv_ests[1]
  phi_npv_u <- phi_npv_ests[3]
  npv_l <- (1 - prev)/(1 - prev + prev * phi_npv_u)
  npv_u <- (1 - prev)/(1 - prev + prev * phi_npv_l)
  return(list(method = "mover-j",
              sensitivity = sen,
              specificity = spe,
              phi_ppv = c(phi_ppv_est = phi_ppv_ests[2],
                          phi_ppv_l = phi_ppv_l,
                          phi_ppv_u = phi_ppv_u),
              ppv=c(ppv_est = ppv_est,
                    ppv_l = ppv_l,
                    ppv_u = ppv_u),
              phi_npv = c(phi_npv_est = phi_npv_ests[2],
                          phi_npv_l = phi_npv_l,
                          phi_npv_u = phi_npv_u),
              npv=c(npv_est = npv_est,
                    npv_l = npv_l,
                    npv_u = npv_u)))
}

get_zhou_ci <- function(x10, x11, x01, x00, n1, n0, prev, conf.level, continuity.correction) {
  z <- qnorm(1-(1-conf.level)/2)
  if (continuity.correction) {
    k = qnorm(1-(1-conf.level)/2)
    c.correct = k^2/2
    x11 <- x11 + c.correct
    x10 <- x10 + c.correct
    x01 <- x01 + c.correct
    x00 <- x00 + c.correct
    n1 <- x11 + x01
    n0 <- x10 + x00
  }

  ## sen, spe and vars
  sen <- x11 / n1
  spe <- x00 / n0
  sen_var <- (sen*(1-sen))/n1
  spe_var <- (spe*(1-spe))/n0

  ## delta method
  ppv_est <-(sen*prev)/(sen*prev + (1-prev)*(1-spe))
  npv_est <- ((1-prev)*spe)/(prev * (1-sen) + (1-prev) * spe)
  ppv_var <- ((prev*(1-spe)*(1-prev))^2*sen_var + (prev*sen*(1-prev))^2*spe_var)/(sen*prev + (1-spe)*(1-prev))^4
  npv_var <- ((spe*(1-prev)*prev)^2*sen_var + ((1-sen)*(1-prev)*prev)^2*spe_var)/((1-sen)*prev+spe*(1-prev))^4
  ppv_u <- ppv_est + z * sqrt(ppv_var)
  ppv_l <- ppv_est - z * sqrt(ppv_var)
  npv_u <- npv_est + z * sqrt(npv_var)
  npv_l <- npv_est - z * sqrt(npv_var)

  # logit transformation
  logit_ppv_est <- log((sen*prev)/((1-spe)*(1-prev)))
  logit_npv_est <- log((spe*(1-prev))/((1-sen)*prev))
  logit_ppv_var <- ((1-sen)/sen)/n1 + (spe/(1-spe))/n0
  logit_npv_var <- (sen/(1-sen))/n1 + ((1-spe)/spe)/n0
  ppv_logit_transformed_u <- exp(logit_ppv_est + z * sqrt(logit_ppv_var))/(1 + exp(logit_ppv_est + z * sqrt(logit_ppv_var)))
  ppv_logit_transformed_l <- exp(logit_ppv_est - z * sqrt(logit_ppv_var))/(1 + exp(logit_ppv_est - z * sqrt(logit_ppv_var)))
  npv_logit_transformed_u <- exp(logit_npv_est + z * sqrt(logit_npv_var))/(1 + exp(logit_npv_est + z * sqrt(logit_npv_var)))
  npv_logit_transformed_l <- exp(logit_npv_est - z * sqrt(logit_npv_var))/(1 + exp(logit_npv_est - z * sqrt(logit_npv_var)))

  return(list(method = "zhou",
              sensitivity = sen,
              specificity = spe,
              ppv=c(ppv_est = ppv_est,
                    ppv_l = ppv_l,
                    ppv_u = ppv_u),
              npv=c(npv_est = npv_est,
                    npv_l = npv_l,
                    npv_u = npv_u),
              ppv_logit_transformed=c(ppv_est = ppv_est,
                                      ppv_l = ppv_logit_transformed_l,
                                      ppv_u = ppv_logit_transformed_u),
              npv_logit_transformed=c(npv_est = npv_est,
                                      npv_l = npv_logit_transformed_l,
                                      npv_u = npv_logit_transformed_u)))
}

get_delta_ci <- function(x10, x11, x01, x00, n1, n0, prev, conf.level, continuity.correction) {
  z <- qnorm(1-(1-conf.level)/2)

  if (continuity.correction) {
    x11 <- x11 + 0.5
    x10 <- x10 + 0.5
    x01 <- x01 + 0.5
    x00 <- x00 + 0.5
    n1 <- x11 + x01
    n0 <- x10 + x00
  }

  ## sen, spe and vars
  sen <- x11 / n1
  spe <- x00 / n0
  sen_var <- (sen*(1-sen))/n1
  spe_var <- (spe*(1-spe))/n0

  ## delta method
  ppv_est <-(sen*prev)/(sen*prev + (1-prev)*(1-spe))
  npv_est <- ((1-prev)*spe)/(prev * (1-sen) + (1-prev) * spe)
  ppv_var <- ((prev*(1-spe)*(1-prev))^2*sen_var + (prev*sen*(1-prev))^2*spe_var)/(sen*prev + (1-spe)*(1-prev))^4
  npv_var <- ((spe*(1-prev)*prev)^2*sen_var + ((1-sen)*(1-prev)*prev)^2*spe_var)/((1-sen)*prev+spe*(1-prev))^4
  ppv_u <- ppv_est + z * sqrt(ppv_var)
  ppv_l <- ppv_est - z * sqrt(ppv_var)
  npv_u <- npv_est + z * sqrt(npv_var)
  npv_l <- npv_est - z * sqrt(npv_var)

  return(list(method = "delta",
              sensitivity = sen,
              specificity = spe,
              ppv=c(ppv_est = ppv_est,
                    ppv_l = ppv_l,
                    ppv_u = ppv_u),
              npv=c(npv_est = npv_est,
                    npv_l = npv_l,
                    npv_u = npv_u)))
}

walter_ci <- function(x1, n1, x0, n0, conf.level) {
  z <- qnorm(1-(1-conf.level)/2)
  p0 <- x0 / n0
  p1 <- x1 / n1
  q0 <- 1 - p0
  q1 <- 1 - p1
  phi <- p1 / p0

  log_phi <- log((x1 + 0.5) / (n1 + 0.5)) - log((x0 + 0.5) / (n0 + 0.5))
  u <- 1 / (x1 + 0.5) - 1 / (n1 + 0.5) + 1 / (x0 + 0.5) - 1 / (n0 + 0.5)
  lower <- exp(log_phi) * exp(-z * sqrt(u))
  upper <- exp(log_phi) * exp(z * sqrt(u))
  return(c(phi = phi, lower = lower, upper = upper))
}

get_walter_ci <- function(x10, x11, x01, x00, n1, n0, prev, conf.level) {
  sen <- x11/n1
  spe <- x00/n0
  phi_ppv <- (1-spe)/sen
  phi_npv <- (1-sen)/spe
  ppv_est <- prev/(prev + (1 - prev) * phi_ppv)
  npv_est <- (1 - prev)/(1 - prev + prev * phi_npv)
  ##calculate the confidence interval for phi_ppv = (1-spe)/sen
  fit_phi_ppv <- walter_ci(x1 = x10, n1 = n0, x0 = x11, n0 = n1, conf.level = conf.level)
  ppv_u <- prev/(prev + (1 - prev) * fit_phi_ppv[2])
  ppv_l <- prev/(prev + (1 - prev) * fit_phi_ppv[3])

  ##calculate the confidence interval for phi_npv = (1-sen)/spe
  fit_phi_npv <- walter_ci(x1 = x01, n1 = n1, x0 = x00, n0 = n0, conf.level = conf.level)
  npv_u <- (1 - prev)/(1 - prev + prev * fit_phi_npv[2])
  npv_l <- (1 - prev)/(1 - prev + prev * fit_phi_npv[3])

  return(list(method = "walter",
              sensitivity = sen,
              specificity = spe,
              phi_ppv = c(fit_phi_ppv[1],
                          fit_phi_ppv[2],
                          fit_phi_ppv[3]),
              ppv=c(ppv_est = ppv_est,
                    ppv_l = unname(ppv_l),
                    ppv_u = unname(ppv_u)),
              phi_npv = c(fit_phi_npv[1],
                          fit_phi_npv[2],
                          fit_phi_npv[3]),
              npv=c(npv_est = npv_est,
                    npv_l = unname(npv_l),
                    npv_u = unname(npv_u))))
}

## return NA if there are no solutions
solve_quadratic <- function(a, b, c) {
  tryCatch({
    # Calculate the discriminant
    D <- b^2 - 4*a*c

    # Check the discriminant to determine the nature of the roots
    if (D > 0) {
      # Two distinct real roots
      root1 <- (-b + sqrt(D)) / (2 * a)
      root2 <- (-b - sqrt(D)) / (2 * a)
      return(c(root1, root2))
    } else if (D == 0) {
      # One real root (roots are real and equal)
      root <- -b / (2 * a)
      return(c(root, root))
    } else {
      # No real roots, issue a warning
      warning("The equation has no real roots.")
      return(c(NA, NA))
    }
  }, error = function(e) {
    # Handle unexpected errors
    message("Error: ", e$message)
    return(c(NA, NA))
  })
}

get_fieller_ci <- function(x10, x11, x01, x00, n1, n0, prev, conf.level, continuity.correction) {
  z <- qnorm(1-(1-conf.level)/2)
  if (continuity.correction) {
    x11 <- x11 + 0.5
    x10 <- x10 + 0.5
    x01 <- x01 + 0.5
    x00 <- x00 + 0.5
    n1 <- x11 + x01
    n0 <- x10 + x00
  }

  ## sen, spe and vars
  sen <- x11 / n1
  spe <- x00 / n0
  sen_var <- (sen*(1-sen))/n1
  spe_var <- (spe*(1-spe))/n0

  phi_ppv <- (1-spe)/sen
  phi_npv <- (1-sen)/spe
  ppv_est <- prev/(prev + (1 - prev) * phi_ppv)
  npv_est <- (1 - prev)/(1 - prev + prev * phi_npv)

  ## fieller method for phi_s
  a1 <- sen^2 - z^2*sen_var
  b1 <- -2*sen*(1-spe)
  c1 <- (1-spe)^2-z^2*spe_var

  a2 <- spe^2 - z^2*spe_var
  b2 <- -2*spe*(1-sen)
  c2 <- (1-sen)^2-z^2*sen_var

  if (a1 <= 0 || b1^2 - 4 * a1 * c1 < 0 || a2 <= 0 || b2^2 - 4 * a2 * c2 < 0) {
    warning("Conditions for Fieller's method are not met.")
    return(list(method = "fieller",
                sensitivity = sen,
                specificity = spe,
                phi_ppv = c(phi_ppv_est = phi_ppv,
                            phi_ppv_l = NA,
                            phi_ppv_u = NA),
                ppv=c(ppv_est = ppv_est,
                      ppv_l = NA,
                      ppv_u = NA),
                phi_npv = c(phi_npv_est = phi_npv,
                            phi_npv_l = NA,
                            phi_npv_u = NA),
                npv=c(npv_est = npv_est,
                      npv_l = NA,
                      npv_u = NA)))
  } else {
    # phi_ppv_l <- (-b1 - sqrt(b1^2 - 4 * a1 * c1)) / (2 * a1)
    # phi_ppv_u <- (-b1 + sqrt(b1^2 - 4 * a1 * c1)) / (2 * a1)
    #
    # phi_npv_l <- (-b2 - sqrt(b2^2 - 4 * a2 * c2)) / (2 * a2)
    # phi_npv_u <- (-b2 + sqrt(b2^2 - 4 * a2 * c2)) / (2 * a2)
    #
    # # Calculate confidence intervals for PPV and NPV
    # ppv_l <- prev / (prev + (1 - prev) * phi_ppv_u)
    # ppv_u <- prev / (prev + (1 - prev) * phi_ppv_l)
    #
    # npv_l <- (1 - prev) / (1 - prev + prev * phi_npv_u)
    # npv_u <- (1 - prev) / (1 - prev + prev * phi_npv_l)

    ps1 <- min(solve_quadratic(a1, b1, c1))
    ps2 <- max(solve_quadratic(a1, b1, c1))
    phi_ppv_l <- ifelse(ps1 < 0, 0, ps1)
    phi_ppv_u <- ifelse(ps2 < 0, 0, ps2)
    ppv_u <- ifelse(is.na(phi_ppv_l), NA, prev/(prev + (1 - prev) * phi_ppv_l))
    ppv_l <- ifelse(is.na(phi_ppv_u), NA, prev/(prev + (1 - prev) * phi_ppv_u))


    ns1 <- min(solve_quadratic(a2, b2, c2))
    ns2 <- max(solve_quadratic(a2, b2, c2))
    phi_npv_l <- ifelse(ns1 < 0, 0, ns1)
    phi_npv_u <- ifelse(ns2 < 0, 0, ns2)
    npv_u <- ifelse(is.na(phi_npv_l), NA, (1 - prev)/(1 - prev + prev * phi_npv_l))
    npv_l <- ifelse(is.na(phi_npv_u), NA, (1 - prev)/(1 - prev + prev * phi_npv_u))

    return(list(method = "fieller",
                sensitivity = sen,
                specificity = spe,
                phi_ppv = c(phi_ppv_est = phi_ppv,
                            phi_ppv_l = phi_ppv_l,
                            phi_ppv_u = phi_ppv_u),
                ppv=c(ppv_est = ppv_est,
                      ppv_l = ppv_l,
                      ppv_u = ppv_u),
                phi_npv = c(phi_npv_est = phi_npv,
                            phi_npv_l = phi_npv_l,
                            phi_npv_u = phi_npv_u),
                npv=c(npv_est = npv_est,
                      npv_l = npv_l,
                      npv_u = npv_u)))
  }
}

get_boot_ci <- function(x10, x11, x01, x00, n1, n0, prev, conf.level, continuity.correction, n_bootstraps, type) {

  ## sen, spe and vars
  sen <- x11 / n1
  spe <- x00 / n0
  phi_ppv <- (1-spe)/sen
  phi_npv <- (1-sen)/spe
  ppv_est <- prev/(prev + (1 - prev) * phi_ppv)
  npv_est <- (1 - prev)/(1 - prev + prev * phi_npv)

  # Convert to a four-column binary matrix
  .dat <- data.frame(
    disease = c(rep(1, x11), rep(0, x10), rep(1, x01), rep(0, x00)),
    test = c(rep(1, x11), rep(1, x10), rep(0, x01), rep(0, x00))
  )

  # Function to calculate phi_ppv and phi_npv
  calculate_metrics <- function(data, indices) {
    sample_data <- data[indices, ]

    # Calculate sensitivity, specificity
      x11_corr <- sum(sample_data$disease == 1 & sample_data$test == 1)
      x10_corr <- sum(sample_data$disease == 0 & sample_data$test == 1)
      x01_corr <- sum(sample_data$disease == 1 & sample_data$test == 0)
      x00_corr <- sum(sample_data$disease == 0 & sample_data$test == 0)
      # Calculate sensitivity, specificity with correction
      sensitivity <- x11_corr / (x11_corr + x01_corr)
      specificity <- x00_corr / (x00_corr + x10_corr)

    sensitivity <- sum(sample_data$disease == 1 & sample_data$test == 1) / sum(sample_data$disease == 1)
    specificity <- sum(sample_data$disease == 0 & sample_data$test == 0) / sum(sample_data$disease == 0)

    # Avoid division by zero
    sensitivity <- ifelse(sensitivity == 0, 1e-8, sensitivity)
    specificity <- ifelse(specificity == 0, 1e-8, specificity)

    # Calculate phi_ppv and phi_npv
    phi_ppv <- (1 - specificity) / sensitivity
    phi_npv <- (1 - sensitivity) / specificity

    return(c(phi_ppv, phi_npv))
  }

  # Function to calculate phi_ppv and phi_npv with continuity correction
  calculate_metrics_cc <- function(data, indices) {
    sample_data <- data[indices, ]

    # Calculate sensitivity, specificity
    x11_corr <- sum(sample_data$disease == 1 & sample_data$test == 1) + 0.5
    x10_corr <- sum(sample_data$disease == 0 & sample_data$test == 1) + 0.5
    x01_corr <- sum(sample_data$disease == 1 & sample_data$test == 0) + 0.5
    x00_corr <- sum(sample_data$disease == 0 & sample_data$test == 0) + 0.5
    # Calculate sensitivity, specificity with correction
    sensitivity <- x11_corr / (x11_corr + x01_corr)
    specificity <- x00_corr / (x00_corr + x10_corr)

    sensitivity <- sum(sample_data$disease == 1 & sample_data$test == 1) / sum(sample_data$disease == 1)
    specificity <- sum(sample_data$disease == 0 & sample_data$test == 0) / sum(sample_data$disease == 0)

    # Avoid division by zero
    sensitivity <- ifelse(sensitivity == 0, 1e-8, sensitivity)
    specificity <- ifelse(specificity == 0, 1e-8, specificity)

    # Calculate phi_ppv and phi_npv
    phi_ppv <- (1 - specificity) / sensitivity
    phi_npv <- (1 - sensitivity) / specificity

    return(c(phi_ppv, phi_npv))
  }

  # Bootstrap with the boot package
  if(continuity.correction) {
    bootstrap_results <- boot(.dat, calculate_metrics_cc, R = n_bootstraps)
  } else {
    bootstrap_results <- boot(.dat, calculate_metrics, R = n_bootstraps)
  }

  # Calculating the BCa confidence interval for both metrics
  safe_boot_ci <- function(boot_obj, type = type, conf.level) {
    tryCatch({
      boot_ci_phi_ppv <- boot.ci(boot_obj, type = type, conf = conf.level, index = 1)
      boot_ci_phi_npv <- boot.ci(boot_obj, type = type, conf = conf.level, index = 2)
      # Displaying the BCa confidence intervals
      if (type =="bca"){
        phi_ppv_l <- boot_ci_phi_ppv$bca[4]
        phi_ppv_u <- boot_ci_phi_ppv$bca[5]

        phi_npv_l <- boot_ci_phi_npv$bca[4]
        phi_npv_u <- boot_ci_phi_npv$bca[5]
      }

      if (type =="norm"){
        phi_ppv_l <- boot_ci_phi_ppv$normal[2]
        phi_ppv_u <- boot_ci_phi_ppv$normal[3]

        phi_npv_l <- boot_ci_phi_npv$normal[2]
        phi_npv_u <- boot_ci_phi_npv$normal[3]
      }

      if (type =="perc"){
        phi_ppv_l <- boot_ci_phi_ppv$percent[4]
        phi_ppv_u <- boot_ci_phi_ppv$percent[5]

        phi_npv_l <- boot_ci_phi_npv$percent[4]
        phi_npv_u <- boot_ci_phi_npv$percent[5]
      }

      ppv_u <- prev/(prev + (1 - prev) * phi_ppv_l)
      ppv_l <- prev/(prev + (1 - prev) * phi_ppv_u)

      npv_l <- (1 - prev)/(1 - prev + prev * phi_npv_u)
      npv_u <- (1 - prev)/(1 - prev + prev * phi_npv_l)
      return(c(phi_ppv_l,
             phi_ppv_u,
             phi_npv_l,
             phi_npv_u,
             ppv_l,
             ppv_u,
             npv_l,
             npv_u))
    }, error = function(e) {
      # Handle the error gracefully
      warning("Failed to compute CI: ", e$message)
      return(rep(NA, 8))
    })
  }

  boot_ci <- safe_boot_ci(boot_obj = bootstrap_results, type = type, conf.level = conf.level)

  return(list(method = "boot",
              sensitivity = sen,
              specificity = spe,
              phi_ppv = c(phi_ppv_est = phi_ppv,
                          phi_ppv_l = boot_ci[1],
                          phi_ppv_u = boot_ci[2]),
              ppv=c(ppv_est = ppv_est,
                    ppv_l = boot_ci[5],
                    ppv_u = boot_ci[6]),
              phi_npv = c(phi_npv_est = phi_npv,
                          phi_npv_l = boot_ci[3],
                          phi_npv_u = boot_ci[4]),
              npv=c(npv_est = npv_est,
                    npv_l = boot_ci[7],
                    npv_u = boot_ci[8])))
}

