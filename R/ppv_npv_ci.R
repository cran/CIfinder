#' Estimate the confidence intervals for positive predictive value (PPV) and negative predictive value (NPV) based on different methods
#'
#' @param x1 number of positives for both reference(true) marker and testing marker.
#' @param n1 number of positives for reference (true) marker.
#' @param x0 number of negatives for both reference(true) marker and testing marker.
#' @param n0 number of positives for reference (true) marker.
#' @param prevalence disease prevalence.
#' @param method current support "gart and nam", "walter", "mover-j", "pepe", "zhou", or "delta"; Default is "gart and nam"; Check the \code{Details} for additional information for each method.
#' @param conf.level confidence level. default 0.95.
#' @param bias_correction Logical, indicating whether to apply bias correction in the score denominator. default FALSE. This argument can be used only for `gart and nam` method.
#' @param continuity.correction logical. default FALSE. 0.5 will be applied if TRUE except the \code{zhou}'s method where \eqn{\frac{z_{\alpha/2}^2}{2}} is used.
#' @param ... Other arguments passed on to method (e.g., defining the `Beta(ai,bi)` prior distributions in mover-j method for each group (default `ai = bi = 0.5` for Jeffreys method))
#' @details
#' Six methods are supported in current version: "gart and nam", "walter", "mover-j", "pepe", "zhou", and "delta".
#'
#' Among those, \strong{gart and nam}, \strong{walter}, and \strong{mover-j} construct the confidence intervals for PPV and NPV by converting the confidence intervals for the ratio of two binomial proportions (\eqn{\phi=\frac{p_1}{p_0}}) where
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
#' \item The \strong{walter} method constructs the confidence interval for \eqn{\phi} based on \eqn{log(\phi)}. 0.5 is added to `x1`, `x0`, `n1`, and `n0`. Thus, no continuity correct should be applied additionally. This method has shown skewness concerns for
#' small ratios and sample sizes.
#'
#' \item The \strong{mover-j} method constructs the confidence interval for \eqn{\phi} from separate intervals for the individual group rates (i.e., \eqn{p_1} and \eqn{p_0}).
#' By applying the equal-tailed Jeffreys method (default 0.5 to each group), it may achieve a skewness-corrected interval for \eqn{\phi}. \cr
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
#'}
#' @return A list object contains the method and the estimates of sensitivity, specificity, PPV, NPV and their confidence intervals.
#' @importFrom ratesci scoreci
#' @import Rdpack kableExtra
#' @references
#' \enumerate{
#' \item \insertRef{gn_1988}{CIfinder}
#'
#' \item \insertRef{zhou_2007}{CIfinder}
#'
#' \item \insertRef{laud_2017}{CIfinder}
#'
#' \item \insertRef{pepe_2003}{CIfinder}
#'}
#' @importFrom ratesci moverci
#' @importFrom stats pnorm quantile
#' @export
#' @examples
#' ppv_npv_ci(60, 65, 113, 113, prevalence = 0.02)
ppv_npv_ci <- function (x1, n1, x0, n0, prevalence, method = "gart and nam", conf.level = 0.95, bias_correction = FALSE, continuity.correction = FALSE, ...) {

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
    out <- get_delta_ci(x10 = x10, x11 = x11, x01 = x01, x00 = x00, n1 = n1, n0 = n0, conf.level = conf.level, prev = prev)
  }

  return(out)
}

get_gn_ci <- function (x10, x11, x01, x00, n1, n0, conf.level, prev, bias_correction, continuity.correction) {

  if ((x11==0 | x00==0) & !continuity.correction) {
    stop("x1 or x0 must be greater than zero for Gart and Nam method, or consider to use continuity correction")
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

get_delta_ci <- function(x10, x11, x01, x00, n1, n0, prev, conf.level) {
  z <- qnorm(1-(1-conf.level)/2)
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
