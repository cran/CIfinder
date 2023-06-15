scorephi <- function(phi, x1, n1, x0, n0, skew = TRUE, cc = FALSE, conf.level = 0.95, alternative = "two.sided") {

  if (alternative == "two.sided") {
    alpha <- 1- conf.level
  } else {
    alpha <- (1 - conf.level) * 2
  }

  p1hat <- x1 / n1
  p0hat <- x0 / n0
  x <- x1 + x0
  N <- n1 + n0
  bias <- 0

  lambda <- ifelse(cc, N/(N-1), 1)

  #see page 346 LAUD, 2017
  Sphi <- p1hat - p0hat * phi
  # Ap0^2+bp0+c=0
  A <- N * phi
  B <- (-(n1 * phi + x1 + n0 + x0 * phi))
  C <- x
  num <- (-B - sqrt(pmax(0, B^2 - 4 * A * C)))
  p0t <- ifelse(A == 0, -C / B,
                ifelse((num == 0 | C == 0), 0, num / (2 * A)))
  p1t <- p0t * phi

  V <- pmax(0, (p1t * (1 - p1t) / n1 + (phi^2) * p0t *
                  (1 - p0t) / n0)*lambda)
  mu3 <- (p1t * (1 - p1t) * (1 - 2 * p1t) / (n1^2) -
            (phi^3) * p0t * (1 - p0t) * (1 - 2 * p0t) / (n0^2))
  V[is.na(V)] <- Inf

  ##apply to equation(3), page 336, LAUD 2017
  scterm <- mu3 / (6 * V^(3 / 2))
  scterm[mu3 == 0] <- 0
  score1 <- Sphi / sqrt(V)
  score1[Sphi == 0] <- 0
  A <- scterm
  B <- 1
  C <- -(score1 + scterm)
  num <- (-B + sqrt(pmax(0, B^2 - 4 * A * C)))
  dsct <- B^2 - 4 * A * C
  score <- ifelse((skew == FALSE | scterm == 0), score1,
                  num / (2 * A)
  )
  if (skew == TRUE) {
    qtnorm <- stats::qnorm(1 - alpha / 2)
    scoresimp <- score1 - (qtnorm^2 - 1) * scterm
    score[!is.na(dsct) & dsct < 0] <- scoresimp[!is.na(dsct) &
                                                  dsct < 0]
  }

  pval <- stats::pnorm(score)

  outlist <- list(
    score = score, p1t = p1t, Sphi = Sphi,
    num = num, V = V, p0t = p0t, mu3 = mu3, pval = pval,
    dsct = dsct
  )

  return(outlist)
}
