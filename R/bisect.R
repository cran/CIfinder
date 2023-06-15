bisect <- function (ftn, precis = 6, max.iter = 100, uplow = "low") {
  tiny <- (10^-(precis))/2
  hi <- 1
  lo <- -1
  dp <- 2
  niter <- 1
  while (niter <= max.iter && any(dp > tiny | is.na(hi))) {
    dp <- 0.5 * dp
    mid <- pmax(-1, pmin(1, round((hi + lo)/2, 10)))
    scor <- ftn(round(tan(pi * (mid + 1)/4), 10))
    check <- (scor <= 0) | is.na(scor)
    hi[check] <- mid[check]
    lo[!check] <- mid[!check]
    niter <- niter + 1
  }
  if (uplow == "low") {
    best <- lo
  }
  else {
    best <- hi
  }
  return(tan((best + 1) * pi/4))
}
