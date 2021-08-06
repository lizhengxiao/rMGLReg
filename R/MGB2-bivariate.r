


#' @name MGB2
#' @rdname  MGB2
#' @title Bivarite MGB2 copula
#' @description Density, distribution function, and random generation for the 2-dimensional MGB2.
#' @param u1,u2 numeric vectors of equal length with values in \eqn{\left\[0,1\right\]}.
#' @param pars1,pars2,pars3 copula parameters, denoted by p_1, p_2, q
#' @importFrom stats qbeta
#' @importFrom stats pbeta
NULL

#' @rdname  MGB2
#' @export
#' @examples
#' dcMGB2.bivar(u1 = c(0.2, 0.5), u2 = c(0.9, 0.2), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)
dcMGB2.bivar <- function(u1, u2, pars1, pars2, pars3) {
  dim <- 2
  p1 <- pars1
  p2 <- pars2
  q <- pars3

  x1 <- qbeta(u1, shape1 = p1, shape2 = q) / (1 - qbeta(u1, shape1 = p1, shape2 = q))
  x2 <- qbeta(u2, shape1 = p2, shape2 = q) / (1 - qbeta(u2, shape1 = p2, shape2 = q))

  if (u1 == 0 | u2 == 0) {
    dc <- 0
  } else {
    dc <- gamma(q)^(dim - 1) * gamma(p1 + p2 + q) / (gamma(p1 + q) * gamma(p2 + q)) * (1 + x1)^(p1 + q) * (1 + x2)^(p2 + q) / (1 + x1 + x2)^(p1 + p2 + q)
  }
  dc
}
dcMGB2.bivar <- Vectorize(dcMGB2.bivar)



#' @rdname MGB2
#' @export
#' @examples
#' pcMGB2.bivar(u1 = c(0.5, 0.1), u2 = c(0.9, 0.1), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)
pcMGB2.bivar <- function(u1, u2, pars1, pars2, pars3) {
  p1 <- pars1
  p2 <- pars2
  q <- pars3

  x1 <- qbeta(u1, shape1 = p1, shape2 = q) / (1 - qbeta(u1, shape1 = p1, shape2 = q))
  x2 <- qbeta(u2, shape1 = p2, shape2 = q) / (1 - qbeta(u2, shape1 = p2, shape2 = q))
  fin <- function(theta) { # rely on a
    m1 <- x1 / theta
    m2 <- x2 / theta
    z1 <- pgamma(m1, shape = p1, scale = 1) * pgamma(m2, shape = p2, scale = 1)
    out <- z1 * theta^(-(q + 1)) * exp(-1 / theta) / gamma(q)
    out
  }
  if (u1 == 0 | u2 == 0) {
    z <- 0
  } else {
    z <- as.numeric(pracma::integral(fun = fin, xmin = 0, xmax = Inf, method = "Clenshaw"))
    # z <- as.numeric(integrate(f = fin, lower = 0, upper = Inf, rel.tol = 1e-4)$value)
  }
  return(z) # rely on k
}
pcMGB2.bivar <- Vectorize(pcMGB2.bivar)


#' @rdname MGB2
#' @export
hcMGB2.bivar <- function(u1, u2, pars1, pars2, pars3) {
  p1 <- pars1
  p2 <- pars2
  q <- pars3

  x1 <- qbeta(u1, shape1 = p1, shape2 = q) / (1 - qbeta(u1, shape1 = p1, shape2 = q))
  x2 <- qbeta(u2, shape1 = p2, shape2 = q) / (1 - qbeta(u2, shape1 = p2, shape2 = q))

  index1 <- u1 == 1 & u2 != 1
  index2 <- u1 != 1 & u2 == 1
  index3 <- u1 == 1 & u2 == 1
  z1 <- x2 / (x1 + x2 + 1)
  z2 <- x1 / (x2 + x1 + 1)
  z1[index1] <- 0
  z2[index1] <- 1
  z1[index2] <- 1
  z2[index2] <- 0
  z1[index3] <- 1
  z2[index3] <- 1
  hfunc1 <- pbeta(z1, shape1 = p2, shape2 = p1 + q)
  hfunc2 <- pbeta(z2, shape1 = p1, shape2 = p2 + q)
  list(hfunc1 = hfunc1, hfunc2 = hfunc2)
}
