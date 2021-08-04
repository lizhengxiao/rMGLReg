


#'
#' Bivarite MGL copula
#'
#' @description Density (dCopula), distribution function (pCopula), and random generation (rCopula) for a copula object.
#' @param u_1,u_2 numeric vectors of equal length with values in [0,1].
#' @param pars1,pars2,pars3 copula parameters, denoted by p_1, p_2, q
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @return Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGB2.bivar(u1 = c(0, 0.5), u2 = c(0.9, 0), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)
#' pcMGB2.bivar(u1 = c(0.5, 0.1), u2 = c(0.9, 0.1), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)
#'
dcMGB2.bivar <- function(u1, u2, pars1, pars2, pars3){
  dim <- 2
  p1 <- pars1
  p2 <- pars2
  q <- pars3

  x1 <- qbeta(u1, shape1 = p1, shape2 = q)/( 1 - qbeta(u1, shape1 = p1, shape2 = q))
  x2 <- qbeta(u2, shape1 = p2, shape2 = q)/( 1 - qbeta(u2, shape1 = p2, shape2 = q))

  if(u1 == 0|u2 == 0) dc <- 0
  else {
    dc <- gamma(q)^(dim - 1)*gamma(p1 + p2 + q)/(gamma(p1 + q)*gamma(p2 + q))*(1 + x1)^(p1 + q)*(1 + x2)^(p2 + q)/(1 + x1 + x2)^(p1 + p2 + q)
  }
  dc
}
dcMGB2.bivar <- Vectorize(dcMGB2.bivar)


#'
#' Bivarite MGL copula
#' @description Density (dCopula), distribution function (pCopula), and random generation (rCopula) for a copula object.
#' @param u_1,u_2 numeric vectors of equal length with values in [0,1].
#' @param pars1,pars2,pars3 copula parameters, denoted by p_1, p_2, q
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @return Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGB2.bivar(u1 = c(0, 0.5), u2 = c(0.9, 0), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)
#' pcMGB2.bivar(u1 = c(0.5, 0.1), u2 = c(0.9, 0.1), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)
#'
pcMGB2.bivar <- function(u1, u2, pars1, pars2, pars3){
  p1 <- pars1
  p2 <- pars2
  q <- pars3

  x1 <- qbeta(u1, shape1 = p1, shape2 = q)/( 1 - qbeta(u1, shape1 = p1, shape2 = q))
  x2 <- qbeta(u2, shape1 = p2, shape2 = q)/( 1 - qbeta(u2, shape1 = p2, shape2 = q))
  fin <- function(theta){ # rely on a
    m1 <- x1/theta
    m2 <- x2/theta
    z1 <- pgamma(m1, shape = p1, scale = 1)*pgamma(m2, shape = p2, scale = 1)
    out <- z1*theta^(-(q + 1))*exp(-1/theta)/gamma(q)
    out
  }
  if(u1==0|u2==0){
    z <- 0
  } else {
    z <- as.numeric(pracma::integral(fun = fin, xmin = 0, xmax = Inf, method = 'Clenshaw'))
    #z <- as.numeric(integrate(f = fin, lower = 0, upper = Inf, rel.tol = 1e-4)$value)
  }
  return(z) # rely on k
}
pcMGB2.bivar <- Vectorize(pcMGB2.bivar)


#' Conditional Distribution Function of a Bivariate Copula
#'
#' @param u1,u2 numeric vectors of equal length with values in [0,1].
#' @param pars1,pars2,pars3 copula parameters, denoted by p_1, p_2, q
#'
#' @details
#' The h-function is defined as the conditional distribution function of a bivariate copula, i.e.,
#' h_1(u_2|u_1,θ) := P(U_2 ≤ u_2 | U_1 = u_1) = \partial C(u_1,u_2) / \partial u_1,
#'
#' h_2(u_1|u_2,θ) := P(U_1 ≤ u_1 | U_2 = u_2) := \partial C(u_1,u_2) / \partial u_2,
#'
#'where (U_1, U_2) \sim C, and C is a bivariate copula distribution function with parameter(s) θ. For more details see Aas et al. (2009).
#' @Value BiCopHfunc returns a list with
#'
#'
#'
#' @export
#'
#' @References \
#' Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence. Insurance: Mathematics and Economics 44 (2), 182-198.
#' @examples
#' hMGB2.bivar(u1 = c(0.3, 0.4), u2 = c(0.4, 0.5), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)
#'
hcMGB2.bivar <- function(u1, u2, pars1, pars2, pars3){
  p1 <- pars1
  p2 <- pars2
  q <- pars3

  x1 <- qbeta(u1, shape1 = p1, shape2 = q)/( 1 - qbeta(u1, shape1 = p1, shape2 = q))
  x2 <- qbeta(u2, shape1 = p2, shape2 = q)/( 1 - qbeta(u2, shape1 = p2, shape2 = q))

  index1 <- u1 == 1 & u2 != 1
  index2 <- u1 != 1& u2 == 1
  index3 <- u1 == 1 & u2 == 1
  z1 <- x2/(x1 + x2 + 1)
  z2 <- x1/(x2 + x1 + 1)
  z1[index1] <- 0; z2[index1] <- 1
  z1[index2] <- 1; z2[index2] <- 0
  z1[index3] <- 1; z2[index3] <- 1
  hfunc1 <- pbeta(z1 ,shape1 = p2, shape2 = p1 + q)
  hfunc2 <- pbeta(z2 ,shape1 = p1, shape2 = p2 + q)
  list(hfunc1 = hfunc1, hfunc2 = hfunc2)
}

