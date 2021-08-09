

#' @name MGB2
#' @rdname  MGB2
#' @aliases dMGB2  pMGB2  hMGB2
#' @title Bivarite MGB2 copula
#' @description Density, distribution function, and h-functions for the bivariate MGB2 copula proposed in Yang et al.,(2011).
#' @param u1,u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param pars1,pars2,pars3 copula parameters, denoted by \eqn{p_1, p_2, q}.
#' @importFrom stats qbeta
#' @importFrom stats pbeta
#' @references Xipei Yang, Edward W Frees, and Zhengjun Zhang. A generalized beta copula with applications in
#' modeling multivariate long-tailed data. Insurance: Mathematics and Economics, 49(2):265-284, 2011.
#' @details
#' The MGB2 copula with parameters (p1, p2, q) has joint density
#'		\deqn{c(u_1,u_2;p_1,p_2,q)=\frac{\Gamma(q)\Gamma(\sum_{i=1}^2 p_i + q)}{\prod_{i=1}^2\Gamma(p_i+q)}\frac{\prod_{i=1}^2 (1 + x(u_i;p_i,q))^{p_i+q}}{(1 + \sum_{i=1}^{2}x(u_i;p_i,q))^{\sum_{i=1}^{2}p_i+q}},}
#' for \eqn{p_1, p_2>0, q>0}.
#' (Here Gamma(a) is the function implemented by R's \code{\link[base]{gamma}} and defined in its help.
#'
#'
#' The joint cdf of the MGB2 copula is
#'
#' \deqn{C(u_1,u_2;p_1,p_2,q)=\int_{0}^{+\infty}\prod_{i=1}^{2}G_p(\frac{I_{p_i,q}^{-1}(u_i)}{(1-I_{p_i,q}^{-1}(u_i))\theta})\times \frac{\theta^{-(q+1)}e^{-1/\theta}}{\Gamma(q)}d\theta, }
#' where \eqn{I_{m,n}^{-1}()} denotes the inverse of the beta cumulative distribution function (or regularized incomplete beta function)
#' with parameters shape1 = m and shape2 = n
#' implemented by R's \code{\link[stats]{qbeta}}.
#'
#' The h-function is defined as the conditional distribution function of a bivariate copula, i.e.,
#' \deqn{h_1(u_2|u_1,p_1,p_2,q) := P(U_2 \leq u_2 | U_1 = u_1) = \partial C(u_1,u_2) / \partial u_1,}
#'
#' \deqn{h_2(u_1|u_2,p_1,p_2,q) := P(U_1 \leq u_1 | U_2 = u_2) := \partial C(u_1,u_2) / \partial u_2,}
#'
#' where \eqn{(U_1, U_2) \sim C}, and \eqn{C} is a bivariate copula distribution function with parameter(s) \eqn{p_1,p_2,q}.
#'
#'
#' @return
#' \code{dcMGB2.bivar} gives the density
#'
#' \code{pcMGB2.bivar} gives the distribution function
#'
#' \code{hcMGB2.bivar} gives the h-functions (hfunc1, hfunc2).
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
#' @examples
#' hcMGB2.bivar(u1 = c(0.5, 0.1), u2 = c(0.9, 0.1), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)$hfunc1
#' hcMGB2.bivar(u1 = c(0.5, 0.1), u2 = c(0.9, 0.1), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)$hfunc2
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
