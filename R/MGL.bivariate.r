#'
#' Bivarite MGL copula
#'
#' @param u_1,u_2 numeric vectors of equal length with values in [0,1].
#' @param pars copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @return Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGL.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
#' pcMGL.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
dcGLMGA.multi <- function(u, pars, log = FALSE){
  # 1. coding as a matrix
  dim <- ncol(u)
  a <- 1/pars[1]
  q <- (qbeta(1-u, shape1 = 0.5, shape2 = a)/(1 - qbeta(1-u, shape1 = 0.5, shape2 = a)))
  logdc <- 0
  for(i in 1:nrow(u)){
    logdc[i] <- (dim - 1)*lgamma(a) + lgamma(a + dim/2) - dim*lgamma(a + 0.5) + (a + 0.5)*sum(log(q[i,] + 1)) - (a + dim/2)*log(sum(q[i,]) + 1)
  }
  dc <- exp(logdc)
  dc[which(q == Inf)] <- 0
  if(log == TRUE) {logdc} else {dc}
}

#'
#' Bivarite MGL copula
#'
#' @param u_1,u_2 numeric vectors of equal length with values in [0,1].
#' @param pars copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @return Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGL.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
#' pcMGL.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
pcMGL.bivar <- function(u1, u2, pars) {
  a <- 1/pars
  q1 <- qbeta(1-u1, shape1 = 0.5, shape2 = a)/(1 - qbeta(1-u1, shape1 = 0.5, shape2 = a))
  q2 <- qbeta(1-u2, shape1 = 0.5, shape2 = a)/(1 - qbeta(1-u2, shape1 = 0.5, shape2 = a))
  fin <- function(theta){ # rely on a
    m1 <- pracma::erfc((q1*theta)^0.5)
    m2 <- pracma::erfc((q2*theta)^0.5)
    out <- m1*m2*theta^(a-1)*exp(-theta)/gamma(a)
    out
  }
  if(u1==0|u2==0){
    z <- 0
  } else if (u1==1&u2!=0){
    z <- u2
  } else if (u1!=0&u2==1){
    z <- u1
  } else {
    z <- as.numeric(pracma::integral(fun = fin, xmin = 0, xmax = Inf, method = 'Clenshaw'))
  }
  return(z) # rely on k
}
pcMGL.bivar <- Vectorize(pcMGL.bivar)




#' Conditional Distribution Function of a Bivariate Copula
#'
#' @param u1 numeric vectors of equal length with values in [0,1].
#' @param u2 numeric vectors of equal length with values in [0,1].
#' @param pars numeric; single number or vector of size length(u1); copula parameter > 0.
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
#' @References Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence. Insurance: Mathematics and Economics 44 (2), 182-198.
#' @examples
#' hcMGL.bivar(u1 = c(0.1, 0.001, 0.3), u2 = c(0, 0.9999, 0.88), pars = 2)
hcMGL.bivar <- function(u1, u2, pars) {
  a <- 1/pars
  # another method
  index1 <- u1 == 0 & u2 != 0
  index2 <- u1 != 0& u2 == 0
  index3 <- u1 == 0 & u2 == 0
  q1 <- qbeta(1 - u1, shape1 = 0.5, shape2 = a)/(1 - qbeta(1 - u1, shape1 = 0.5, shape2 = a))
  q2 <- qbeta(1 - u2, shape1 = 0.5, shape2 = a)/(1 - qbeta(1 - u2, shape1 = 0.5, shape2 = a))
  z1 <- q2/(q1 + q2 + 1)
  z2 <- q1/(q1 + q2 + 1)
  z1[index1] <- 0.0000000001; z2[index1] <- 1
  z1[index2] <- 1; z2[index2] <- 0.0000000001
  z1[index3] <- 0.5 ; z2[index3] <- 0.5

  hfunc1 <- 1 - pbeta(z1 ,shape1 = 0.5, shape2 = a + 0.5) # Pr(U_2 <= u_2 | U_1 = u_1) = \partial C(u_1,u_2) / \partial u_1,
  hfunc2 <- 1 - pbeta(z2 ,shape1 = 0.5, shape2 = a + 0.5) # C(u1|u2)
  list(hfunc1 = hfunc1, hfunc2 = hfunc2)

}
