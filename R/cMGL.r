#' d-dimensional MGL copula
#'
#' @param u d-dimensional matrix
#' @param pars copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param d d-dimensions
#'
#' @returnDensity, Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGL.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 10, log = FALSE)
#' dcMGL.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = T)
#' pcMGL.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
#'
#' Usim <- rcMGL.multi(n = 1000, d = 2, param = 1)
#' Usim
dcMGL.multi <- function(u, pars, log = FALSE){
  # coding as a matrix
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

#' d-dimensional MGL copula
#'
#' @param u d-dimensional matrix
#' @param pars copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param d d-dimensions
#'
#' @returnDensity, Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGL.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 10, log = FALSE)
#' dcMGL.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = T)
#' pcMGL.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
#'
#' Usim <- rcMGL.multi(n = 1000, d = 2, param = 1)
#' Usim
pcMGL.multi <- function(u, pars) {
  dim <- ncol(u)
  a <- 1/pars[1]

  q <- (qbeta(1-u, shape1 = 0.5, shape2 = a)/(1 - qbeta(1-u, shape1 = 0.5, shape2 = a)))
  z <- 0
  for(i in 1:nrow(u)){
    fin <- function(theta){ # rely on a
      m <- pracma::erfc((q[i,]*theta)^0.5)
      k <- prod(m)
      out <- k*theta^(a-1)*exp(-theta)/gamma(a)
      out
    }
    z[i] <- as.numeric(pracma::integral(fun = fin,
                                        xmin = 0,
                                        xmax = Inf, method = 'Clenshaw'))
  }
  return(z) # rely on k
}

#' d-dimensional MGL copula
#'
#' @param u d-dimensional matrix
#' @param pars copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param d d-dimensions
#'
#' @returnDensity, Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGL.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 10, log = FALSE)
#' dcMGL.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = T)
#' pcMGL.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
#'
#' Usim <- rcMGL.multi(n = 1000, d = 2, param = 1)
#' Usim
rcMGL.multi <- function(n, d, param){
  a <- 1/param
  dseq <- seq(1:(d-1))
  anew <- c(a, a + dseq/2)
  Iinv_trans <- function(u, a) {
    qbeta(1 - u, shape1 = 0.5, shape2 = a)/(1 - qbeta(1 - u, shape1 = 0.5, shape2 = a))
  }
  umat <- matrix(data = runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  mmat <- qmat <- zmat <- Usim <- matrix(data = 0, nrow = n, ncol = d)

  for(j in 1:ncol(umat)){
    qmat[,j] <- Iinv_trans(umat[,j], a = anew[j]) # column j denots the qj
  }

  mmat[,1] <- qmat[, 1]
  mmat[,2] <- (1 + mmat[,1])*qmat[, 2]
  if(d >2){
    for(j in 3:ncol(umat)){
      tempmat <- mmat[, (1):(ncol(umat)-1)]
      if(n == 1 ){
        msum <- sum(tempmat)
      } else {
        msum <- apply(tempmat, MARGIN = 1, sum)
      }
      mmat[,j] <- (1 + msum)*qmat[, j]
    }
  }

  for(j in 1:ncol(zmat)){
    k <- mmat[,j]/(mmat[,j] + 1)
    Usim[, j] <- 1 - pbeta(k, shape1 = 0.5, shape2 = a)
  }

  return(Usim)
}

#'
#' Bivarite MGL copula
#'
#' @param u_1,u_2
#' @param pars copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @returnDensity, Density, distribution function, and h function for the bivariate MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGL.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
#' pcMGL.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
dcMGL.bivar <- function(u1, u2, pars){
  dim <- 2
  a <- 1/pars
  q1 <- qbeta(1-u1, shape1 = 0.5, shape2 = a)/(1 - qbeta(1-u1, shape1 = 0.5, shape2 = a))
  q2 <- qbeta(1-u2, shape1 = 0.5, shape2 = a)/(1 - qbeta(1-u2, shape1 = 0.5, shape2 = a))
  # if(u1 == 0|u2 == 0|q1 == Inf|q2 == Inf) dc <- 0
  if(q1 == Inf|q2 == Inf) {
    dc <- 0
  } else {
    q <- cbind(q1, q2)
    #dc <- gamma(a)^(dim - 1)*gamma(a + dim/2)/(gamma(a + 0.5)^dim)*(q1 + 1)^(a + 0.5)*(q2 + 1)^(a + 0.5)/((q1 + q2) + 1)^(a + dim/2)
    dc <- gamma(a)^(dim - 1)/(gamma(a + 0.5)^(dim - 1))/(gamma(a + 0.5))*gamma(a + dim/2)*(q1 + 1)^(a + 0.5)*(q2 + 1)^(a + 0.5)/((q1 + q2) + 1)^(a + dim/2)
  }
  return(dc)
}
dcMGL.bivar <- Vectorize(dcMGL.bivar)



#'
#' Bivarite MGL copula
#'
#' @param u_1,u_2
#' @param pars copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @returnDensity, Density, distribution function, and h function for the bivariate MGL copula with copula parameter delta.
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
#' h_2(u_1|u_2,θ) := P(U_1 ≤ u_1 | U_2 = u_2) := \partial C(u_1,u_2) / \partial u_2,
#'  where (U_1, U_2) \sim C, and C is a bivariate copula distribution function with parameter(s) θ. For more details see Aas et al. (2009).
#' @Value BiCopHfunc returns a list with
#'
#'
#'
#' @export
#'
#' @References Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence. Insurance: Mathematics and Economics 44 (2), 182-198.
#' @examples
#' hMGL.bivar(u1 = c(0.1, 0.001, 0.3), u2 = c(0, 0.9999, 0.88), pars = 2)
hMGL.bivar <- function(u1, u2, pars) {
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
