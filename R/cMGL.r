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
#' dcGLMGA.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 10, log = FALSE)
#' dcGLMGA.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = T)
#' pcGLMGA.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
#'
#' Usim <- rcGLMGA.multi(n = 1000, d = 2, param = 1)
#' Usim
dcGLMGA.multi <- function(u, pars, log = FALSE){
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
#' dcGLMGA.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 10, log = FALSE)
#' dcGLMGA.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = T)
#' pcGLMGA.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
#'
#' Usim <- rcGLMGA.multi(n = 1000, d = 2, param = 1)
#' Usim
pcGLMGA.multi <- function(u, pars) {
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
#' dcGLMGA.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 10, log = FALSE)
#' dcGLMGA.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = T)
#' pcGLMGA.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
#' Usim <- rcGLMGA.multi(n = 1000, d = 2, param = 1)
#' Usim
rcGLMGA.multi <- function(n, d, param){
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



