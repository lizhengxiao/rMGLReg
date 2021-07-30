#' Density, distribution function, quantile function and random generation for the GLMGA distribution with parameters sigma, a and b.
#'
#' @param y
#' @param sigma
#' @param a
#' @param b
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#'
#' @return
#' @export
#'
#' @examples dLMGA(c(0.5,0.1), sigma = 2, a = 2, b = 3, log = FALSE) # density function at value 0.5 and 0.1
#'
dLMGA <- function(y, sigma, a, b, log = FALSE) {
  if( log == FALSE){
    exp(-0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b))
  } else if (log == TRUE){
    -0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  }
}



#' Density, distribution function, quantile function and random generation for the GLMGA distribution with parameters sigma, a and b.
#'
#' @param y
#' @param sigma
#' @param a
#' @param b
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#'
#' @return
#' @export
#'
#' @examples pLMGA(c(10,20), sigma = 2, a = 2, b = 3) # cdf at value 10 and 20.
pLMGA <- function(y, sigma, a, b) {
  z <- y^(-1/sigma)/(y^(-1/sigma) + 2*b)
  p <- 1- pbeta(z, shape1 = 0.5, shape2 = a)
  p
}


#' Density, distribution function, quantile function and random generation for the GLMGA distribution with parameters sigma, a and b.
#'
#' @param u
#' @param sigma
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples qLMGA(c(0.5,0.1), sigma = 2, a = 2, b = 3) # quantile function at level 50% and 10%
qLMGA <- function(u, sigma, a, b) {
  c <- (2*b)^(-sigma)
  #I <- pbeta(u, shape1 = 0.5, shape2 = a)
  Iinv <- qbeta(1- u, shape1 = 0.5, shape2 = a)
  c*(Iinv/(1 - Iinv))^(-sigma)
}



#' Density, distribution function, quantile function and random generation for the GLMGA distribution with parameters sigma, a and b.
#'
#' @param n
#' @param sigma
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples rLMGA(n = 10, sigma = 2, a = 2, b = 3) # simulate the data from GLMGA distribution
rLMGA <- function(n, sigma, a, b) {
  u <- runif(n, min = 0, max = 1)
  qLMGA <- Vectorize(qLMGA)
  r <- qLMGA(u, sigma = sigma, a = a, b = b)
  r
}


