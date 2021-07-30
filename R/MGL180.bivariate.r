#'
#' Survival bivarite MGL copula
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
#' dcMGL180.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
#' pcMGL180.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
dcMGL180.bivar <- function(u1, u2, pars) {
  dcMGL.bivar( 1 - u1, 1 - u2, pars)
}





#'
#' Survival bivarite MGL copula
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
#' dcMGL180.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
#' pcMGL180.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
pcMGL180.bivar <- function(u1, u2, pars) {
  u1 + u2 - 1 + pcMGL.bivar( 1 - u1, 1 - u2, pars)
}



#'
#' Survival bivarite MGL copula
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
#' dcMGL180.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
#' pcMGL180.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
#' Usim <- rcMGL180.bivar(5000, param = 3)
#' plot(Usim)
#'
rcMGL180.bivar <- function(n, param){
  Usim <- rcMGL.multi(n = n, param = param, d = 2)
  1 - Usim
}

