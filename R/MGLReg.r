



#' Fitting bivariate MGL copula regression models
#'
#' @description MGL.reg is used to fit bivariate MGL copula regression models.
#' @param U two-dimenstional matrix with values in [0,1].
#' @param X design matrix
#' @param copula 'MGL', 'MGL180', "MGL-EV", "MGL-EV180"
#' @param initpar Initial values for the parameters to be optimized over.
#' @param ... 	Further arguments to be passed to fn and gr in optiom.
#'
#' @return For optim, a list with components:
#' @export
#'
#' @examples
MGL.reg <- function(U, X, copula = c('MGL', 'MGL180', "MGL-EV", "MGL-EV180"), initpar, ...){

  dcMGL.reg <- function(U, param) {
    dim <- length(U)
    a <- 1/param[1]
    q <- qbeta(1 - U, shape1 = 0.5, shape2 = a)/(1 - qbeta(1 - U, shape1 = 0.5, shape2 = a))
    logdc <- (dim - 1)*lgamma(a) + lgamma(a + dim/2) - dim*lgamma(a + 0.5) + (a + 0.5)*sum(log(q + 1)) - (a + dim/2)*log(sum(q) + 1)
    out <- exp(logdc)
    return(out)
  }

  dcMGL180.reg <- function(U, param) {
    dcMGL.reg(1 - U, param = param)
  }

  dcMGLEV180.reg <- function(U, param) {
    u1 <- U[1];
    u2 <- U[2];
    as.numeric(dcMGLEV180.bivar(u1, u2, param = param[1]))
  }

  dcMGLEV.reg <- function(U, param) {
    u1 <- U[1];
    u2 <- U[2];
    as.numeric(dcMGLEV.bivar(u1, u2, param = param[1]))
  }


  if(copula == 'MGL'){
    dcop <- dcMGL.reg
  } else if(copula == 'MGL180'){
    dcop = dcMGL180.reg
  } else if(copula == "MGL-EV"){
    dcop = dcMGLEV.reg
  } else if(copula == "MGL-EV180"){
    dcop  <- dcMGLEV180.reg
  }

  copLogL <- function(pars, X) {
    ll <- 0
    delta <- exp(X%*%pars)

    for (i in 1:nrow(X)){
      ll[i] <- dcop(U[i,], param = as.vector(delta[i]))
    }
    res <- sum((log(ll)))
    return(res)
  }

  resopt <- optim(par = initpar,
                  fn = copLogL,
                  X = X,
                  hessian = T, ...)

  resopt
}
