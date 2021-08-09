#' @title the squared fit error
#'
#' @param XY the 2-dim loss data
#' @param copula  copula = c("MGL", "MGL180", "MGL-EV", "MGL-EV180", "Gumbel", "Normal", "MGB2", "t"
#' @param lowerLimit The lower limit of integration, a vector for hypercubes.
#' @param upperLimit The upper limit of integration, a vector for hypercubes.
#' @param pars copula parameters
#' @param ... All other arguments passed to the function f.
#' @importFrom copula pellipticalCopula
#' @importFrom copula pevCopula
#' @importFrom copula C.n
#' @importFrom cubature cuhre
#' @return vector of length nComp; the integral of integrand over the hypercube
#' @export
#'
error.copula <- function(XY, copula = c(
  "MGL", "MGL180", "MGL-EV",
  "MGL-EV180",
  "Gumbel",
  "Normal", "MGB2", "t"
),
lowerLimit,
upperLimit,
pars, ... ){

  pnormcop <- function(U, param)
    as.numeric(pellipticalCopula(U, rho = param[1], type = "norm")) # normal copula
  ptcop <- function(U, param)
    as.numeric(pellipticalCopula(U, rho = param[1], type = "t",
                                 param = param[2]))  # t copula
  pgumcop <- function(U, param)
    as.numeric(pevCopula(U, type = "gumbel", param = param[1])) # Bivariate Extreme Value Copulae - gumbel copula
  pMGLcop <- function(U, param){
    as.numeric(pcMGL.bivar(u1 = U[1], u2 = U[2], pars = param[1]))
  }

  pMGL180cop <- function(U, param){
    as.numeric(U[1] + U[2] - 1 + pcMGL.bivar(u1 = 1 - U[1], u2 = 1 - U[2], pars = param[1]))
  }

  pMGB2cop <- function(U, param){
    as.numeric(pcMGB2.bivar(u1 = U[1], u2 = U[2], pars1 = param[1], pars2 = param[2], pars3 = param[3]))
  }

  pMGLEVcop <- function(U, param){
    as.numeric(pcMGLEV.bivar(u1 = U[1], u2 = U[2], param = param[1]))
  }

  pMGLEV180cop <- function(U, param){
    as.numeric(pcMGLEV180.bivar(u1 = U[1], u2 = U[2], param = param[1]))
  }

  if (copula == "MGL") {
    pcop <- pMGLcop
  } else if (copula == "MGL180") {
    pcop <- pMGL180cop
  } else if (copula == "MGL-EV") {
    pcop <- pMGLEVcop
  } else if (copula == "MGL-EV180") {
    pcop <- pMGLEV180cop
  } else if (copula == "Gumbel") {
    pcop <- pgumcop
  } else if (copula == "Normal") {
    pcop <- pnormcop
  } else if (copula == "t") {
    pcop <- ptcop
  } else if (copula == "MGB2") {
    pcop <- pMGB2cop
  }


  error_c <- function(U){
    u1 <- U[1]
    u2 <- U[2]
    empc <- C.n(u = U, X = XY)
    fitc <- do.call(pcop, args = list(U = U, param = pars))
    out <- (fitc - empc)^2
    out
  }

  error <- cubature::cuhre(f = error_c,
                           lowerLimit = lowerLimit,
                           upperLimit = upperLimit,
                           ...)



  error$integral

}
