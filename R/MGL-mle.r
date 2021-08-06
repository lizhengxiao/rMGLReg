
#' @title Fitting bivariate MGL copula models
#' @description MGL.reg is used to fit bivariate MGL copula regression models.
#' @param U two-dimenstional matrix with values in \eqn{\left\[0,1\right\]}.
#' @param copula copula 'MGL', 'MGL180', "MGL-EV", "MGL-EV180", "MGB2", "Normal" , "Student-t"
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param initpar Initial values for the parameters to be optimized over.
#' @param ... additional arguments to be passed to f.
#' @importFrom stats nlm
#' @return A list containing the following components
#' @export
#'
MGL.mle <- function(U, copula = c(
                      "MGL", "MGL180", "MGL-EV",
                      "MGL-EV180",
                      "Gumbel",
                      "Normal", "MGB2", "t"
                    ),
                    hessian = TRUE,
                    initpar, ...) {
  dnormcop <- function(U, param) {
    as.numeric(fCopulae::dellipticalCopula(U, rho = param[1], type = "norm"))
  } # normal copula

  dtcop <- function(U, param) {
    as.numeric(fCopulae::dellipticalCopula(U,
      rho = param[1], type = "t",
      param = param[2]
    ))
  } # t copula
  dgumcop <- function(U, param) {
    as.numeric(fCopulae::devCopula(U, type = "gumbel", param = param[1]))
  } # Bivariate Extreme


  dMGL <- function(U, param) {
    as.numeric(dcMGL.bivar(u1 = U[, 1], u2 = U[, 2], pars = param[1]))
  }

  dMGL180 <- function(U, param) {
    as.numeric(dcMGL180.bivar(u1 = U[, 1], u2 = U[, 2], pars = param[1]))
  }

  dMGB2 <- function(U, param) {
    as.numeric(dcMGB2.bivar(u1 = U[, 1], u2 = U[, 2], pars1 = param[1], pars2 = param[2], pars3 = param[3]))
  }

  dMGLEV <- function(U, param) {
    as.numeric(dcMGLEV.bivar(u1 = U[, 1], u2 = U[, 2], param = param[1])) # Bivariate Extreme
  }

  dMGLEV180 <- function(U, param) {
    as.numeric(dcMGLEV180.bivar(u1 = U[, 1], u2 = U[, 2], param = param[1])) # Bivariate Extreme
  }

  argnorm <- list(length = 1, lower = 0, upper = 1, name = "Gaussian")
  argt <- list(length = 2, lower = c(0, 0), upper = c(1, 100), name = "Student")
  arggum <- list(length = 1, lower = 1, upper = 50, name = "Gumbel")
  argMGB2 <- list(
    length = 3, lower = c(0, 0, 0), upper = c(50, 50, 50),
    name = "MGB2"
  )
  argMG180 <- list(length = 1, lower = 0, upper = 10, name = "MGL180")
  argMGLEV180 <- list(length = 1, lower = 0, upper = 10, name = "MGL-EV180")
  argMG <- list(length = 1, lower = 0, upper = 10, name = "MGL")
  argMGLEV <- list(length = 1, lower = 0, upper = 10, name = "MGL-EV")
  if (copula == "MGL") {
    dcop <- dMGL
    arg.cop <- argMG
  } else if (copula == "MGL180") {
    dcop <- dMGL180
    arg.cop <- argMG180
  } else if (copula == "MGL-EV") {
    dcop <- dMGLEV
    arg.cop <- argMGLEV
  } else if (copula == "MGL-EV180") {
    dcop <- dMGLEV180
    arg.cop <- argMGLEV180
  } else if (copula == "Gumbel") {
    dcop <- dgumcop
    arg.cop <- arggum
  } else if (copula == "Normal") {
    dcop <- dnormcop
    arg.cop <- argnorm
  } else if (copula == "t") {
    dcop <- dtcop
    arg.cop <- argt
  } else if (copula == "MGB2") {
    dcop <- dMGB2
    arg.cop <- argMGB2
  }

  # loglike.copula <- function(U, initpar, ...){
  copLogL <- function(x) {
    if (all(arg.cop$lower < x && arg.cop$upper > x)) {
      logL <- log(dcop(U, param = x))
      # res <- -sum((logL))
      remove.naninf <- function(x) x[!is.nan(x) & is.finite(x)]
      res <- -sum(remove.naninf(logL))
    } else {
      res <- 10000000000000
    }
    return(res)
  }
  # }

  resopt <- nlm(
    f = copLogL,
    # U = U,
    p = initpar,
    # dcop = dcop, arg.cop = arg.cop,
    hessian = hessian
  )


  # resopt
  list(
    loglike = -resopt$minimum,
    copula = list(name = arg.cop$name),
    estimates = resopt$estimate,
    se = sqrt(diag(solve(resopt$hessian))),
    AIC = 2 * length(resopt$estimate) + 2 * resopt$minimum,
    BIC = log(nrow(U)) * length(resopt$estimate) + 2 * resopt$minimum
  )
}
