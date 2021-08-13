



#' Fitting bivariate mixed MGL and MGL-EV copula regression models
#'
#' @description \code{MGL.reg.mixed} is used to fit bivariate MGL and MGL-EV copula regression models for continuous and semi-continuous variables.
#' @param obs two-dimenstional matrix for observations
#' @param U two-dimenstional matrix for pseudo copula data with values in \eqn{[0,1]} for (F(y1), F(y2)).
#' @param U_ two-dimenstional matrix for pseudo copula data for the data (F(y1), F(y2-1)).
#' @param X design matrix.
#' @param copula 'MGL', 'MGL180', "MGL-EV", "MGL-EV180"
#' @param umin Threshold
#' @param f values of the density function for marginal distribution.
#' @param initpar Initial values for the parameters to be optimized over.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param ... additional arguments, see \code{\link[stats]{optim}} for more details.
#' @importFrom stats qbeta
#' @importFrom stats optim
#' @return A list containing the following components:
#' * loglike: the value of the estimated maximum of the loglikelihood function.
#' * copula: the name of the fitted copula. "MGL180" and "MGL-EV180" denote the survival MGL and MGL-EV copula respectively.
#' * estimates: the point at which the maximum value of the loglikelihood is obtained.
#' * se: the standard errors of the estimators.
#' * AIC, BIC: the goodness fit of the regression models.
#' * hessian: the hessian at the estimated maximum of the loglikelihood (if requested).
#' @details
#' * Y1: continuous data.
#' * Y2: semi-continuous data where Y2>umin is continuous and Y2<=umin is discrete.
#' * copula: "MGL180" and "MGLEV180" denote the survival MGL and survival MGL-EV copula respectively.
#' * For "Gumbel" regression model, the copula parameter \deqn{\delta_i = \exp(X\beta) + 1.}
#' * For "MGL", "MGL180", "MGL-EV", "MGL-EV180" regression model, the copula parameter \deqn{\delta_i = \exp(X\beta),} where \eqn{\beta} is the vector of coefficients to be estimated in the copula regression.
#'
#'
#'
#' @export
#'
#'
MGL.reg.mixed <- function(obs, U, U_, f, X, copula = c(
                            "MGL", "MGL180", "MGL-EV",
                            "MGL-EV180",
                            "Gumbel"
                          ),
                          umin = 0,
                          hessian = TRUE, initpar, ...) {
  dcMGL.reg <- function(U, param) {
    dim <- length(U)
    a <- 1 / param[1]
    q <- qbeta(1 - U, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - U, shape1 = 0.5, shape2 = a))
    logdc <- (dim - 1) * lgamma(a) + lgamma(a + dim / 2) - dim * lgamma(a + 0.5) + (a + 0.5) * sum(log(q + 1)) - (a + dim / 2) * log(sum(q) + 1)
    out <- exp(logdc)
    return(out)
  }

  dcMGL180.reg <- function(U, param) {
    dcMGL.reg(1 - U, param = param)
  }
  dcMGLEV180.reg <- function(U, param) {
    u1 <- U[1]
    u2 <- U[2]
    as.numeric(dcMGLEV180.bivar(u1, u2, param = param[1]))
  }

  dcMGLEV.reg <- function(U, param) {
    u1 <- U[1]
    u2 <- U[2]
    as.numeric(dcMGLEV.bivar(u1, u2, param = param[1]))
  }

  dgumcop.reg <- function(U, param) {
    as.numeric(fCopulae::devCopula(u = U[1], v = U[2], type = "gumbel", param = param[1])) # Bivariate Extreme
  }

  hcMGL.reg <- function(U, param) {
    hfunc1 <- hcMGL.bivar(u1 = U[1], u2 = U[2], pars = param[1])$hfunc1
    hfunc2 <- hcMGL.bivar(u1 = U[1], u2 = U[2], pars = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out
  }

  hcMGL180.reg <- function(U, param) {
    hfunc1 <- 1 - hcMGL.bivar(u1 = 1 - U[1], u2 = 1 - U[2], pars = param[1])$hfunc1
    hfunc2 <- 1 - hcMGL.bivar(u1 = 1 - U[1], u2 = 1 - U[2], pars = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out
  }

  hcMGLEV180.reg <- function(U, param) {
    hfunc1 <- hcMGLEV180.bivar(u1 = U[1], u2 = U[2], param = param[1])$hfunc1
    hfunc2 <- hcMGLEV180.bivar(u1 = U[1], u2 = U[2], param = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out
  }
  hcMGLEV.reg <- function(U, param) {
    hfunc1 <- 1 - hcMGLEV180.bivar(u1 = 1 - U[1], u2 = 1 - U[2], param = param[1])$hfunc1
    hfunc2 <- 1 - hcMGLEV180.bivar(u1 = 1 - U[1], u2 = 1 - U[2], param = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out
  }

  hgumcop.reg <- function(U, param) {
    VineCopula::BiCopHfunc(u1 = U[1], u2 = U[2], family = 4, par = param[1])
  }

  if (copula == "MGL") {
    dcop <- dcMGL.reg
    hcop <- hcMGL.reg
  } else if (copula == "MGL180") {
    dcop <- dcMGL180.reg
    hcop <- hcMGL180.reg
  } else if (copula == "MGL-EV") {
    dcop <- dcMGLEV.reg
    hcop <- hcMGLEV.reg
  } else if (copula == "MGL-EV180") {
    dcop <- dcMGLEV180.reg
    hcop <- hcMGLEV180.reg
  } else if (copula == "Gumbel") {
    dcop <- dgumcop.reg
    hcop <- hgumcop.reg
  }
  # Obs1 <- obs[, 1] # y1 - the first vector of observations
  Obs2 <- obs[, 2] # y2 - the second vector of observations
  f1 <- f[, 1]
  f2 <- f[, 2]
  copLogL <- function(pars, X) {
    ll <- 0
    if (copula == "Gumbel") {
      delta <- exp(X %*% pars) + 1
    } else {
      delta <- exp(X %*% pars)
    }
    for (i in seq_len(nrow(X))) {
      if (Obs2[i] <= umin) {
        ll[i] <- f1[i] * (hcop(U[i, ], param = delta[i])$hfunc1 - hcop(U_[i, ], param = delta[i])$hfunc1)
      } else {
        ll[i] <- f1[i] * f2[i] * dcop(U[i, ], param = delta[i])
      }
    }
    res <- -sum((log(ll)))
    return(res)
  }

  resopt <- optim(
    par = initpar,
    fn = copLogL,
    X = X,
    hessian = hessian, ...
  )
  # resopt

  list(
    loglike = -resopt$value,
    copula = list(name = copula),
    estimates = resopt$par,
    se = sqrt(diag(solve(resopt$hessian))),
    hessian = -resopt$hessian,
    AIC = 2 * length(resopt$par) + 2 * resopt$value,
    BIC = log(nrow(U)) * length(resopt$par) + 2 * resopt$value
  )
}
