



#' Fitting d-dimensional MGL and survival MGL regression models
#'
#' @description \code{MGL.reg.multi} is used to fit d-dimensional MGL and survival MGL regression models.
#' @param U d-dimenstional matrix with values in \eqn{[0,1]}.
#' @param X design matrix.
#' @param copula 'MGL', 'MGL180'.
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
#' * copula: "MGL180" denotes the survival MGL and survival MGL-EV copula respectively.
#' * For "MGL" and "MGL180" regression model, the copula parameter \deqn{\delta_i = \exp(X\beta),} where \eqn{\beta} is the vector of coefficients to be estimated in the copula regression.
#' @export
#'
#'
MGL.reg.multi <- function(U, X, copula = c(
                            "MGL", "MGL180"
                          ),
                          hessian = TRUE, initpar, ...) {
  dcMGL.reg <- function(U, param, log = TRUE) {
    dim <- length(U)
    a <- 1 / param[1]
    q <- qbeta(1 - U, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - U, shape1 = 0.5, shape2 = a))
    logdc <- (dim - 1) * lgamma(a) + lgamma(a + dim / 2) - dim * lgamma(a + 0.5) + (a + 0.5) * sum(log(q + 1)) - (a + dim / 2) * log(sum(q) + 1)

    if (log == TRUE) {
      out <- logdc
    } else {
      out <- exp(logdc)
    }
    return(out)
  }

  dcMGL180.reg <- function(U, param, log = TRUE) {
    dcMGL.reg(1 - U, param = param, log = log)
  }


  if (copula == "MGL") {
    dcop <- dcMGL.reg
  } else if (copula == "MGL180") {
    dcop <- dcMGL180.reg
  }


  copLogL <- function(pars, X) {
    ll <- 0
    delta <- exp(X %*% pars)

    for (i in 1:nrow(X)) {
      ll[i] <- dcop(U[i, ], param = as.vector(delta[i]), log = FALSE)
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
