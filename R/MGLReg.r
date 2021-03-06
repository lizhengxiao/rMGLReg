



#' Fitting bivariate/multivariate MGL and MGL-EV copula regression models.
#'
#' @description \code{MGL.reg} is used to fit bivariate MGL and MGL-EV copula regression models for two continuous variables.
#' @param U two-dimensional matrix with values in \eqn{[0,1]}.
#' @param X design matrix
#' @param copula 'MGL', 'MGL180', "MGL-EV", "MGL-EV180", "Gumbel", "MGB2".
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
#' The estimation method is performed via \code{\link[stats]{optim}} function. Y1 and Y2 are both continuous variables.
#'
#' copula: "MGL180" and "MGLEV180" denote the survival MGL and survival MGL-EV copula respectively.
#' * For "Gumbel" regression model, the copula parameter \deqn{\delta_i = \exp(X_i^T\beta) + 1,} where \eqn{{X}_{i}=(1,x_{i1}...,x_{ik})} denotes the vector of covariates and \eqn{\beta} is the vector of coefficients to be estimated in the copula regression.
#' * For "MGL", "MGL180", "MGL-EV", "MGL-EV180" regression model, the copula parameter \deqn{\delta_i = \exp(X_i^T\beta).}
#' * For "MGB2", the copula parameter \deqn{\q_i = \exp(X_i^T\beta)} and \eqn{(p_1,p_2)} remain to be constant.
#' * Note that the regression modelling can be extended to the high-dimensional case when copula is "MGL180" and "MGL".
#'
#' @examples
#' # 10-dimensional regression models
#' set.seed(111)
#' d <- 10
#' n <- 1000 # sample size
#' beta.true <- c(-0.6, 0.5, 0.2) # true regression coefficients
#' x1 <- rnorm(n, 0, 1)
#' x2 <- rnorm(n, 0, 1)
#' X <- model.matrix(~ x1 + x2) # design matrix
#' delta.sim <- as.vector(exp(X%*%beta.true)) # true copula parameters
#' Usim <- matrix(0, nrow = n, ncol = d)
#' for (i in 1:n){
#'   Usim[i, ] <- rcMGL.multi(n = 1, d = d, pars = delta.sim[i])
#' }
#' m.MGLMGA <- MGL.reg(U = Usim, copula = "MGL",
#'                        X = X, method = "Nelder-Mead",
#'                        initpar = c(-0.32, 0.001, 0.001)
#' )
#' m.MGLMGA
#'
#' @export
#'
#'
MGL.reg <- function(U, X, copula = c(
                      "MGL", "MGL180", "MGL-EV",
                      "MGL-EV180",
                      "Gumbel", "MGB2"),
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


  if (copula == "MGL") {
    dcop <- dcMGL.reg
  } else if (copula == "MGL180") {
    dcop <- dcMGL180.reg
  } else if (copula == "MGL-EV") {
    dcop <- dcMGLEV.reg
  } else if (copula == "MGL-EV180") {
    dcop <- dcMGLEV180.reg
  } else if (copula == "Gumbel") {
    dcop <- dgumcop.reg
  }

  copLogL <- function(pars, X) {
    ll <- 0

    if (copula == "Gumbel") {
      delta <- exp(X %*% pars) + 1
      for (i in seq_len(nrow(X))) {ll[i] <- dcop(U[i, ], param = as.vector(delta[i]))}
    } else if (copula == "MGB2"){
      p1 <- exp(pars[ncol(X) + 1])
      p2 <- exp(pars[ncol(X) + 2])
      q <- exp(X%*%pars[1:(ncol(X))])
      for (i in seq_len(nrow(X))) {ll[i] <- dcMGB2.bivar(u1 = U[i,1], u2 = U[i,2], pars1 = p1, pars2 = p2, pars3 = q[i])}
    } else {
      delta <- exp(X %*% pars)
      for (i in seq_len(nrow(X))) {ll[i] <- dcop(U[i, ], param = as.vector(delta[i]))}
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


  if (hessian == TRUE){
  out <- list(
          loglike = -resopt$value,
          copula = list(name = copula),
          estimates = resopt$par,
          se = sqrt(diag(solve(resopt$hessian))),
          hessian = -resopt$hessian,
          AIC = 2 * length(resopt$par) + 2 * resopt$value,
          BIC = log(nrow(U)) * length(resopt$par) + 2 * resopt$value)
  } else {
    out <- list(
      loglike = -resopt$value,
      copula = list(name = copula),
      estimates = resopt$par,
      AIC = 2 * length(resopt$par) + 2 * resopt$value,
      BIC = log(nrow(U)) * length(resopt$par) + 2 * resopt$value)
  }
  out

}
