
#' @title Fitting bivariate MGL copula models
#' @description \code{MGL.mle.mixed} is used to fit bivariate mixed copula regression models via maximum likelihood (ML) method for continuous and semi-continuous variables..
#' @param U two-dimenstional matrix for pseudo copula data with values in \eqn{[0,1]} for (F(y1), F(y2)).
#' @param copula copula 'MGL', 'MGL180', "MGL-EV", "MGL-EV180", "MGB2", "Normal" , "Student-t".
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param initpar Initial values for the parameters to be optimized over.
#' @param U_ two-dimenstional matrix for pseudo copula data for the data (F(y1), F(y2-1)).
#' @param obs two-dimenstional matrix for loss observations (y1, y2).
#' @param f values of the density function for marginal distribution.
#' @param umin threshold value used in the semi-continuous data.
#' @param ... additional arguments, see \code{\link[stats]{nlm}} for more details.
#' @importFrom stats nlm
#' @md
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
#'
#' For a portfolio of \eqn{n} observations \eqn{(y_{i1},y_{i2}; \; i=1,\ldots,n)}, the joint density function of \eqn{(Y_1,Y_2)} can be written as
#' \deqn{
#' 	f_{Y_{1},Y_2}(y_{i1},y_{i2})=\begin{cases}
#' 	f_{Y_1}(y_{i1})[
#' 	h_{2|1}(F_{Y_{1}}(y_{i1}),F_{Y_{2}}(y_{i2})) - h_{2|1}(F_{Y_{1}}(y_{i1}),F_{Y_{2}}(y_{i2}-1))
#' 	], & y_{i2}\le umin,\\
#' 	f_{Y_1}(y_{i1})f_{Y_2}(y_{i2})c(F_{Y_{1}}(y_{i1}), F_{Y_{2}}(y_{i2})), & y_{i2} > umin,
#' 	\end{cases}
#' }
#' where the density \eqn{f_{Y_j}(\cdot)} and cdf \eqn{F_{Y_j}(\cdot)} of the marginal distributions  (\eqn{i=1,2}) are specified respectively. Here
#' \eqn{h_{2|1}(u_1, u_2)=\partial C(u_1,u_2)/\partial u_1} is the \eqn{h}-function of bivariate copula.
#'
#'
#' @export
#'
MGL.mle.mixed <- function(obs, U, U_, f, copula = c(
                            "MGL", "MGL180", "MGL-EV",
                            "MGL-EV180",
                            "Gumbel",
                            "Normal", "MGB2", "t"
                          ), umin,
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

  hnormcop <- function(U, param) {
    VineCopula::BiCopHfunc(U[, 1], U[, 2], family = 1, par = param[1])
  }
  htcop <- function(U, param) {
    VineCopula::BiCopHfunc(U[, 1], U[, 2], family = 2, par = param[1], par2 = param[2])
  }
  hgumcop <- function(U, param) {
    VineCopula::BiCopHfunc(U[, 1], U[, 2], family = 4, par = param[1])
  }
  hMGL <- function(U, param) {
    hcMGL.bivar(u1 = U[, 1], u2 = U[, 2], pars = param[1])
  }
  hMGL180 <- function(U, param) {
    hfunc1 <- 1 - hcMGL.bivar(u1 = 1 - U[, 1], u2 = 1 - U[, 2], pars = param[1])$hfunc1
    hfunc2 <- 1 - hcMGL.bivar(u1 = 1 - U[, 1], u2 = 1 - U[, 2], pars = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out
  }
  hMGLEV180 <- function(U, param) {
    hfunc1 <- hcMGLEV180.bivar(u1 = U[, 1], u2 = U[, 2], param = param[1])$hfunc1
    hfunc2 <- hcMGLEV180.bivar(u1 = U[, 1], u2 = U[, 2], param = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out
  }
  hMGLEV <- function(U, param) {
    hfunc1 <- 1 - hcMGLEV180.bivar(u1 = 1 - U[, 1], u2 = 1 - U[, 2], param = param[1])$hfunc1
    hfunc2 <- 1 - hcMGLEV180.bivar(u1 = 1 - U[, 1], u2 = 1 - U[, 2], param = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out
  }
  hMGB2 <- function(U, param) {
    hcMGB2.bivar(u1 = U[, 1], u2 = U[, 2], pars1 = param[1], pars2 = param[2], pars3 = param[3])
  }

  argnorm <- list(length = 1, lower = 0, upper = 1, name = "Gaussian")
  argt <- list(length = 2, lower = c(0, 0), upper = c(1, 100), name = "Student")
  arggum <- list(length = 1, lower = 1, upper = 50, name = "Gumbel")
  argMGB2 <- list(length = 3, lower = c(0, 0, 0), upper = c(50, 50, 50), name = "MGB2")
  argMG180 <- list(length = 1, lower = 0, upper = 10, name = "MGL180")
  argMGLEV180 <- list(length = 1, lower = 0, upper = 10, name = "MGL-EV180")
  argMG <- list(length = 1, lower = 0, upper = 10, name = "MGL")
  argMGLEV <- list(length = 1, lower = 0, upper = 10, name = "MGL-EV")
  if (copula == "MGL") {
    dcop <- dMGL
    arg.cop <- argMG
    hcop <- hMGL
  } else if (copula == "MGL180") {
    dcop <- dMGL180
    arg.cop <- argMG180
    hcop <- hMGL180
  } else if (copula == "MGL-EV") {
    dcop <- dMGLEV
    arg.cop <- argMGLEV
    hcop <- hMGLEV
  } else if (copula == "MGL-EV180") {
    dcop <- dMGLEV180
    arg.cop <- argMGLEV180
    hcop <- hMGLEV180
  } else if (copula == "Gumbel") {
    dcop <- dgumcop
    arg.cop <- arggum
    hcop <- hgumcop
  } else if (copula == "Normal") {
    dcop <- dnormcop
    arg.cop <- argnorm
    hcop <- hnormcop
  } else if (copula == "t") {
    dcop <- dtcop
    arg.cop <- argt
    hcop <- htcop
  } else if (copula == "MGB2") {
    dcop <- dMGB2
    arg.cop <- argMGB2
    hcop <- hMGB2
  }


  # Obs1 <- obs[, 1] # y1
  Obs2 <- obs[, 2] # y2
  f1 <- f[, 1]
  f2 <- f[, 2]
  copLogL <- function(x) {
    if (all(arg.cop$lower < x && arg.cop$upper > x)) {
      index1 <- which(Obs2 <= umin)
      index2 <- which(Obs2 > umin)
      m1 <- f1[index1] * (hcop(U[index1, ], param = x)$hfunc1 - hcop(U_[index1, ], param = x)$hfunc1)
      logL1 <- (log(m1))
      m2 <- f1[index2] * f2[index2] * dcop(U[index2, ], param = x)
      logL2 <- (log(m2))
      ll <- c(logL1, logL2)
      res <- -sum((ll)) # define the loglikelihood
    } else {
      res <- 100000000000000000
    }
    return(res)
  }

  resopt <- nlm(
    f = copLogL,
    p = initpar,
    hessian = hessian
  )


  list(
    loglike = -resopt$minimum,
    copula = list(name = arg.cop$name),
    estimates = resopt$estimate,
    se = sqrt(diag(solve(resopt$hessian))),
    hessian = -resopt$hessian,
    AIC = 2 * length(resopt$estimate) + 2 * resopt$minimum,
    BIC = log(nrow(U)) * length(resopt$estimate) + 2 * resopt$minimum
  )
}
