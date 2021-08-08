



#' Fitting bivariate MGL copula regression models
#'
#' @description MGL.reg is used to fit bivariate MGL copula regression models.
#' @param U two-dimenstional matrix with values in \eqn{\left\[0,1\right\]}.
#' @param X design matrix
#' @param copula 'MGL', 'MGL180', "MGL-EV", "MGL-EV180", "Gumbel"
#' @param initpar Initial values for the parameters to be optimized over.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param ... 	Further arguments to be passed to fn and gr in optiom.
#' @importFrom stats qbeta
#' @importFrom stats optim
#' @return For MGL.reg, a list with components:
#' @export
#'
#'
MGL.reg <- function(U, X, copula = c(
                      "MGL", "MGL180", "MGL-EV",
                      "MGL-EV180",
                      "Gumbel"
                    ),
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
    } else {
      delta <- exp(X %*% pars)
    }
    for (i in 1:nrow(X)) {
      ll[i] <- dcop(U[i, ], param = as.vector(delta[i]))
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

  # resopt
  # list(
  #   loglike = -resopt$minimum,
  #   copula = list(name = copula),
  #   estimates = resopt$estimate,
  #   se = sqrt(diag(solve(resopt$hessian))),
  #   AIC = 2 * length(resopt$estimate) + 2 * resopt$minimum,
  #   BIC = log(nrow(U)) * length(resopt$estimate) + 2 * resopt$minimum
  # )
  # modout <- function(m) {
  #   Hessian <- -m$hessian
  #   se <- sqrt(diag(solve(Hessian)))                  ## stardard error
  #   Z <- m$par/se                                             ##  Z statistic
  #   p <- ifelse(Z>=0, pnorm(Z, lower=F)*2, pnorm(Z)*2)           ## p value
  #   summarytable <- data.frame(m$par, se, Z, p)
  #   LL <- m$value
  #   NLL <- -m$value  # - loglikelihood value
  #   AIC <- 2*length(m$par) + 2*NLL
  #   BIC <- log(nrow(U))*length(m$par) + 2*NLL
  #   list(summary = summarytable, ll =  m$value, AIC = AIC, BIC = BIC)
  # }
  #
  # modout(resopt)
}
