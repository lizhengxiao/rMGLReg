
#' @title Fitting bivariate MGL copula models
#' @description MGL.reg is used to fit bivariate MGL copula regression models.
#' @param U two-dimenstional matrix with values in \eqn{\left\[0,1\right\]}.
#' @param copula copula 'MGL', 'MGL180', "MGL-EV", "MGL-EV180", "MGB2", "Normal" , "Student-t"
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param initpar Initial values for the parameters to be optimized over.
#' @param ... additional arguments to be passed to f.
#' @param U_ asd
#' @param obs two-dimenstional matrix for observations
#' @param f values of the density function for marginal distribution
#' @param umin sd
#' @importFrom stats nlm
#' @return A list containing the following components
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

  hnormcop <- function(U, param){
    VineCopula::BiCopHfunc(U[,1], U[,2], family = 1, par = param[1])}
  htcop <- function(U, param)
    VineCopula::BiCopHfunc(U[,1], U[,2], family = 2, par = param[1], par2 = param[2])
  hgumcop <- function(U, param)
    VineCopula::BiCopHfunc(U[,1], U[,2], family = 4, par = param[1])
  hMGL <- function(U, param){
    hcMGL.bivar(u1 = U[,1], u2 = U[,2], pars = param[1])}
  hMGL180 <- function(U, param){
    hfunc1 <- 1 - hcMGL.bivar(u1 = 1 - U[,1], u2 = 1 - U[,2], pars = param[1])$hfunc1
    hfunc2 <- 1 - hcMGL.bivar(u1 = 1 - U[,1], u2 = 1 - U[,2], pars = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out}
  hMGLEV180 <- function(U, param){
    hfunc1 <- hMGLEV180.bivar(u1 = U[,1], u2 = U[,2], param = param[1])$hfunc1
    hfunc2 <- hMGLEV180.bivar(u1 = U[,1], u2 = U[,2], param = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out}
  hMGLEV <- function(U, param){
    hfunc1 <- 1 - hMGLEV180.bivar(u1 = 1 - U[,1], u2 = 1 - U[,2], param = param[1])$hfunc1
    hfunc2 <- 1 - hMGLEV180.bivar(u1 = 1 - U[,1], u2 = 1 - U[,2], param = param[1])$hfunc2
    out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
    out}
  hMGB2 <- function(U, param){
    hcMGB2.bivar(u1 = U[,1], u2 = U[,2], pars1 = param[1], pars2 = param[2], pars3 = param[3])}

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
    hcop <-   hMGL
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
    hcop <-   hgumcop
  } else if (copula == "Normal") {
    dcop <- dnormcop
    arg.cop <- argnorm
    hcop <- hnormcop
  } else if (copula == "t") {
    dcop <- dtcop
    arg.cop <- argt
    hcop <-    htcop
  } else if (copula == "MGB2") {
    dcop <- dMGB2
    arg.cop <- argMGB2
    hcop <-   hMGB2
  }

  # loglike.copula <- function(U, initpar, ...){
  Obs1 <- obs[,1] # y1
  Obs2 <- obs[,2] # y2
  f1 <- f[,1]
  f2 <- f[,2]
  copLogL <- function(x) {
    if(all(arg.cop$lower < x && arg.cop$upper > x)){
      index1 <- which(Obs2 <= umin)
      index2 <- which(Obs2 > umin)
      m1 <- f1[index1]*(hcop(U[index1,], param = x)$hfunc1 - hcop(U_[index1,], param = x)$hfunc1)
      logL1 <- (log(m1))
      m2 <- f1[index2]*f2[index2]*dcop(U[index2,], param = x)
      logL2 <- (log(m2))
      ll <- c(logL1, logL2)
      res <- -sum((ll)) # 定义似然函数
    } else { res <- 100000000000000000 }
    return(res)}
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
