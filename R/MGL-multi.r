
#' @name cMGL.multi
#' @rdname  cMGL.multi
#' @title d-dimensional MGL and survival MGL copula
#' @description
#' Density, distribution function, and random generation for the MGL and survival MGL copula.
#'
#' @param u d-dimensional matrix
#' @param d d-dimensional
#' @param pars copula parameter, denoted by \eqn{\delta>0}.
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @importFrom stats qbeta
#' @importFrom stats pbeta
#' @importFrom stats runif
#' @importFrom stats pgamma
#' @md
#' @return
#' * \code{dcMGL.multi}, \code{pcMGL.multi} and \code{rcMGL.multi} gives values of Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter \eqn{\delta>0}.
#' * \code{dcMGL180.multi}, \code{pcMGL180.multi} and \code{rcMGL180.multi} gives values of Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter \eqn{\delta>0}.
NULL


#' @rdname cMGL.multi
#' @export
#' @examples
#'
#'
#' dcMGL.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 2, log = FALSE)
#'
#' dcMGL.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = TRUE)
dcMGL.multi <- function(u, pars, log = FALSE) {
  # coding as a matrix
  dim <- ncol(u)
  a <- 1 / pars[1]
  q <- (qbeta(1 - u, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u, shape1 = 0.5, shape2 = a)))
  logdc <- 0
  for (i in 1:nrow(u)) {
    logdc[i] <- (dim - 1) * lgamma(a) + lgamma(a + dim / 2) - dim * lgamma(a + 0.5) + (a + 0.5) * sum(log(q[i, ] + 1)) - (a + dim / 2) * log(sum(q[i, ]) + 1)
  }
  dc <- exp(logdc)
  dc[which(q == Inf)] <- 0
  if (log == TRUE) {
    out <- logdc
  } else {
    out <- dc
  }
  return(out)
}


#' @rdname cMGL.multi
#' @export
#' @examples
#' \dontrun{
#'
#' dcMGL180.multi(u = cbind(c(0.6, 0.1, 0.5), c(0.3, 0.9, 0.2)), pars = 2, log = FALSE)
#'
#' dcMGL180.multi(u = cbind(c(0.6, 0.1), c(0.3, 0.9), c(0.5, 0.6)), pars = 2, log = TRUE)}
dcMGL180.multi <- function(u, pars, log = FALSE){
  dcMGL.multi(u = 1 - u, pars = pars, log = log)
}

#' @rdname cMGL.multi
#' @export
#' @examples
#' # 2-dim MGL copula
#' pcMGL.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
#' # 3-dim MGL copula
#' pcMGL.multi(u = cbind(c(0, 0.2, 0.5), c(0.5, 0.2, 0.5), c(0.01, 0.5, 0.9)), pars = 3)
pcMGL.multi <- function(u, pars) {
  dim <- ncol(u)
  a <- 1 / pars[1]

  q <- (qbeta(1 - u, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u, shape1 = 0.5, shape2 = a)))
  z <- 0
  for (i in 1:nrow(u)) {
    if(any(q[i,] %in% Inf)) {
      z[i] <- 0
    } else {
      temp <- as.vector(q[i, ])
      fin <- function(theta) { # rely on a
          # theta <- as.vector()
          m <- pracma::erfc((temp * theta)^0.5)
          k <- base::prod(m)
          out <- k * theta^(a - 1) * exp(-theta) / gamma(a)
          out
      }
      fin <- Vectorize(fin)
      z[i] <- as.numeric(pracma::integral(
        fun = fin,
        xmin = 0,
        xmax = Inf, method = "Clenshaw"
      ))
    }
  }
  return(z) # rely on k
}

#' @rdname cMGL.multi
#' @export
#' @examples
#' pcMGL180.multi(u = cbind(c(0.5, 0.5), c(0.01, 0.9)), pars = 3)
pcMGL180.multi <- function(u, pars){
  dim <- ncol(u)
  delta <- pars[1]

  q <- (qbeta(u, shape1 = 0.5, shape2 = delta) / (1 - qbeta(u, shape1 = 0.5, shape2 = delta)))
  z <- 0
  for (i in 1:nrow(u)) {
    # if(any(q[i,] %in% Inf)) {
    #   z[i] <- 0
    # } else {
    temp <- as.vector(q[i, ])
    fin <- function(theta) { # rely on theta
        m <- pgamma(q = temp/theta, shape = 0.5, scale = 1)
        k <- base::prod(m)
        out <- k*theta^(-(delta+1))*exp(-1/theta)/gamma(delta)
        out
        }
    fin <- Vectorize(fin)

    z[i] <- as.numeric(pracma::integral(
      fun = fin,
      xmin = 0,
      xmax = Inf, method = "Clenshaw"
    ))
  }
  return(z)

}
# pcMGL180.multi(u = cbind(c(0.5, 0.1), c(0.9, 0.1)), pars = 1.5)
# pcMGB2.bivar(u1 = c(0.5, 0.1), u2 = c(0.9, 0.1), pars1 = 0.5, pars2 = 0.5, pars3 = 1.5)


#' @rdname cMGL.multi
#' @export
#' @examples
#' Usim <- rcMGL.multi(n = 1000, d = 2, pars = 1)
#' plot(Usim)
rcMGL.multi <- function(n, d, pars) {
  a <- 1 / pars
  dseq <- seq(1:(d - 1))
  anew <- c(a, a + dseq / 2)
  Iinv_trans <- function(u, a) {
    qbeta(1 - u, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u, shape1 = 0.5, shape2 = a))
  }
  umat <- matrix(data = runif(n * d, min = 0, max = 1), nrow = n, ncol = d)
  mmat <- qmat <- zmat <- Usim <- matrix(data = 0, nrow = n, ncol = d)

  for (j in 1:ncol(umat)) {
    qmat[, j] <- Iinv_trans(umat[, j], a = anew[j]) # column j denots the qj
  }

  mmat[, 1] <- qmat[, 1]
  mmat[, 2] <- (1 + mmat[, 1]) * qmat[, 2]
  if (d > 2) {
    for (j in 3:ncol(umat)) {
      tempmat <- mmat[, (1):(ncol(umat) - 1)]
      if (n == 1) {
        msum <- sum(tempmat)
      } else {
        msum <- apply(tempmat, MARGIN = 1, sum)
      }
      mmat[, j] <- (1 + msum) * qmat[, j]
    }
  }

  for (j in 1:ncol(zmat)) {
    k <- mmat[, j] / (mmat[, j] + 1)
    Usim[, j] <- 1 - pbeta(k, shape1 = 0.5, shape2 = a)
  }

  return(Usim)
}


#' @rdname cMGL.multi
#' @export
#' @examples
#' Usim <- rcMGL180.multi(n = 1000, d = 2, pars = 1)
#' plot(Usim)
rcMGL180.multi <- function(n, d, pars) {
  1 - rcMGL.multi(n = n, d = d, pars = pars)
}
