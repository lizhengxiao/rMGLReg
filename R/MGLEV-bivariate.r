
#' @name BMGL-EV
#' @rdname  BMGL-EV
#' @title Bivarite MGL-EV and survival MGL-EV copula
#' @description  Density, distribution function, and h-function the bivariate MGL-EV and survival MGL-EV copula with copula parameter \eqn{\delta}.
#' @param u1,u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param param copula parameter, denoted by \eqn{\delta>0}.
#' @md
#' @importFrom stats qbeta
#' @importFrom stats pbeta
#' @details
#' - The extreme value copula \eqn{\bar{C}^{MGL-EV}} of the  survival MGL copula is given by
#' 	\deqn{
#' 		\bar{C}^{MGL-EV}(u_1,u_2;\delta)=\exp\left[\log\left(u_1u_2\right)A_{\delta}\left(\frac{\log\left(u_2\right)}{\log\left(u_1u_2\right)}\right)\right],
#' }
#' where \eqn{A_{\delta}} is the Pickands dependence function implemented by R's \code{\link[rMGLReg]{Afunction}}.
#'
#' - The cdf of the MGL-EV copula is also given by
#'
#' \deqn{
#' 		C^{MGL-EV}(u_1, u_2; {\delta})=\frac{(u_1+u_2)\left(1-A_{\delta}\left(\frac{u_1}{u_1+u_2}\right)\right)}{2\left(1 - A_{\delta}\left(\frac{1}{2}\right)\right)}.
#' }
#'
#' - The h-function is defined as the conditional distribution function of a bivariate copula, i.e.,
#' \deqn{h_1(u_2|u_1,\delta) := P(U_2 \leq u_2 | U_1 = u_1) = \partial C(u_1,u_2) / \partial u_1,}
#'
#' \deqn{h_2(u_1|u_2,\delta) := P(U_1 \leq u_1 | U_2 = u_2) := \partial C(u_1,u_2) / \partial u_2,}
#'
#' where \eqn{(U_1, U_2) \sim C}, and \eqn{C} is a bivariate copula distribution function with parameter(s) \eqn{\delta}. For more details see Aas et al. (2009).
#'
#' @references Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence. Insurance: Mathematics and Economics 44 (2), 182-198.
#'
#' @return
#' * \code{dcMGLEV.bivar}, \code{pcMGLEV.bivar} and \code{hMGLEV.bivar}  give values of density distribution and h-function for the 2-dimensional MGL copula with copula parameter \eqn{\delta>0}.
#' * \code{dcMGLEV180.bivar}, \code{pcMGLEV180.bivar} and \code{hMGLEV180.bivar}  give values of density and distribution and h-function for the 2-dimensional MGL copula with copula parameter \eqn{\delta>0}.
#' *  \code{hMGLEV.bivar} and \code{hMGLEV180.bivar} return a list with
#'
#'    * hfunc1: \eqn{\partial C(u_1,u_2) / \partial u_1,}.
#'    * hfunc2: \eqn{\partial C(u_1,u_2) / \partial u_2,}.
NULL




#' @rdname  BMGL-EV
#' @export
#' @examples
#' pcMGLEV.bivar(u1 = c(0.3, 0.9), u2 = c(0.5, 0.8), param = 2)
pcMGLEV.bivar <- function(u1, u2, param) {
  out <- u1 + u2 - 1 + pcMGLEV180.bivar(1 - u1, 1 - u2, param = param)
  return(out)
}

#' @rdname  BMGL-EV
#' @export
#' @examples
#' dcMGLEV.bivar(u1 = 0.001, u2 = 0.999, param = 1)
dcMGLEV.bivar <- function(u1, u2, param) {
  dcMGLEV180.bivar(1 - u1, 1 - u2, param = param)
}



#' @rdname  BMGL-EV
#' @export
#' @examples
#' pcMGLEV180.bivar(u1 = c(0.3, 0.9), u2 = c(0.5, 0.8), param = 2)
pcMGLEV180.bivar <- function(u1, u2, param) {
  lf <- function(y1, y2, param) {
    delta <- param
    z2 <- (y2)^(-delta) / (y1^(-delta) + y2^(-delta))
    z1 <- (y1)^(-delta) / (y1^(-delta) + y2^(-delta))
    out <- y1 * pbeta(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) + y2 * pbeta(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5)
    return(out)
  }
  # lf <- Vectorize(lf)
  out <- exp(-lf(-log(u1), -log(u2), param = param))
  index1 <- u1 == 0 | u2 == 0
  index2 <- u1 == 1 & u2 != 0
  index3 <- u1 == 0 & u2 == 1
  out[index1] <- 0
  out[index2] <- u2[index2]
  out[index3] <- u1[index3]

  # if(u1==0|u2==0){
  #   out <- 0
  # } else if (u1==1|u2!=0){
  #   out <- u2
  # } else if (u1==0|u2==1){
  #   out <- u1
  # } else {
  #   out <- exp(-lf(-log(u1), -log(u2), param = param))
  # }
  return(out)
}

#' @rdname  BMGL-EV
#' @export
#' @examples
#' dcMGLEV180.bivar(u1 = 0.5, u2 = 0.78, param = 1.8)
dcMGLEV180.bivar <- function(u1, u2, param) {
  lf <- function(y1, y2, param) {
    delta <- param
    z2 <- (y2)^(-delta) / (y1^(-delta) + y2^(-delta))
    z1 <- (y1)^(-delta) / (y1^(-delta) + y2^(-delta))
    out <- y1 * pbeta(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) + y2 * pbeta(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5)
    return(out)
  }

  pbeta_dev <- function(y, shape1, shape2) {
    (1 - y)^(shape2 - 1) * y^(shape1 - 1) / beta(shape1, shape2)
  }

  pbeta_dev2 <- function(y, shape1, shape2) {
    -(1 - y)^(shape2 - 2) * y^(shape1 - 2) * ((y - 1) * shape1 + (shape2 - 2) * y + 1) / beta(shape1, shape2)
  }
  z2_x1 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta * x1^(-delta - 1) * x2^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z1_x1 <- function(x1, x2, param) {
    delta <- param
    part1 <- -delta * x1^(-delta - 1) * x2^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z2_x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- -delta * x2^(-delta - 1) * x1^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z1_x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta * x2^(-delta - 1) * x1^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z2_x1x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta^2 * (x1 * x2)^(-delta - 1) * (x2^(-delta) - x1^(-delta))
    part2 <- (x2^(-delta) + x1^(-delta))^3
    part1 / part2
  }

  z1_x1x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta^2 * (x1 * x2)^(-delta - 1) * (x1^(-delta) - x2^(-delta))
    part2 <- (x2^(-delta) + x1^(-delta))^3
    part1 / part2
  }



  lf_dev_x1 <- function(x1, x2, param) {
    delta <- param
    z2 <- x2^(-delta) / (x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    out <- pbeta(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) + x1 * pbeta_dev(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z2_x1(x1 = x1, x2 = x2, param = delta) + x2 * pbeta_dev(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z1_x1(x1 = x1, x2 = x2, param = delta)
    out
  }


  lf_dev_x2 <- function(x1, x2, param) {
    delta <- param
    z2 <- x2^(-delta) / (x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    out <- pbeta(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5) + x2 * pbeta_dev(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z1_x2(x1 = x1, x2 = x2, param = delta) + x1 * pbeta_dev(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z2_x2(x1 = x1, x2 = x2, param = delta)
    out
  }

  lf_dev2 <- function(x1, x2, param) {
    delta <- param
    z2 <- x2^(-delta) / (x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    If <- function(z) {
      pbeta(z, shape1 = 0.5, shape2 = 1 / delta + 0.5)
    }
    If_dev <- function(z) {
      pbeta_dev(y = z, shape1 = 0.5, shape2 = 1 / delta + 0.5)
    }
    If_dev2 <- function(z) {
      pbeta_dev2(y = z, shape1 = 0.5, shape2 = 1 / delta + 0.5)
    }
    p1 <- If_dev(z2) * z2_x2(x1 = x1, x2 = x2, param = delta)
    p2 <- x1 * (If_dev2(z2) * z2_x2(x1 = x1, x2 = x2, param = delta) * z2_x1(x1 = x1, x2 = x2, param = delta)
      + If_dev(z2) * z2_x1x2(x1 = x1, x2 = x2, param = delta))
    p3 <- ((If_dev(z1) + x2 * If_dev2(z1) * z1_x2(x1 = x1, x2 = x2, param = delta)) * z1_x1(x1 = x1, x2 = x2, param = delta)
      + x2 * If_dev(z1) * z1_x1x2(x1 = x1, x2 = x2, param = delta))
    out <- p1 + p2 + p3
    out
  }
  # final results
  p1 <- pcMGLEV180.bivar(u1 = u1, u2 = u2, param = param) / (u1 * u2)
  # p1 <- exp(-lf(y1 = -log(u1), y2 = -log(u2), param = param))/(u1*u2)
  p2 <- lf_dev_x1(-log(u1), -log(u2), param = param) * lf_dev_x2(-log(u1), -log(u2), param = param) - lf_dev2(-log(u1), -log(u2), param = param)
  out <- p1 * p2
  return(out)
}



#' @rdname  BMGL-EV
#' @export
#' @examples
#' hcMGLEV180.bivar(u1 = 0.5, u2 = 0.78, param = 1.8)
hcMGLEV180.bivar <- function(u1, u2, param) {
  delta <- param
  lf <- function(y1, y2, param) {
    delta <- param
    z2 <- (y2)^(-delta) / (y1^(-delta) + y2^(-delta))
    z1 <- (y1)^(-delta) / (y1^(-delta) + y2^(-delta))
    out <- y1 * pbeta(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) + y2 * pbeta(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5)
    return(out)
  }

  pbeta_dev <- function(y, shape1, shape2) {
    (1 - y)^(shape2 - 1) * y^(shape1 - 1) / beta(shape1, shape2)
  }

  pbeta_dev2 <- function(y, shape1, shape2) {
    -(1 - y)^(shape2 - 2) * y^(shape1 - 2) * ((y - 1) * shape1 + (shape2 - 2) * y + 1) / beta(shape1, shape2)
  }
  z2_x1 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta * x1^(-delta - 1) * x2^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z1_x1 <- function(x1, x2, param) {
    delta <- param
    part1 <- -delta * x1^(-delta - 1) * x2^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z2_x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- -delta * x2^(-delta - 1) * x1^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z1_x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta * x2^(-delta - 1) * x1^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1 / part2
  }

  z2_x1x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta^2 * (x1 * x2)^(-delta - 1) * (x2^(-delta) - x1^(-delta))
    part2 <- (x2^(-delta) + x1^(-delta))^3
    part1 / part2
  }

  z1_x1x2 <- function(x1, x2, param) {
    delta <- param
    part1 <- delta^2 * (x1 * x2)^(-delta - 1) * (x1^(-delta) - x2^(-delta))
    part2 <- (x2^(-delta) + x1^(-delta))^3
    part1 / part2
  }



  lf_dev_x1 <- function(x1, x2, param) {
    delta <- param
    z2 <- x2^(-delta) / (x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    out <- pbeta(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) + x1 * pbeta_dev(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z2_x1(x1 = x1, x2 = x2, param = delta) + x2 * pbeta_dev(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z1_x1(x1 = x1, x2 = x2, param = delta)
    out
  }


  lf_dev_x2 <- function(x1, x2, param) {
    delta <- param
    z2 <- x2^(-delta) / (x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    out <- pbeta(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5) + x2 * pbeta_dev(z1, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z1_x2(x1 = x1, x2 = x2, param = delta) + x1 * pbeta_dev(z2, shape1 = 0.5, shape2 = 1 / delta + 0.5) * z2_x2(x1 = x1, x2 = x2, param = delta)
    out
  }

  # final results
  Cu1u2 <- pcMGLEV180.bivar(u1 = u1, u2 = u2, param = delta)
  index1 <- u1 == 0 & u2 != 0
  index2 <- u1 != 0 & u2 == 0
  # index3 <- u1 == 0 & u2 == 0
  hfunc1 <- Cu1u2 / u1 * lf_dev_x1(x1 = -log(u1), x2 = -log(u2), param = delta)
  hfunc2 <- Cu1u2 / u2 * lf_dev_x2(x1 = -log(u1), x2 = -log(u2), param = delta)

  # utemp <- cbind(u1, u2)
  u10 <- rep(1e-10, length(u1[index1]))
  u20 <- rep(1e-10, length(u2[index2]))
  hfunc1[index1] <- pcMGLEV180.bivar(u1 = u10, u2 = u2[index1], param = delta) / u10 * lf_dev_x1(x1 = -log(u10), x2 = -log(u2[index1]), param = delta)
  hfunc2[index1] <- 0.0000000001

  hfunc1[index2] <- 0.0000000001
  hfunc2[index2] <- pcMGLEV180.bivar(u1 = u1[index2], u2 = u20, param = delta) / u20 * lf_dev_x2(x1 = -log(u1[index2]), x2 = -log(u20), param = delta)

  list(hfunc1 = hfunc1, hfunc2 = hfunc2) # h_1(u_2|u_1,??)
}




#' @rdname  BMGL-EV
#' @export
#' @examples
#' hcMGLEV.bivar(u1 = 0.5, u2 = 0.78, param = 1.8)
hcMGLEV.bivar <- function(u1, u2, param) {
  hfunc1 <- 1 - hcMGLEV180.bivar(u1 = 1 - u1, u2 = 1 - u2, param = param)$hfunc1
  hfunc2 <- 1 - hcMGLEV180.bivar(u1 = 1 - u1, u2 = 1 - u2, param = param)$hfunc2
  list(hfunc1 = hfunc1, hfunc2 = hfunc2)
}





#' @title  the Pickands dependence function \eqn{A} in survival MGL-EV copula
#'
#' @param w \eqn{w \in [0,1]}.
#' @param param copula parameter in survival MGL-EV copula.
#'
#' @return the value of the Pickands dependence function in \eqn{[0.5,1]}
#' @md
#' @details
#' - In the bivariate case
#' \eqn{d=2} the stable tail dependence function \eqn{\ell} can be represented in terms of the Pickands dependence function \eqn{A}:
#'
#' \deqn{
#' 		\begin{align}
#' 		A_{\delta}\left(w\right)=w I_{{\frac{1}{2}, \frac{1}{\delta}+\frac{1}{2}}}\left[\frac{\left(1-w \right)^{-\delta}}{\left(1-w \right)^{-\delta} + w^{-\delta} }\right]  + \left(1-w\right) I_{{\frac{1}{2}, \frac{1}{\delta}+\frac{1}{2}}}\left[\frac{w^{-\delta}}{\left(1-w \right)^{-\delta} + w^{-\delta} }\right].
#' 		\end{align}
#' }
#' where \eqn{A_{\delta}} is the Pickands dependence function of  survival MGL-EV copula.
#' Here \eqn{I_{m,n}^{-1}()} denotes the inverse of the beta cumulative distribution function \eqn{I_{m,n}()} (or regularized incomplete beta function)
#' with parameters shape1 = m and shape2 = n
#' implemented by R's \code{\link[stats]{qbeta}} and \code{\link[stats]{pbeta}} respectively.
#'
#'
#' - The extreme value copula  \eqn{\bar{C}^{MGL-EV}} of the  survival MGL copula is given by
#' 	\deqn{
#' 		\bar{C}^{MGL-EV}(u_1,u_2;\delta)=\exp\left[\log\left(u_1u_2\right)A_{\delta}\left(\frac{\log\left(u_2\right)}{\log\left(u_1u_2\right)}\right)\right].
#' }
#'
#' @importFrom stats qbeta
#' @importFrom stats pbeta
#'
#' @export
#' @examples
#' Afunction(w = 0.7, param = 1.2)
#'
#' delta <- 0.7
#' w.vector <- seq(0.0001, 0.9999, length.out = 100)
#' y.vector <- Afunction(w = w.vector, param = delta)
#' plot(w.vector, y.vector, lwd = 2, col = 3, ylim = c(0.5, 1))
#'
#' delta <- 1
#' w.vector <- seq(0.0001, 0.9999, length.out = 100)
#' y.vector <- Afunction(w = w.vector, param = delta)
#' lines(w.vector, y.vector, lwd = 2, col = 4)
#'
#' delta <- 3
#' w.vector <- seq(0.0001, 0.9999, length.out = 100)
#' y.vector <- Afunction(w = w.vector, param = delta)
#' lines(w.vector, y.vector, lwd = 2, col = 5)
#'
#' delta <- 5
#' w.vector <- seq(0.0001, 0.9999, length.out = 100)
#' y.vector <- Afunction(w = w.vector, param = delta)
#' lines(w.vector, y.vector, lwd = 2, col = 1)
#'
Afunction <- function(w, param){
  delta <- param
  k <- - delta
  z2 <- (1 - w)^k/(w^k + (1 - w)^k)
  z1 <- (w)^k/(w^k + (1 - w)^k)
  out <- w*pbeta(z2, shape1 = 0.5, shape2 = 1/delta + 0.5) + (1-w)*pbeta(z1, shape1 = 0.5, shape2 = 1/delta + 0.5)
  return(out)
}



