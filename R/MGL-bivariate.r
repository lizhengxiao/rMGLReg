#' @name BMGL
#' @rdname  BMGL
#' @title Bivarite MGL and survival MGL copula
#' @description Density, distribution function, and random generation for the bivariate MGL and survival MGL copula with copula parameter \eqn{\delta}.
#' @param u1,u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param pars copula parameter, denoted by \eqn{\delta >0}.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @md
#' @details
#' - For \eqn{d=2}, the MGL copula density function is given by
#'
#' \eqn{
#' c^{\text{MGL}}(u_1, \dots, u_d;a)=\frac{\Gamma(a)^{d-1}\Gamma(a+\frac{d}{2})}{\Gamma(a+\frac{1}{2})^d}\frac{\prod_{j=1}^{d}(t(u_j;a)+1)^{a+\frac{1}{2}}}{\left(\sum_{j=1}^{d}t(u_j;a)+1\right)^{a+\frac{d}{2}}},
#' \quad (u_1, ..., u_d)\in \left[0,1\right]^d,
#' }
#'
#' where \eqn{a = \frac{1}{\delta}}, and
#' \eqn{I_{m,n}()}  denotes the beta cumulative distribution function (or regularized incomplete beta function)
#' with parameters shape1 = m and shape2 = n
#' implemented by R's \code{\link[stats]{pbeta}}.
#'
#' - The MGL copula distribution function has the form of
#'
#' \deqn{
#'  C^{\text{MGL}}(u_1, \dots, u_d;a) =\mathbb{E}_{\Theta}\left[\prod_{j=1}^{d}\text{erfc}\left(
#'  \sqrt{t(u_j;a)\Theta}
#'  \right)\right]
#'  }
#' with \eqn{\Theta \sim \text{gamma}(a,1)} and
#' \deqn{
#'  t(u_j;a)=\frac{I^{-1}_{\frac{1}{2},{a}}(1-u_j)}{1-I^{-1}_{\frac{1}{2},{a}}(1-u_j)}.
#'  }
#'
#' - For \eqn{d=2}, the density and distribution function of survival MGL copula is
#'  \deqn{
#'  	\bar{c}^{MGL}(u_1,u_2;\delta)=c^{MGL}(1-u_1,1-u_2;\delta),
#'  }
#'  and
#'  \deqn{
#'  	\bar{C}^{MGL}(u_1,u_2;\delta)=u_1+u_2-1+C^{MGL}(1-u_1,1-u_2;\delta).
#'  }
#'
#'
#' @return
#' - \code{dcMGL.bivar}, \code{pcMGL.bivar} and \code{rcMGL.bivar} gives values of Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter \eqn{\delta>0}.
#' - \code{dcMGL180.bivar}, \code{pcMGL180.bivar} and \code{rcMGL180.bivar} gives values of Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter \eqn{\delta>0}.
NULL




#' @rdname  BMGL
#' @export
#' @examples
#' dcMGL.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
dcMGL.bivar <- function(u1, u2, pars) {
  dim <- 2
  a <- 1 / pars
  q1 <- qbeta(1 - u1, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u1, shape1 = 0.5, shape2 = a))
  q2 <- qbeta(1 - u2, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u2, shape1 = 0.5, shape2 = a))
  # if(u1 == 0|u2 == 0|q1 == Inf|q2 == Inf) dc <- 0
  if (q1 == Inf | q2 == Inf) {
    dc <- 0
  } else {
    # q <- cbind(q1, q2)
    # dc <- gamma(a)^(dim - 1)*gamma(a + dim/2)/(gamma(a + 0.5)^dim)*(q1 + 1)^(a + 0.5)*(q2 + 1)^(a + 0.5)/((q1 + q2) + 1)^(a + dim/2)
    dc <- gamma(a)^(dim - 1) / (gamma(a + 0.5)^(dim - 1)) / (gamma(a + 0.5)) * gamma(a + dim / 2) * (q1 + 1)^(a + 0.5) * (q2 + 1)^(a + 0.5) / ((q1 + q2) + 1)^(a + dim / 2)
  }
  return(dc)
}
dcMGL.bivar <- Vectorize(dcMGL.bivar)



#' @rdname BMGL
#' @export
#' @examples
#' pcMGL.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
pcMGL.bivar <- function(u1, u2, pars) {
  a <- 1 / pars
  q1 <- qbeta(1 - u1, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u1, shape1 = 0.5, shape2 = a))
  q2 <- qbeta(1 - u2, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u2, shape1 = 0.5, shape2 = a))
  fin <- function(theta) { # rely on a
    m1 <- pracma::erfc((q1 * theta)^0.5)
    m2 <- pracma::erfc((q2 * theta)^0.5)
    out <- m1 * m2 * theta^(a - 1) * exp(-theta) / gamma(a)
    out
  }
  if (u1 == 0 | u2 == 0) {
    z <- 0
  } else if (u1 == 1 & u2 != 0) {
    z <- u2
  } else if (u1 != 0 & u2 == 1) {
    z <- u1
  } else {
    z <- as.numeric(pracma::integral(fun = fin, xmin = 0, xmax = Inf, method = "Clenshaw"))
  }
  return(z) # rely on k
}
pcMGL.bivar <- Vectorize(pcMGL.bivar)

#' @rdname BMGL
#' @export
#' @examples
#' head(rcMGL.bivar(n = 200, pars = 0.8))
rcMGL.bivar <- function(n, pars) {
  Usim <- rcMGL.multi(n = n, pars = pars, d = 2)
  Usim
}

#' @rdname  BMGL
#' @export
#' @examples
#' dcMGL180.bivar(u1 = 0.001, u2 = 0.999, pars = 1)
dcMGL180.bivar <- function(u1, u2, pars) {
  dcMGL.bivar(1 - u1, 1 - u2, pars)
}

#' @rdname BMGL
#' @export
#' @examples pcMGL180.bivar(u1 = c(0.6, 0.1, 0.7), u2 = c(0.3, 0.6, 0.9), pars = 5)
pcMGL180.bivar <- function(u1, u2, pars) {
  u1 + u2 - 1 + pcMGL.bivar(1 - u1, 1 - u2, pars)
}



#' @rdname BMGL
#' @export
#' @examples rcMGL180.bivar(n = 200, pars = 0.8)
rcMGL180.bivar <- function(n, pars) {
  Usim <- rcMGL.multi(n = n, pars = pars, d = 2)
  1 - Usim
}


#' @name hBMGL
#' @rdname  hBMGL
#' @title Conditional Distribution Function of a Bivariate MGL and survival Copula
#'
#' @param u1 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param u2 numeric vectors of equal length with values in \eqn{[0,1]}.
#' @param pars numeric; single number or vector of size length(u1); copula parameter > 0.
#' @md
#' @details
#' The h-function is defined as the conditional distribution function of a bivariate copula, i.e.,
#' \deqn{h_1(u_2|u_1,\delta) := P(U_2 \leq u_2 | U_1 = u_1) = \partial C(u_1,u_2) / \partial u_1,}
#'
#' \deqn{h_2(u_1|u_2,\delta) := P(U_1 \leq u_1 | U_2 = u_2) := \partial C(u_1,u_2) / \partial u_2,}
#'
#' where \eqn{(U_1, U_2) \sim C}, and \eqn{C} is a bivariate copula distribution function with parameter(s) \eqn{\delta}. For more details see Aas et al. (2009).
#'
#' - The h-function corresponding to the bivariate MGL copula is given by
#' 	\deqn{
#' 	\begin{align}
#' 	h_{2|1}^{\text{MGL}}(u_2|u_1;\delta)&= 1-I_{\frac{1}{2},\frac{1}{\delta}+\frac{1}{2}}\left[
#' 	  \frac{t(u_2;\delta)}{t(u_1;\delta) + t(u_2;\delta)+1}
#' 	  \right],\\
#' 	h_{1|2}^{\text{MGL}}(u_1|u_2;\delta) &= 1-I_{\frac{1}{2},\frac{1}{\delta}+\frac{1}{2}}\left[
#' 	  \frac{t(u_1;\delta)}{t(u_1;\delta) + t(u_2;\delta)+1}
#' 	  \right],
#' 	\end{align}
#' 	}
#' 	for all \eqn{(u_1, u_2)\in\left[0,1\right]^2}.
#' 	Here
#' 		\deqn{
#' 		t(u_j;a)=\frac{I^{-1}_{\frac{1}{2},{a}}(1-u_j)}{1-I^{-1}_{\frac{1}{2},{a}}(1-u_j)}.
#' 		}
#'
#'
#' - The h-function corresponding to the survival bivariate MGL copula is given by
#'  \deqn{
#'  \begin{align}
#' 	h_{2|1}^{\text{MGL180}}(u_2|u_1;\delta)&=1-h_{2|1}^{\text{MGL}}(1-u_2|1-u_1;\delta),\\
#'  h_{1|2}^{\text{MGL180}}(u_1|u_2;\delta)&=1-h_{1|2}^{\text{MGL}}(1-u_1|1-u_2;\delta).
#'  \end{align}
#'  }
#'
#'
#'
#'
#'
#'
#' @return hcMGL.bivar returns a list with
#' * hfunc1: \eqn{\partial C(u_1,u_2) / \partial u_1,}
#' * hfunc2: \eqn{\partial C(u_1,u_2) / \partial u_2,}
#' @references Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009). Pair-copula constructions of multiple dependence. Insurance: Mathematics and Economics 44 (2), 182-198.
NULL


#' @rdname hBMGL
#' @export
#' @examples
#' hcMGL.bivar(u1 = c(0.1, 0.001, 0.3), u2 = c(0, 0.9999, 0.88), pars = 2)
hcMGL.bivar <- function(u1, u2, pars) {
  a <- 1 / pars
  # another method
  index1 <- u1 == 0 & u2 != 0
  index2 <- u1 != 0 & u2 == 0
  index3 <- u1 == 0 & u2 == 0
  q1 <- qbeta(1 - u1, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u1, shape1 = 0.5, shape2 = a))
  q2 <- qbeta(1 - u2, shape1 = 0.5, shape2 = a) / (1 - qbeta(1 - u2, shape1 = 0.5, shape2 = a))
  z1 <- q2 / (q1 + q2 + 1)
  z2 <- q1 / (q1 + q2 + 1)
  z1[index1] <- 0.0000000001
  z2[index1] <- 1
  z1[index2] <- 1
  z2[index2] <- 0.0000000001
  z1[index3] <- 0.5
  z2[index3] <- 0.5

  hfunc1 <- 1 - pbeta(z1, shape1 = 0.5, shape2 = a + 0.5) # Pr(U_2 <= u_2 | U_1 = u_1) = \partial C(u_1,u_2) / \partial u_1,
  hfunc2 <- 1 - pbeta(z2, shape1 = 0.5, shape2 = a + 0.5) # C(u1|u2)
  list(hfunc1 = hfunc1, hfunc2 = hfunc2)
}

#' @rdname hBMGL
#' @export
#' @examples
#' hcMGL180.bivar(u1 = c(0.1, 0.001, 0.3), u2 = c(0, 0.9999, 0.88), pars = 2)
hcMGL180.bivar <- function(u1, u2, pars) {
  hfunc1 <- 1 - hcMGL.bivar(u1 = 1 - u1, u2 = 1 - u2, pars = pars[1])$hfunc1
  hfunc2 <- 1 - hcMGL.bivar(u1 = 1 - u1, u2 = 1 - u2, pars = pars[1])$hfunc2
  out <- list(hfunc1 = hfunc1, hfunc2 = hfunc2)
  out
}
