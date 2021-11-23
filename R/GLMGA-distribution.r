
#' @name GLMGA
#' @rdname  GLMGA
#' @title The GLMGA distribution
#'
#' @param y vector of quantiles.
#' @param u vector of probabilities.
#' @param sigma parameter of GLMGA distribution.
#' @param a parameter of GLMGA distribution.
#' @param b parameter of GLMGA distribution.
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @importFrom stats runif
#' @importFrom stats pbeta
#' @description Density (\code{dGLMGA}), distribution function (\code{pGLMGA}), quantile function (\code{qGLMGA}) and random generation (\code{rGLMGA}) for the GLMGA distribution with parameters sigma, a and b.
#' @references Zhengxiao Li, Jan Beirlant, Shengwang Meng. Generalizing The Log-Moyal Distribution And Regression Models For Heavy-Tailed Loss Data. ASTIN Bulletin: The Journal of the IAA, 11(1):57-99, 2021.
#' @details
#' The GLMGA distribution with parameters (sigma, a, b) has density
#' \deqn{f(y)=\frac{(2b)^a}{\sigma B(a,\frac{1}{2})}\frac{y^{-(\frac{1}{2\sigma}+1)}}{(y^{-\frac{1}{\sigma}}+2b)^{a + \frac{1}{2}}},}
#' for \eqn{y>0,\sigma>0, a>0, b>0}.
#'
#' The cumulative distribution function \eqn{F(y)} is
#'
#' \deqn{F(y)=1-I_{\frac{1}{2},a}(\frac{y^{-1/\sigma}}{y^{-1/\sigma}+2b}).}
#'
#' Here \eqn{I_{m,n}()} is the beta cumulative distribution function (or regularized incomplete beta function) with parameters shape1 = m and shape2 = n
#' implemented by R's \code{\link[stats]{pbeta}} and defined in its help.
#'
#' The quantile function \eqn{F^{-1}(u)} is
#'
#' \deqn{(2b)^{-\sigma}[\frac{I^{-1}_{\frac{1}{2},a}(1-p)}{1-I^{-1}_{\frac{1}{2},a}(1-p)}]^{-\sigma},}
#' where \eqn{u \in (0,1)}, and \eqn{I_{m,n}^{-1}()} denotes the inverse of the beta cumulative distribution function (or regularized incomplete beta function)
#' with parameters shape1 = m and shape2 = n
#' implemented by R's \code{\link[stats]{qbeta}}.
#'
#'
#'
#' @return
#' \code{dGLMGA} gives the density, \code{pGLMGA} gives the distribution function, \code{qGLMGA} gives the quantile function, and \code{rGLMGA} generates random deviates.
#'
#' Invalid arguments will result in return value NaN, with a warning.
#'
#' The length of the result is determined by n for rgamma, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#'
NULL




#' @rdname  GLMGA
#' @export
#' @examples
#' # density function at value 0.5 and 0.1
#' dGLMGA(c(0.5, 0.1), sigma = 2, a = 2, b = 3, log = FALSE)
dGLMGA <- function(y, sigma, a, b, log = FALSE) {
  if (log == FALSE) {
    exp(-0.5 * log(2 * pi) - log(sigma) + a * log(b) + lgamma(a + 0.5) - lgamma(a) - (1 / (2 * sigma) + 1) * log(y) - (a + 0.5) * log(0.5 * (1 / y)^(1 / sigma) + b))
  } else if (log == TRUE) {
    -0.5 * log(2 * pi) - log(sigma) + a * log(b) + lgamma(a + 0.5) - lgamma(a) - (1 / (2 * sigma) + 1) * log(y) - (a + 0.5) * log(0.5 * (1 / y)^(1 / sigma) + b)
  }
}


#' @rdname GLMGA
#' @export
#'
#' @examples
#' # cdf at value 10 and 20.
#' pGLMGA(c(10, 20), sigma = 2, a = 2, b = 3)
pGLMGA <- function(y, sigma, a, b) {
  z <- y^(-1 / sigma) / (y^(-1 / sigma) + 2 * b)
  p <- 1 - pbeta(z, shape1 = 0.5, shape2 = a)
  p
}


#' @rdname GLMGA
#' @export
#'
#' @examples
#' # quantile function at level 50% and 10%
#' qGLMGA(c(0.5, 0.1), sigma = 2, a = 2, b = 3)
qGLMGA <- function(u, sigma, a, b) {
  c <- (2 * b)^(-sigma)
  # I <- pbeta(u, shape1 = 0.5, shape2 = a)
  Iinv <- qbeta(1 - u, shape1 = 0.5, shape2 = a)
  c * (Iinv / (1 - Iinv))^(-sigma)
}


#' @rdname GLMGA
#' @export
#' @examples
#' # simulate 10 samples from GLMGA distribution with parameters (2, 2, 3)
#' rGLMGA(n = 10, sigma = 2, a = 2, b = 3)
rGLMGA <- function(n, sigma, a, b) {
  u <- runif(n, min = 0, max = 1)
  qGLMGA <- Vectorize(qGLMGA)
  r <- qGLMGA(u, sigma = sigma, a = a, b = b)
  r
}


#' #' Sum of vector elements
#' #'
#' #' \code{sum} returns the sum of all the values present in its arguments.
#' #'
#' #' This is a generic function: methods can be defined for it directly
#' #' or via the \code{\link{Summary}} group generic. For this to work properly,
#' #' the arguments \code{...} should be unnamed, and dispatch is on the
#' #' first argument.
#' #'
#' #' @param ... Numeric, complex, or logical vectors.
#' #' @param na.rm A logical scalar. Should missing values (including NaN)
#' #'   be removed?
#' #' @return If all inputs are integer and logical, then the output
#' #'   will be an integer. If integer overflow
#' #'   \url{https://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#' #'   will be NA with a warning. Otherwise it will be a length-one numeric or
#' #'   complex vector.
#' #'
#' #'   Zero-length vectors have sum 0 by definition. See
#' #'   \url{https://en.wikipedia.org/wiki/Empty_sum} for more details.
#' #' @examples
#' #' sum(1:10)
#' #' sum(1:5, 6:10)
#' #' sum(F, F, F, T, T)
#' #'
#' #' sum(.Machine$integer.max, 1L)
#' #' sum(.Machine$integer.max, 1)
#' #'
#' #' \dontrun{
#' #' sum("a")
#' #' }
#' sum <- function(..., na.rm = TRUE) {}
