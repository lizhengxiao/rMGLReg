
#' The bivariate GLMGA distribution
#'
#' @param y1 vector of quantiles.
#' @param y2 vector of quantiles.
#' @param sigma vector of parameters with length = 2.
#' @param a  common parameter with length = 1.
#' @param b vector of parameters with length = 2.
#' @md
#' @details
#'
#' - The multivariate GLMGA density is given by
#' 	\deqn{
#' 	f(y_1,\dots,y_d)
#' 	  =\frac{\Gamma(a+\frac{d}{2})}{\Gamma(a)\Gamma(\frac{1}{2})^d\prod_{j=1}^{d}\sigma_j y_j}
#' 	\frac{\prod_{j=1}^{d}\left[(2b_j)^{\sigma_j}y_j\right]^{-\frac{1}{2\sigma_j}}}
#' 	{\left[\sum_{j=1}^{d}\left((2b_j)^{\sigma_j}y_j\right)^{-\frac{1}{\sigma_j}}+1\right]^{a+\frac{d}{2}}}
#' 	}
#' 	for \eqn{y_{j} >0}, being \eqn{\sigma_j>0, a>0, b_j>0}.
#'
#' 	- Since \eqn{Y_j} (\eqn{j=1,\ldots,d}) are conditionally independent given \eqn{\Theta}, the marginal distributions are obtained  by setting $d=1$ which leads to the densities in
#' 	\code{\link[rMGLReg]{dGLMGA}}.
#'
#' @return Density function of the bivariate GLMGA distribution with parameters sigma, a and b.
#' @export
#'
#' @examples
#' # density function of MGLMGA distribution
#' n.grid <- 20
#' sigma <- c(0.5, 0.5)
#' a <- 20
#' b <- c(5, 5)
#' xgrid <- ygrid <- seq(0.01, 20, length.out = n.grid)
#' grid <- expand.grid("u1" = xgrid, "u2" = ygrid)
#' mtrx3d <- matrix(0, nrow = nrow(grid), ncol = 3)
#' mtrx3d[, 1] <- grid[, 1]
#' mtrx3d[, 2] <- grid[, 2]
#' for (i in 1:nrow(mtrx3d)) {
#'   mtrx3d[i, 3] <- dMGLMGA(
#'     y1 = grid[i, 1], y2 = grid[i, 2],
#'     sigma = sigma,
#'     a = a, b = b
#'   )
#' }
#' head(mtrx3d)
dMGLMGA <- function(y1, y2, sigma, a, b) {
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  b1 <- b[1]
  b2 <- b[2]
  d <- 2
  out1 <- gamma(a + d / 2) / (gamma(a) * gamma(0.5)^d * (sigma1 * y1 * sigma2 * y2)) * ((2 * b1)^sigma1 * y1)^(-1 / (2 * sigma1)) * ((2 * b2)^sigma2 * y2)^(-1 / (2 * sigma2))
  out2 <- (((2 * b1)^(sigma1) * y1)^(-1 / sigma1) + ((2 * b2)^(sigma2) * y2)^(-1 / sigma2) + 1)^(a + d / 2)

  z <- out1 / out2
  return(z)
}
