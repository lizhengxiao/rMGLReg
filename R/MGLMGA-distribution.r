
#' The bivariate GLMGA distribution
#'
#' @param y1 vector of quantiles.
#' @param y2 vector of quantiles.
#' @param sigma vector of parameters
#' @param a  common parameter
#' @param b vector of parameters
#'
#' @return Density, for the GLMGA distribution with parameters sigma, a and b.
#' @export
#'
#' @examples
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
