
#'
#' Bivarite MGL-EV copula
#'
#' @param u_1,u_2 numeric vectors of equal length with values in [0,1].
#' @param param copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @return Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGLEV.bivar(u1 = 0.001, u2 = 0.999, param = 1)
#' pcMGLEV.bivar(u1 = c(0.3, 0.9), u2 = c(0.5, 0.8), param = 2)
pcMGLEV.bivar <- function(u1, u2, param){
  out <- u1 + u2 - 1 + pcMGLEV180.bivar(1 - u1, 1 - u2, param = param)
  return(out)
}


#'
#' Bivarite MGL-EV copula
#'
#' @param u_1,u_2 numeric vectors of equal length with values in [0,1].
#' @param param copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @return Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGLEV.bivar(u1 = 0.001, u2 = 0.999, param = 1)
#' pcMGLEV.bivar(u1 = c(0.3, 0.9), u2 = c(0.5, 0.8), param = 2)
dcMGLEV.bivar <- function(u1, u2, param){
  dcMGLEV180.bivar(1 - u1, 1 - u2, param = param)
}



#'
#' Bivarite surivival MGL-EV copula
#'
#' @param u_1,u_2 numeric vectors of equal length with values in [0,1].
#' @param param copula parameter, denoted by delta
#' @param log logical; if TRUE, probabilities/densities p are returned as log(p).
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#' @return Density, distribution function, and random generation for the d-dimensional MGL copula with copula parameter delta.
#' @export
#'
#' @examples
#' dcMGLEV180.bivar(u1 = 0.5, u2 = 0.78, param = 1.8)
#' pcMGLEV180.bivar(u1 = c(0.3, 0.9), u2 = c(0.5, 0.8), param = 2)
pcMGLEV180.bivar <- function(u1, u2, param){
  lf <- function(y1, y2, param){
    delta <- param
    z2 <- (y2)^(- delta)/(y1^(- delta) + y2^(- delta))
    z1 <- (y1)^(- delta)/(y1^(- delta) + y2^(- delta))
    out <- y1*pbeta(z2, shape1 = 0.5, shape2 = 1/delta + 0.5) + y2*pbeta(z1, shape1 = 0.5, shape2 = 1/delta + 0.5)
    return(out)
  }
  # lf <- Vectorize(lf)
  out <- exp(-lf(-log(u1), -log(u2), param = param))
  index1 <- u1 == 0|u2 == 0
  index2 <- u1==1&u2!=0
  index3 <- u1==0&u2==1
  out[index1] = 0
  out[index2] = u2[index2]
  out[index3] = u1[index3]

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


dcMGLEV180.bivar <- function(u1, u2, param){
  lf <- function(y1, y2, param){
    delta <- param
    z2 <- (y2)^(- delta)/(y1^(- delta) + y2^(- delta))
    z1 <- (y1)^(- delta)/(y1^(- delta) + y2^(- delta))
    out <- y1*pbeta(z2, shape1 = 0.5, shape2 = 1/delta + 0.5) + y2*pbeta(z1, shape1 = 0.5, shape2 = 1/delta + 0.5)
    return(out)
  }

  pbeta_dev <- function(y, shape1, shape2){
    (1 - y)^(shape2 - 1)*y^(shape1 - 1)/beta(shape1, shape2)
  }

  pbeta_dev2 <- function(y, shape1, shape2){
    -(1 - y)^(shape2-2)*y^(shape1-2)*((y-1)*shape1 + (shape2-2)*y+1)/beta(shape1, shape2)
  }
  z2_x1 <- function(x1, x2, param){
    delta <- param
    part1 <- delta*x1^(-delta-1)*x2^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1/part2
  }

  z1_x1 <- function(x1, x2, param){
    delta <- param
    part1 <- -delta*x1^(-delta-1)*x2^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1/part2
  }

  z2_x2 <- function(x1, x2, param){
    delta <- param
    part1 <- -delta*x2^(-delta-1)*x1^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1/part2
  }

  z1_x2 <- function(x1, x2, param){
    delta <- param
    part1 <- delta*x2^(-delta-1)*x1^(-delta)
    part2 <- (x2^(-delta) + x1^(-delta))^2
    part1/part2
  }

  z2_x1x2 <- function(x1, x2, param){
    delta <- param
    part1 <- delta^2*(x1*x2)^(-delta-1)*(x2^(-delta) - x1^(-delta))
    part2 <- (x2^(-delta) + x1^(-delta))^3
    part1/part2
  }

  z1_x1x2 <- function(x1, x2, param){
    delta <- param
    part1 <- delta^2*(x1*x2)^(-delta-1)*(x1^(-delta) - x2^(-delta))
    part2 <- (x2^(-delta) + x1^(-delta))^3
    part1/part2
  }



  lf_dev_x1 <- function(x1, x2, param){
    delta <- param
    z2 <- x2^(-delta)/(x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    out <- pbeta(z2, shape1 = 0.5, shape2 = 1/delta + 0.5) + x1*pbeta_dev(z2, shape1 = 0.5, shape2 = 1/delta + 0.5)*z2_x1(x1 = x1, x2 = x2, param = delta) + x2*pbeta_dev(z1, shape1 = 0.5, shape2 = 1/delta + 0.5)*z1_x1(x1 = x1, x2 = x2, param = delta)
    out
  }


  lf_dev_x2 <- function(x1, x2, param){
    delta <- param
    z2 <- x2^(-delta)/(x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    out <- pbeta(z1, shape1 = 0.5, shape2 = 1/delta + 0.5) + x2*pbeta_dev(z1, shape1 = 0.5, shape2 = 1/delta + 0.5)*z1_x2(x1 = x1, x2 = x2, param = delta) + x1*pbeta_dev(z2, shape1 = 0.5, shape2 = 1/delta + 0.5)*z2_x2(x1 = x1, x2 = x2, param = delta)
    out
  }

  lf_dev2 <- function(x1, x2, param){
    delta <- param
    z2 <- x2^(-delta)/(x2^(-delta) + x1^(-delta))
    z1 <- 1 - z2
    If <- function(z){pbeta(z, shape1 = 0.5, shape2 = 1/delta + 0.5)}
    If_dev <- function(z){pbeta_dev(y = z, shape1 = 0.5, shape2 = 1/delta + 0.5)}
    If_dev2 <- function(z){pbeta_dev2(y = z, shape1 = 0.5, shape2 = 1/delta + 0.5)}
    p1 <- If_dev(z2)*z2_x2(x1 = x1, x2 = x2, param = delta)
    p2 <- x1*(If_dev2(z2)*z2_x2(x1 = x1, x2 = x2, param = delta)*z2_x1(x1 = x1, x2 = x2, param = delta)
              + If_dev(z2)*z2_x1x2(x1 = x1, x2 = x2, param = delta))
    p3 <-  ((If_dev(z1) + x2*If_dev2(z1)*z1_x2(x1 = x1, x2 = x2, param = delta))*z1_x1(x1 = x1, x2 = x2, param = delta)
            +  x2*If_dev(z1)*z1_x1x2(x1 = x1, x2 = x2, param = delta))
    out <- p1 + p2 + p3
    out
  }
  # final results
  p1 <- pcMGLEV180.bivar(u1 = u1, u2 = u2, param = param)/(u1*u2)
  # p1 <- exp(-lf(y1 = -log(u1), y2 = -log(u2), param = param))/(u1*u2)
  p2 <- lf_dev_x1(-log(u1), -log(u2), param = param)*lf_dev_x2(-log(u1), -log(u2), param = param) - lf_dev2(-log(u1), -log(u2), param = param)
  out <- p1*p2
  return(out)
}


