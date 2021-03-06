


#' @name earthqCHI
#' @title Chinese earthquake loss data set
#' @description The data set is collected from the ``China earthquake yearbook" .
#' This data set concerns the Chinese mainland, which contains risk information on 291 earthquake events with magnitude greater than 4.0 from 1990 to 2015.
#' @docType data
#' @usage data(earthqCHI)
#' @format earthqCHI contains six columns:
#' \itemize{
#'   \item year: occurrence year.
#'   \item y1: Economic Losses (adjusted). It is expressed in millions of Chinese Yuan (CNY) and are adjusted for inflation to reflect values in 2015.
#'   \item yraw: Economic Losses (original).
#'   \item death: the number of deaths from an earthquake event.
#'   \item inj: the number of injures from an earthquake event.
#'   \item magnitude, intensity : the magnitude and maximum intensity of an earthquake event.
#'   \item y2: casualties are defined as fatalities and injured people, which are due to damage to occupied buildings.
#'   \item u1: cdf of the fitted GLMGA regression model used in Li et al.,(2021) evaluated at y1.
#'   \item u2: cdf of the fitted composite Negative Binomial-Pareto model evaluated at y2.
#'   \item u2_: cdf of the fitted composite Negative Binomial-Pareto model evaluated at y2-1.
#'   \item f1: pdf of the fitted GLMGA regression model used in Li et al.,(2021) evaluated at y1.
#'   \item f2: pdf of the fitted composite Negative Binomial-Pareto model evaluated at y2.
#' }
#' @author Zhengxiao Li, 2021-08-03
#' @references Zhengxiao Li, Jan Beirlant, Shengwang Meng. Generalizing The Log-Moyal Distribution And Regression Models For Heavy-Tailed Loss Data. ASTIN Bulletin: The Journal of the IAA, 11(1):57-99, 2021.
#'
#'
"earthqCHI"



#' #' @name earth_model
#' #' @title Chinese earthquake loss data set
#' #' @description The data set is collected from the ``China earthquake yearbook" .
#' #' This data set concerns the Chinese mainland, which contains risk information on 291 earthquake events with magnitude greater than 4.0 from 1990 to 2015.
#' #' @docType data
#' #' @usage data(earth_model)
#' #' @format earth_model contains six columns:
#' #' \itemize{
#' #'   \item y1: Economic Losses (adjusted). It is expressed in millions of Chinese Yuan (CNY) and are adjusted for inflation to reflect values in 2015.
#' #'   \item yraw: Economic Losses (original)
#' #'   \item death: the number of deaths from an earthquake event
#' #'   \item inj: the number of injures from an earthquake event
#' #'   \item magnitude, intensity : the magnitude and maximum intensity of an earthquake event
#' #'   \item y2: casualties are defined as fatalities and injured people, which are due to damage to occupied buildings.
#' #'   \item u1: 12
#' #'   \item u2: 121
#' #'   \item u2_: 121
#' #'   \item f1: 121
#' #'   \item f2: 121
#' #' }
#' #' @source UCSC Table Browser
#' #' @author Zhengxiao Li, 2021-08-03
#' #'
#' #'
#' #'
#' "earth_model"
