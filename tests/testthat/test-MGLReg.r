#
# # library(testthat)
# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
#
#
# test_that("multiplication works", {
#   data("danishmulti")
#   dt <- data.table::data.table(danishmulti)
#   dt[, year := as.numeric(substr(Date, start = 1, stop = 4))]
#   dtnew <- dt[Building>0&Contents>0]
#   y1 <- dtnew$Building
#   y2 <- dtnew$Contents
#   XY <- cbind(y1, y2)
#
#
#   y <- cbind(y1, y2)
#   m1 <- sROC::kCDF(y[,1], bw = 0.2, xgrid = sort(y[,1]))
#   m2 <- sROC::kCDF(y[,2], bw = 0.2, xgrid = sort(y[,2]))
#   x <- cbind(m1$x, m2$x)
#   u1 <- m1$Fhat;y1 <- y[,1]
#   u1[order(y1)] <- m1$Fhat
#
#   u2 <- m2$Fhat; y2 <- y[,2]
#   u2[order(y2)] <- m2$Fhat
#
#   U <- cbind(u1, u2) # empirical cdf
#   Usample <- U
#
#   U <- Usample
#   dtnew[, x1 := (year - mean(year))/(sd(year))]
#   # dtnew[, Datef := as.Date(Date)]
#   dtnew[, day := difftime(as.Date(Date), as.Date("1980-01-03"), units="days")]
#
#
#   X <- splines::ns(dtnew$year, knots = quantile(dtnew$year, c(0.5)), intercept = T)
#
#   # m.MGL180 <- MGL.reg(U = U, copula = "MGL180",
#   #                                method = "Nelder-Mead",
#   #                                X = X,
#   #                                control = list(maxit = 100000,
#   #                                               fnscale = -1),
#   #                                #initpar = c(-0.32, 0.001, 0.001, 0.001),
#   #                                initpar = c(-0.32, 0.001, 0.001)
#   # )
#   #
#   m.MGLEV180 <- MGL.reg(U = U, copula = "MGL-EV180",
#                          method = "Nelder-Mead",
#                          X = X,
#                          control = list(maxit = 100000,
#                                         fnscale = -1),
#                          #initpar = c(-0.32, 0.001, 0.001, 0.001),
#                          initpar = c(-0.32, 0.001, 0.001)
#   )
# })
