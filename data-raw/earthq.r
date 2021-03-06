# earthquakic data set -------------------------
# Number of applicants -------------------------


library(data.table)
library(rMGLReg)
devtools::load_all()
dtnew <- fread("D:/Rpackages/rMGLReg/data-raw/dtnew_copula.csv")
dtnew$ynew <- NULL
dtnew <- dtnew[year %in% seq(1990, 2015)]
dtadjust <- fread("D:/Rpackages/rMGLReg/data-raw/dtcpi.csv", header = T)
dtadjust <- dtadjust[year %in% seq(1990, 2015)]
dtadjust
dtadjust[, CPI_factor := 1 / (CPI / 615.2)]

setkey(dtadjust, year)
setkey(dtnew, year)
dt <- merge(dtnew, dtadjust)
dt[is.na(y)]$y <- median(dt$y, na.rm = T)
dt[, ynew := y * CPI_factor]

dtnew <- dt[, .(year,
  yraw = y,
  y1 = ynew, magnitude, death, inj, intensity,
  y2 = death + inj
)]
dtnew

# marginal distribution for y1
y <- dtnew$y1
Xsigma <- model.matrix(~1, data = dtnew)
Xb <- model.matrix(~1, data = dtnew)
LLlogmoyalGA3 <- function(y, pars, Xsigma, Xb) {
  sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]])
  b <- exp(Xb %*% pars[(dim(Xsigma)[2] + 1):(dim(Xsigma)[2] + dim(Xb)[2])])
  a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2] + 1])
  ll <- -0.5 * log(2 * pi) - log(sigma) + a * log(b) + lgamma(a + 0.5) - lgamma(a) - (1 / (2 * sigma) + 1) * log(y) - (a + 0.5) * log(0.5 * (1 / y)^(1 / sigma) + b)

  loglike <- -sum(ll)
  return(loglike)
}
mlogmoyalGA3 <- optim(
  fn = LLlogmoyalGA3, Xsigma = Xsigma,
  Xb = Xb,
  y = y, hessian = T,
  control = list(maxit = 50000),
  # method = 'Nelder-Mead',
  # par = c(-1,-0.1,0.5,-2,-1))
  par = c(-1, -1, -1)
)
# modout(mlogmoyalGA3)
#
pars <- mlogmoyalGA3$par
sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]])
b <- exp(Xb %*% pars[(dim(Xsigma)[2] + 1):(dim(Xsigma)[2] + dim(Xb)[2])])
a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2] + 1])
ufit <- c()
for (i in 1:length(y)) {
  ufit[i] <- pGLMGA(y[i], sigma = sigma[i], a = a, b = b[i])
}
resLMGA3 <- qnorm(ufit)
u1 <- ufit
f1 <- c()
for (i in 1:length(y)) {
  f1[i] <- dGLMGA(y[i], sigma = sigma[i], a = a, b = b[i])
}


## ----message=TRUE, warning=FALSE------------------------
library(gamlss)
library(gamlss.tr)
library(evd)
library(extRemes)
library(ggplot2)
library(qqplotr)
library(patchwork)
umin <- 20
gen.trun(par = umin + 1, type = "right", name = "tr", family = "NBI")
gen.trun(par = umin + 1, type = "right", name = "tr", family = "PO")
gen.trun(par = umin + 1, type = "right", name = "tr", family = "ZINBI")
gen.trun(par = umin + 1, type = "right", name = "tr", family = "ZIP")
dt2 <- dtnew[y2 <= umin]
m.NB <- gamlss(y2 ~ 1, family = NBItr, data = dt2, method = mixed(2, 10000))
summary(m.NB)
pars <- c(coef(m.NB, what = "mu"), coef(m.NB, what = "sigma"))
summary1 <- summary(m.NB)
summary(m.NB)
logLik(m.NB)
mu.NB <- exp(coef(m.NB))
sigma.NB <- exp(coefficients(m.NB, what = "sigma"))


newtheme <- theme_bw() + theme(
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(margin = margin(t = 0.25, unit = "cm")),
  axis.text.y = element_text(
    margin = margin(r = 0.25, unit = "cm"),
    size = 10,
    vjust = 0.5,
    hjust = 0.5
  ),
  axis.ticks.length = unit(-0.1, "cm"),
  plot.title = element_text(hjust = 0.5),
  legend.direction = "vertical",
  legend.position = c(.95, .99),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.text = element_text(size = 10)
)

# par(mfrow = c(1, 2), tcl = 0.3, mgp = c(1.5,0,0))
# qqnorm(res.death, ylim = c(-3,3), xlim = c(-3,3), main = '')
# abline(0,1, col = 'red')

res.death <- residuals(m.NB, what = "z-scores")
df <- data.table(y = res.death)
v1 <- ggplot(df, aes(sample = y)) +
  stat_qq_point() +
  # stat_qq_band(bandType = "ks", fill = "#8DA0CB", alpha = 0.4, conf = 0.95) +
  stat_qq_line(col = "red", identity = TRUE) +
  newtheme +
  labs(title = "Causalities below the threshold", x = "Standard Normal Quantiles", y = "Sample Quantiles") +
  scale_x_continuous(
    limits = c(-3, 3),
    breaks = seq(-3, 3, by = 1)
  ) +
  scale_y_continuous(
    limits = c(-3, 3),
    breaks = seq(-3, 3, by = 1)
  ) +
  geom_abline(slope = 1, col = "red", lwd = 1)

# v1
# =========================================================================
# death > umin
# ==========================================================================
X1 <- model.matrix(~1, data = dtnew[y2 > umin])
X2 <- model.matrix(~1, data = dtnew[y2 > umin])
y2new <- dtnew[y2 > umin]$y2
GPD.MLE <- fevd(
  x = y2,
  scale.fun = ~1,
  shape.fun = ~1,
  threshold = umin, method = "MLE", use.phi = F,
  type = "GP", data = dtnew
)
GPD.MLE
# plot(GPD.MLE)
GPD.par <- GPD.MLE$results$par
pars <- GPD.MLE$results$par
GPD.hessian <- GPD.MLE$results$hessian
GPD.se <- sqrt(diag(solve(GPD.hessian)))
GPD.Z <- GPD.par / GPD.se ##  Z????????????
GPD.p <- ifelse(GPD.Z >= 0, pnorm(GPD.Z, lower = F) * 2, pnorm(GPD.Z) * 2) ## p???
cbind(GPD.par, GPD.se, GPD.p)

scale <- (X1 %*% pars[1:dim(X1)[2]])
shape <- X2 %*% pars[c((dim(X1)[2] + 1):(dim(X1)[2] + dim(X2)[2]))]
ufit <- 0
for (i in 1:length(y2new)) {
  ufit[i] <- evd::pgpd(y2new[i], loc = umin, scale = scale[i], shape = shape[i])
}
resGPD.death <- qnorm(ufit)
# qqnorm(resGPD.death, ylim = c(-3, 3), xlim = c(-3, 3),
#        main = '')
# abline(0, 1, col = 'red')

df <- data.table(y = resGPD.death)
v2 <- ggplot(df, aes(sample = y)) +
  stat_qq_point() +
  # stat_qq_band(bandType = "ks", fill = "#8DA0CB", alpha = 0.4, conf = 0.95) +
  stat_qq_line(col = "red", identity = TRUE) +
  newtheme +
  labs(title = "Causalities above the threshold", x = "Standard Normal Quantiles", y = "Sample Quantiles") +
  scale_x_continuous(
    limits = c(-3, 3),
    breaks = seq(-3, 3, by = 1)
  ) +
  scale_y_continuous(
    limits = c(-3, 3),
    breaks = seq(-3, 3, by = 1)
  ) +
  geom_abline(slope = 1, col = "red", lwd = 1)
# v2

v0 <- v1 + v2 + plot_layout(ncol = 2)
v0


# ======================================================================================
# GPD and negative binom
y <- dtnew$y2
X0 <- model.matrix(~1, data = dtnew)
X1 <- model.matrix(~1, data = dtnew) # GP
X2 <- model.matrix(~1, data = dtnew) # GP
par1 <- GPD.MLE$results$par
par2 <- c(coef(m.NB, what = "mu"), coef(m.NB, what = "sigma"))
p1 <- nrow(dt2) / nrow(dtnew)
scale <- X1 %*% par1[1:dim(X1)[2]]
shape <- X2 %*% par1[c((dim(X1)[2] + 1):(dim(X1)[2] + dim(X2)[2]))]
dtu2 <- data.table(
  mu = as.vector(exp(X0 %*% par2[1:dim(X0)[2]])),
  sigma = exp(par2[length(par2)]),
  scale = as.vector(X1 %*% par1[1:dim(X1)[2]]),
  shape = as.vector(X2 %*% par1[c((dim(X1)[2] + 1):(dim(X1)[2] + dim(X2)[2]))]),
  y = y, y_ = y - 1,
  y__ = y - 2
)
dtu2[y_ < 0]$y_ <- 0
dtu2[y__ < 0]$y__ <- 0
dtu2
u <- with(dtu2, {
  u2 <- u2_ <- u2__ <- u0 <- 0
  for (i in 1:length(y)) {
    if (y[i] <= umin) {
      u2[i] <- p1 * as.vector(pNBItr(y[i], mu = mu[i], sigma = sigma[i]))
    } else {
      u2[i] <- p1 + (1 - p1) * evd::pgpd(y[i], loc = umin, scale = scale[i], shape = shape[i])
    }

    if (y_[i] <= umin) {
      u2_[i] <- p1 * as.vector(pNBItr(y_[i], mu = mu[i], sigma = sigma[i]))
    } else {
      u2_[i] <- p1 + (1 - p1) * evd::pgpd(y_[i], loc = umin, scale = scale[i], shape = shape[i])
    }

    if (y__[i] <= umin) {
      u2__[i] <- p1 * as.vector(pNBItr(y__[i], mu = mu[i], sigma = sigma[i]))
    } else {
      u2__[i] <- p1 + (1 - p1) * evd::pgpd(y__[i], loc = umin, scale = scale[i], shape = shape[i])
    }
    u0 <- pNBItr(umin, mu = mu[i], sigma = sigma[i])
  }
  u <- cbind(u2, u2_, u2__, u0)
  u
})
u2 <- u[, 1]
u2_ <- u[, 2]
u2__ <- u[, 3]
u2_[which(dtnew$y2 - 1 < 0)] <- 0
u2__[which(dtnew$y2 - 2 < 0)] <- 0
dtu2$u2 <- u2
dtu2$u2_ <- u2_
dtu2$u2__ <- u2__
dtu2

# =========================================================================================
# density of y2 --------------------------------------------------------------------------
dtu2$f2 <- with(dtu2, {
  f <- 0
  for (i in 1:length(y)) {
    if (y[i] <= umin) {
      f[i] <- p1 * as.vector(dNBItr(y[i], mu = mu[i], sigma = sigma[i]))
    } else {
      f[i] <- (1 - p1) * evd::dgpd(y[i], loc = umin, scale = scale[i], shape = shape[i])
    }
  }
  f
})
f2 <- dtu2$f2



dtnew$f1 <- f1
dtnew$f2 <- dtu2$f2
dtnew$u1 <- u1
dtnew$u2 <- dtu2$u2
dtnew$u2_ <- dtu2$u2_




# plot
dtnew[, magnitude.f := cut(magnitude, c(3.9, 5.2, 6, 9))]
dtnew[, .(corc = cor(u1, u2)), by = magnitude.f]
summary(dtnew$magnitude)
table(dtnew$magnitude.f)
library(latex2exp)
df <- dtnew
m1 <- ggplot(df, aes(
  x = u1, y = u2, group = magnitude.f,
  col = magnitude.f,
  shape = magnitude.f
)) +
  geom_point() +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(margin = margin(t = 0.25, unit = "cm")),
    axis.text.y = element_text(
      margin = margin(r = 0.25, unit = "cm"),
      size = 10,
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.ticks.length = unit(-0.1, "cm"),
    plot.title = element_text(hjust = 0.5),
    # legend.direction = 'vertical',
    # legend.position = c(.95, .99),
    # legend.justification = c("right", "top"),
    # legend.box.just = "right",
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "Scatterplot of Losses-Casualties",
    x = TeX("$u_1$"),
    y = TeX("$u_2$")
  ) +
  # scale_colour_discrete(breaks = c(1,2,3)) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(col = "magnitude", shape = "magnitude")

m1


# write.csv(dtnew, "data-raw/earth-model.csv")
earthqCHI <- dtnew
earthqCHI$magnitude.f <- NULL
usethis::use_data(earthqCHI, overwrite = TRUE)
