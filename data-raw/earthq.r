# earthquakic data set -------------------------
# Number of applicants -------------------------


library(data.table)
library(rMGLReg)
dtnew <- fread('D:/Rpackages/rMGLReg/data-raw/dtnew_copula.csv')
dtnew$ynew <- NULL
dtnew <- dtnew[year %in% seq(1990, 2015)]
dtadjust <- fread('D:/Rpackages/rMGLReg/data-raw/dtcpi.csv', header = T)
dtadjust <- dtadjust[year %in% seq(1990, 2015)]
dtadjust
dtadjust[, CPI_factor := 1/(CPI/615.2)]

setkey(dtadjust, year)
setkey(dtnew, year)
dt <- merge(dtnew, dtadjust)
dt[is.na(y)]$y <- median(dt$y, na.rm = T)
dt[, ynew := y*CPI_factor]

dtnew <- dt[,.(year,
               yraw = y,
               y1 = ynew, magnitude, death, inj,intensity,
               y2 = death + inj)]
dtnew

# marginal distribution for y1
y <- dtnew$y1
Xsigma <- model.matrix(~ 1, data = dtnew);
Xb <- model.matrix(~ 1, data = dtnew)
LLlogmoyalGA3 <- function(y, pars, Xsigma, Xb) {
  sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]])
  b <- exp(Xb %*%pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])])
  a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
  ll <- -0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b)

  loglike <- -sum(ll)
  return(loglike)
}
mlogmoyalGA3 <- optim(fn = LLlogmoyalGA3, Xsigma = Xsigma,
                      Xb = Xb,
                      y = y, hessian = T,
                      control = list(maxit = 50000),
                      #method = 'Nelder-Mead',
                      #par = c(-1,-0.1,0.5,-2,-1))
                      par = c(-1,-1,-1))
# modout(mlogmoyalGA3)
#
pars <- mlogmoyalGA3$par
sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]])
b <- exp(Xb %*%pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])])
a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pGLMGA(y[i], sigma = sigma[i], a = a, b = b[i])
}
resLMGA3 <- qnorm(ufit)
modout(mlogmoyalGA3)
u1 <- ufit
f1 <- c()
for(i in 1:length(y)){
  f1[i] <- dLMGA(y[i], sigma = sigma[i], a = a, b = b[i])
}

