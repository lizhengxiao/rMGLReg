
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rMGLReg

<!-- badges: start -->
<!-- badges: end -->

The goal of rMGLReg is to

-   provide a nice visualization tool for interpreting the MGL copula
    and MGL-EV copula, along with its survival copulas.

-   show the maximum likelihood (ME) estimation method for the copula
    regression models with/without covariates.

## Installation

You can install the released version of rMGLReg by running the following
lines in R software:

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("lizhengxiao/rMGLReg")
```

## Example - simuation from the proposed copulas.

This is a basic example which shows the normalized scatter plots for
(*δ* = 1.2) with simulated samples of the 3-dimensional MGL copula and
survival MGL copula.

``` r
library(rMGLReg)
## basic example code
set.seed(271)
n <- 1000
delta <- 1.2
d <- 3
U <- rcMGL.multi(n = 1000, d = d, pars = delta)
cor(U, method = "kendall")
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.2322282 0.2335576
#> [2,] 0.2322282 1.0000000 0.2138058
#> [3,] 0.2335576 0.2138058 1.0000000
par(pty = "s")
pairs(U, gap = 0, cex = 0.5)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
set.seed(271)
n <- 1000
delta <- 1.2
d <- 3
U <- rcMGL180.multi(n = 1000, d = d, pars = delta)
cor(U, method = "kendall")
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.2322282 0.2335576
#> [2,] 0.2322282 1.0000000 0.2138058
#> [3,] 0.2335576 0.2138058 1.0000000
par(pty = "s")
pairs(U, gap = 0, cex = 0.5)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## Example - danish fire data set.

As the first example, we fit the bivariate copula and regression models
to the Danish fire insurance data set which was collected from the
Copenhagen Reinsurance Company and comprises 2167 fire losses over the
period 1980-1990.s

``` r
  library(rMGLReg)
  library(fitdistrplus)
#> 载入需要的程辑包：MASS
#> 载入需要的程辑包：survival
  library(splines)
  library(snpar)
  data("danishmulti")
  dt <- data.table::data.table(danishmulti)
  dtnew <- dt[Building>0&Contents>0]
  y1 <- dtnew$Building
  y2 <- dtnew$Contents
  y <- cbind(y1, y2)
  
  # empirical cdf
  u1 <- snpar::kde(y[,1], kernel = "gaus", 
             xgrid = y[,1],
             h = 0.2)$Fhat
  u2 <- snpar::kde(y[,2], kernel = "gaus", 
             xgrid = y[,2],
             h = 0.2)$Fhat
  U <- cbind(u1, u2) # bivariate pseudo copula data.
  plot(U)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

This table reports the estimation results, AIC and BIC values of the
survival MGL and the survival MGL-EV copula, along with four other
families of copulas with positive upper tail indices, the MGB2 copula,
the Gumbel copula, the Student *t* copula, and the Gaussian copula.

``` r
m.norm <- MGL.mle(U = U,
                        copula  = "Normal",
                        initpar = 0.5)

m.t <- MGL.mle(U = U,
                     copula  = "t",
                     initpar = c(0.5, 4))
m.gumbel <- MGL.mle(U = U,
                     copula  = "Gumbel",
                          initpar = c(2))
m.MGLMGA180 <- MGL.mle(U,
                    copula  = "MGL180",
                             initpar = c(1))
m.MGB2 <- MGL.mle(U,
                      copula  = "MGB2",
                        initpar = c(0.1, 2, 0.4))

m.MGLEV180 <- MGL.mle(U,
                         copula  = "MGL-EV180",
                         initpar = c(2))

recap <- function(x){
  res <- c(alpha = x$estimates,
           se = x$se,
           loglike = x$loglike,
           AIC = x$AIC, BIC = x$BIC)
  if(length(res) < 6)
    res <- c(res[1], NA, NA,res[2], NA, NA, res[3:5])
  if (length(res) > 6 & length(res) < 9)
    res <- c(res[1:2], NA, res[3:4], NA, res[5:7])
  res <- as.matrix(res)
  colnames(res) <- x$copula$name
  res}
res.all <- round(cbind(recap(m.norm),
                       recap(m.t),
                       recap(m.gumbel),
                       recap(m.MGLMGA180),
                       recap(m.MGB2),
                       recap(m.MGLEV180)
), 4)
out.com <- t(res.all)
out.com <- out.com[order(out.com[,9], decreasing = T),]
knitr::kable(out.com, digits = 3)
```

|           | alpha |       |       |    se |       |       | loglike |      AIC |      BIC |
|:----------|------:|------:|------:|------:|------:|------:|--------:|---------:|---------:|
| Gaussian  | 0.252 |    NA |    NA | 0.027 |    NA |    NA |  35.602 |  -69.205 |  -63.890 |
| Student   | 0.193 | 3.400 |    NA | 0.032 | 0.457 |    NA |  64.084 | -124.168 | -113.539 |
| Gumbel    | 1.211 |    NA |    NA | 0.022 |    NA |    NA |  79.127 | -156.255 | -150.940 |
| MGL-EV180 | 0.655 |    NA |    NA | 0.040 |    NA |    NA |  81.760 | -161.519 | -156.204 |
| MGL180    | 0.892 |    NA |    NA | 0.067 |    NA |    NA | 115.967 | -229.935 | -224.620 |
| MGB2      | 0.233 | 1.123 | 0.939 | 0.046 | 0.561 | 0.209 | 127.819 | -249.638 | -233.694 |

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!
