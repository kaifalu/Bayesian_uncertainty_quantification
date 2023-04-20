## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("dplyr")
library("tidyr")
library("ggplot2")
library("pscl") 

## ----set_seed, echo=FALSE------------------------------------------------
set.seed(2)

## --------------------------Problem4-------------------------------

## ----define function----------------------------------------------------------------
shortest_interval <- function(df, x.min, x.max, alpha) {
  p <- function(h) {
    g <- function(x) {y <- df(x); ifelse(y > h, y, 0)}
    integrate(g, x.min, x.max)$value - alpha*integrate(df, x.min, x.max)$value
  }
  sample <- seq(x.min, x.max, length.out = 100)
  h = uniroot(p, c(0, max(df(sample))), tol=1e-12)$root # find the threshold of pdf
  g <- function(x) {df(x) - h}
  min_value <- sample[which.max(g(sample))] # find the max value of df(x)-h, should be higher than 0
  interval <- c()
  if (g(x.min) < 0) {interval = append(interval, uniroot(g, c(x.min, min_value), tol=1e-12)$root)}
  else {root = append(root,x.min)}
  if (g(x.max) < 0) {interval = append(interval, uniroot(g, c(min_value, x.max), tol=1e-12)$root)}
  else {interval = append(interval,x.max)}
  interval
}

## ----test function----------------------------------------------------------------

## ----test function 1----------------------------------------------
x.min = 0
x.max = 1
alpha = 0.05
a = 4
b = 8
shortest_interval(function(x) x^(a-1)*(1-x)^(b-1)/exp(lbeta(a, b)), x.min, x.max, 1-alpha)

## ----test function 2----------------------------------------------
x.min = -4
x.max = 1
alpha = 0.05
shortest_interval(function(x) exp(-x^2/2)/sqrt(2*pi), x.min, x.max, 1-alpha)

## ----test function 3----------------------------------------------
x.min = -1
x.max = 1
alpha = 0.05
shortest_interval(function(x) exp(-x^2/2)/sqrt(2*pi), x.min, x.max, 1-alpha)