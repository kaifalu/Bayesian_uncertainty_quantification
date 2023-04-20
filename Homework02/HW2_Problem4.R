## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("dplyr")
library("tidyr")
library("ggplot2")
library("pscl") 

## ----set_seed, echo=FALSE------------------------------------------------
set.seed(2)

## ----HW1: Problem4------------------------------------------------------------

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

## ----HW1:test function----------------------------------------------------------------

## ----test function 1----------------------------------------------
x.min = 0
x.max = 1
alpha = 0.05
a = 4
b = 8
exact_SI_1 = shortest_interval(function(x) x^(a-1)*(1-x)^(b-1)/exp(lbeta(a, b)), x.min, x.max, 1-alpha)
exact_SI_1_length = exact_SI_1[2] - exact_SI_1[1]
exact_SI_1
prob_1 = integrate(function(x) x^(a-1)*(1-x)^(b-1)/exp(lbeta(a, b)), x.min, x.max)$value

## ----test function 2----------------------------------------------
x.min = -4
x.max = 1
alpha = 0.05
exact_SI_2 = shortest_interval(function(x) exp(-x^2/2)/sqrt(2*pi), x.min, x.max, 1-alpha)
exact_SI_2_length = exact_SI_2[2] - exact_SI_2[1]
exact_SI_2
prob_2 = integrate(function(x) exp(-x^2/2)/sqrt(2*pi), x.min, x.max)$value

## ----test function 3----------------------------------------------
x.min = -1
x.max = 1
alpha = 0.05
exact_SI_3 = shortest_interval(function(x) exp(-x^2/2)/sqrt(2*pi), x.min, x.max, 1-alpha)
exact_SI_3_length = exact_SI_3[2] - exact_SI_3[1]
exact_SI_3
prob_3 = integrate(function(x) exp(-x^2/2)/sqrt(2*pi), x.min, x.max)$value


## --------------------------HW2:Problem4-------------------------------

## ----define function----------------------------------------------------------------
shortest_interval <- function(thetaV, alpha) {
  # calculate density behind data samples
  prob = density(thetaV)$y/(sum(density(thetaV)$y))
  xx = density(thetaV)$x
  lower = which(cumsum(prob) >= alpha/2)[1] # lower bound
  upper = which(cumsum(prob) >= 1-alpha/2)[1]-1 # upper bound
  if (prob[lower] > prob[upper]) {
    for (i in 1:lower) {
      if ((prob[lower-i] - prob[upper-i]) < 10^-4) {break}
    }
    lower = lower-i+1
    upper = upper-i+1
    dist = sum(prob[lower:upper]) - (1-alpha)
    for (j in 0:i) {
      if (prob[lower+j] > prob[upper-j]) {dist = dist - prob[upper-j]; lower = lower-1}
      else {dist = dist - prob[lower+j]; upper = upper+1}
      if (dist < 0) {break}
    }
    lower = lower+j+1
    upper = upper-j-1
  }
  if (prob[lower] < prob[upper]) {
    for (i in 1:(length(prob)-upper)) {
      if ((prob[lower+i] - prob[upper+i]) > 10^-4) {break}
    }
    lower = lower+i-1
    upper = upper+i-1
    dist = sum(prob[lower:upper]) - (1-alpha)
    for (j in 0:i) {
      if (prob[lower+j] > prob[upper-j]) {dist = dist - prob[upper-j]; lower = lower-1}
      else {dist = dist - prob[lower+j]; upper = upper+1}
      if (dist < 0) {break}
    }
    lower = lower+j+1
    upper = upper-j-1
  }
  c(xx[lower], xx[upper])
}

## ----HW2: test function----------------------------------------------------------------
M <- 50*2^(0:7)
alpha = 0.05

## ----test function 1----------------------------------------------
set.seed(0)
thetaData1 = rbeta(6400, 4, 8)
print('Test Problem 1')
approximate_SI_1 <- c()
coverage_SI_1 <-c()
for (i in 1:8) {
  print(paste("i =",i-1))
  thetaV1 <- thetaData1[1:M[i]]
  SI_1_i = shortest_interval(thetaV1, alpha)
  approximate_SI_1 = append(approximate_SI_1, SI_1_i)
  coverage_SI_1 = append(coverage_SI_1, integrate(function(x) x^(a-1)*(1-x)^(b-1)/exp(lbeta(a, b)), SI_1_i[1], SI_1_i[2])$value/prob_1)
  print(SI_1_i)
}
approximate_SI_1 <- matrix(approximate_SI_1, ncol = 2, byrow = TRUE)
print('Approximate length')
approximate_SI_1[,2] - approximate_SI_1[,1]
print('ratio')
(approximate_SI_1[,2] - approximate_SI_1[,1])/exact_SI_1_length
print('coverage')
coverage_SI_1

## ----test function 2&3----------------------------------------------
set.seed(0)
rr = rnorm(10^4)

## ----test function 2----------------------------------------------
ii = (rr < 1) & (rr > -4)
thetaData2 = rr[ii]
print('Test Problem 2:')
approximate_SI_2 <- c()
coverage_SI_2 <-c()
for (i in 1:8) {
  print(paste("i =",i-1))
  thetaV2 <- thetaData2[1:M[i]]
  SI_2_i = shortest_interval(thetaV2, alpha)
  approximate_SI_2 = append(approximate_SI_2, SI_2_i)
  print(SI_2_i)
  if (SI_2_i[2] > 1) {coverage_SI_2 = append(coverage_SI_2, integrate(function(x) exp(-x^2/2)/sqrt(2*pi), SI_2_i[1], 1)$value/prob_2)}
  else {coverage_SI_2 = append(coverage_SI_2, integrate(function(x) exp(-x^2/2)/sqrt(2*pi), SI_2_i[1], SI_2_i[2])$value/prob_2)}
}
approximate_SI_2 <- matrix(approximate_SI_2, ncol = 2, byrow = TRUE)
print('Approximate length')
approximate_SI_2[,2] - approximate_SI_2[,1]
print('ratio')
(approximate_SI_2[,2] - approximate_SI_2[,1])/exact_SI_2_length
print('coverage')
coverage_SI_2

## ----test function 3----------------------------------------------
iii = (rr < 1) & (rr > -1)
thetaData3 = rr[iii]
print('Test Problem 3:')
approximate_SI_3 <- c()
coverage_SI_3 <-c()
for (i in 1:8) {
  print(paste("i =",i-1))
  thetaV3 <- thetaData3[1:M[i]]
  SI_3_i = shortest_interval(thetaV3, alpha)
  approximate_SI_3 = append(approximate_SI_3, SI_3_i)
  print(SI_3_i)
  if (SI_3_i[1] < -1) {coverage_SI_3 = append(coverage_SI_3, integrate(function(x) exp(-x^2/2)/sqrt(2*pi), -1, SI_3_i[2])$value/prob_3)}
  else {coverage_SI_3 = append(coverage_SI_3, integrate(function(x) exp(-x^2/2)/sqrt(2*pi), SI_3_i[1], SI_3_i[2])$value/prob_3)}
}
approximate_SI_3 <- matrix(approximate_SI_3, ncol = 2, byrow = TRUE)
print('Approximate length')
approximate_SI_3[,2] - approximate_SI_3[,1]
print('ratio')
(approximate_SI_3[,2] - approximate_SI_3[,1])/exact_SI_3_length
print('coverage')
coverage_SI_3
