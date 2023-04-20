## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("reshape2")
library("plyr")
library("ggplot2")

## ----set_seed and data generation--------------------------------------------
set.seed(0)
ntrials = 100
iMax = 4
nMax = 4*4^iMax # nMax = 1024;
U = as.numeric(runif(ntrials*nMax)>0.5) # uniform distribution 0,1
BM = matrix(U, nrow = ntrials) # 100*1024 0-1 matrix

n = 4*4^(1:4) # sample size

## ----scenario 1--------------------------------------------
theta_min = 0
theta_max = 1

## -----Frequentist----
theta_F <- c()
for (i in 1:iMax) {
  for (j in 1:ntrials) {
    theta_F = append(theta_F, sum(BM[j,1:n[i]])/n[i])
  }
}
theta_F <- matrix(theta_F, ncol=iMax, byrow=FALSE)
F_m1 <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(F_m1)
hist(theta_F[,1], prob=TRUE, col = 'grey', breaks = 20, xlab ="MLE", main = "Sampling distribution (n=16)")
lines(density(theta_F[,1]), col = 'blue', lwd=2)
hist(theta_F[,2], prob=TRUE, col = 'grey', breaks = 20, xlab ="MLE", main = "Sampling distribution (n=64)")
lines(density(theta_F[,2]), col = 'blue', lwd=2)
hist(theta_F[,3], prob=TRUE, col = 'grey', breaks = 20, xlab ="MLE", main = "Sampling distribution (n=256)")
lines(density(theta_F[,3]), col = 'blue', lwd=2)
hist(theta_F[,4], prob=TRUE, col = 'grey', breaks = 20, xlab ="MLE", main = "Sampling distribution (n=1024)")
lines(density(theta_F[,4]), col = 'blue', lwd=2)

## -----Bayesian----
theta_B <- c()
for (i in 1:iMax) {
  theta_B = append(theta_B, sum(BM[1,1:n[i]]))
}
B_m1 <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(B_m1)
fun1 <- function(x) {
  h <- function(x) {x^theta_B[1]*(1-x)^(n[1]-theta_B[1])} # used for integration
  (x^theta_B[1]*(1-x)^(n[1]-theta_B[1]))/integrate(h,theta_min,theta_max)$value
}
curve(fun1, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=16)")
fun2 <- function(x) {
  h <- function(x) {x^theta_B[2]*(1-x)^(n[2]-theta_B[2])} # used for integration
  (x^theta_B[2]*(1-x)^(n[2]-theta_B[2]))/integrate(h,theta_min,theta_max)$value
}
curve(fun2, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=64)")
fun3 <- function(x) {
  h <- function(x) {x^theta_B[3]*(1-x)^(n[3]-theta_B[3])} # used for integration
  (x^theta_B[3]*(1-x)^(n[3]-theta_B[3]))/integrate(h,theta_min,theta_max)$value
}
curve(fun3, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=256)")
fun4 <- function(x) {
  h <- function(x) {x^theta_B[4]*(1-x)^(n[4]-theta_B[4])} # used for integration
  (x^theta_B[4]*(1-x)^(n[4]-theta_B[4]))/integrate(h,theta_min,theta_max)$value
}
curve(fun4, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=1024)")


## ----scenario 2--------------------------------------------
theta_min = 0
theta_max = 0.4

## -----Frequentist----
F_m2 <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(F_m2)
hist(theta_F[,1], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0,theta_max), xlab ="MLE", main = "Sampling distribution (n=16)")
lines(density(theta_F[,1]), col = 'blue', lwd=2)
hist(theta_F[,2], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0.1,theta_max), xlab ="MLE", main = "Sampling distribution (n=64)")
lines(density(theta_F[,2]), col = 'blue', lwd=2)
hist(theta_F[,3], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0.3,theta_max), xlab ="MLE", main = "Sampling distribution (n=256)")
lines(density(theta_F[,3]), col = 'blue', lwd=2)
hist(theta_F[,4], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0.3,theta_max), xlab ="MLE", main = "Sampling distribution (n=1024)")
lines(density(theta_F[,4]), col = 'blue', lwd=2)

## -----Bayesian----
B_m2 <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(B_m2)
curve(fun1, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=16)")
curve(fun2, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=64)")
curve(fun3, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=256)")
curve(fun4, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=1024)")


## ----scenario 3--------------------------------------------
theta_min = 0
theta_max = 0.5

## -----Frequentist----
F_m3 <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(F_m3)
hist(theta_F[,1], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0,theta_max), xlab ="MLE", main = "Sampling distribution (n=16)")
lines(density(theta_F[,1]), col = 'blue', lwd=2)
hist(theta_F[,2], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0.2,theta_max), xlab ="MLE", main = "Sampling distribution (n=64)")
lines(density(theta_F[,2]), col = 'blue', lwd=2)
hist(theta_F[,3], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0.3,theta_max), xlab ="MLE", main = "Sampling distribution (n=256)")
lines(density(theta_F[,3]), col = 'blue', lwd=2)
hist(theta_F[,4], prob=TRUE, col = 'grey', breaks = 20, xlim=c(0.4,theta_max), xlab ="MLE", main = "Sampling distribution (n=1024)")
lines(density(theta_F[,4]), col = 'blue', lwd=2)

## -----Bayesian----
B_m3 <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
layout(B_m3)
curve(fun1, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=16)")
curve(fun2, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=64)")
curve(fun3, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=256)")
curve(fun4, theta_min, theta_max, xlab =expression(theta), ylab = expression(paste("p(", theta,"|y)")), main = "Posterior (n=1024)")