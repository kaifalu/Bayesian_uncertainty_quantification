## ----libraries, message=FALSE, warning=FALSE-----------------------------
library(reshape2)
library(plyr)
library(ggplot2)
library(xtable)
library(MASS)
library(mvtnorm)

## ----set_seed------------------------------------------------------------
set.seed(1)

## ----composition_sampler-------------------------------------------------
composition_bivariate_normal = function(n_points, rho) {
  theta = matrix(c(0,0), nrow=n_points, ncol=2, byrow=TRUE)
  v = sqrt(1-rho^2)
  for (i in 1:n_points) {
    theta[i,1] = rnorm(1, 0, 1)
    theta[i,2] = rnorm(1, rho*theta[i,1], v)
  }
  return(theta)
}

theta_1 = composition_bivariate_normal(n<-1000, rho=0.9)
plot(theta_1[,1], theta_1[,2], xlab = expression(X[1]), ylab = expression(X[2]), pch=".", cex=2)

## ----Metropolized composition_sampler-----------------------------------
metropolized_composition_normal = function(theta0, n_points, rho) {
  theta = matrix(c(0,0), nrow=n_points, ncol=2, byrow=TRUE)
  log_q = function(theta) dnorm(theta,log=TRUE)
  current = theta0
  v = sqrt(1-rho^2)
  for (i in 1:n_points) {
    proposed = rnorm(1,current)
    logr = log_q(proposed)-log_q(current)
    current = ifelse(log(runif(1))<logr, proposed, current)
    theta[i,1] = current
    theta[i,2] = rnorm(1, rho*theta[i,1], v)
  }
  return(theta)
}

theta_2 = metropolized_composition_normal(0, n<-1000, rho=0.9)
plot(theta_2[,1], theta_2[,2], xlab = expression(X[1]), ylab = expression(X[2]), pch=".", cex=2)

## ----random walk metropolis hastings algorithm-------------------------
metropolis_hastings_normal = function(theta0, n_points, rho) {
  s = matrix(c(1,rho,rho,1),nrow=2,ncol=2,byrow=TRUE)
  log_q = function(theta) -(theta[1]^2+theta[2]^2-2*rho*theta[1]*theta[2])/(2*(1-rho^2))
  #initialize
  S = diag(2)  # S_0
  current = theta0
  theta = matrix(current, nrow=n_points, ncol=2, byrow=TRUE)
  for (i in 2:n_points) {
    proposed = mvrnorm(1,current,2.4^2*S/2)
    logr = log_q(proposed)-log_q(current)
    if (log(runif(1)) < logr) current = proposed
    theta[i,] = current
    if (i%%50 == 0) S = var(theta[1:i,])
  }
  return(theta)
}

theta_3 = metropolis_hastings_normal(c(0,1), n<-1000, rho=0.9)
plot(theta_3[,1], theta_3[,2], xlab = expression(X[1]), ylab = expression(X[2]), pch=".", cex=2)


## ----bivariate_normal_mcmc, echo=TRUE------------------------------------
gibbs_bivariate_normal = function(theta0, n_points, rho) {
  theta = matrix(theta0, nrow=n_points, ncol=2, byrow=TRUE)
  v = sqrt(1-rho^2)
  for (i in 2:n_points) {
    theta[i,1] = rnorm(1, rho*theta[i-1,2], v)
    theta[i,2] = rnorm(1, rho*theta[i  ,1], v)
  }
  return(theta)
}

theta_4 = gibbs_bivariate_normal(c(0,1), n<-1000, rho=0.9)
plot(theta_4[,1], theta_4[,2], xlab = expression(X[1]), ylab = expression(X[2]), pch=".", cex=2)

## ----problem 1.5------------------------------------
rho_vector = c(0.25,0.5,0.75,0.9,0.95,0.99)
F_m1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)

# sampler_2
layout(F_m1)
for (i in rho_vector) {
  X_sampler2 = metropolized_composition_normal(0, n<-1000, rho<-i)
  acf(X_sampler2[,1], 100, plot=TRUE, main = paste(expression(rho), "=", i))
  print(length(unique(na.omit(X_sampler2[,1]))))
}

#sampler_3
layout(F_m1)
for (i in rho_vector) {
  X_sampler3 = metropolis_hastings_normal(c(0,1), n<-1000, rho<-i)
  acf(X_sampler3[,1], 100, plot=TRUE, main = paste(expression(rho), "=", i))
  print(length(unique(na.omit(X_sampler3[,1]))))
}

#sampler_4
layout(F_m1)
for (i in rho_vector) {
  X_sampler4 = gibbs_bivariate_normal(c(0,1), n<-1000, rho<-i)
  acf(X_sampler4[,1], 100, plot=TRUE, main = paste(expression(rho), "=", i))
  print(length(unique(na.omit(X_sampler4[,1]))))
}


## ----problem 1.6------------------------------------
rho_vector_new = c(0.25,0.5,0.99)
F_m2 <- matrix(c(1,2,3,4,5,6,7,8,9),nrow=3,ncol=3,byrow=TRUE)

# sampler_2
layout(F_m2)
for (i in rho_vector_new) {
  n.out = 101
  xx <- seq(-3, 3, length=n.out)
  grid <- expand.grid(x=xx, y=xx)
  Sigma = matrix(c(1,i,i,1),nrow=2,ncol=2,byrow=TRUE)
  like <- matrix(apply(grid, 1, function(x) mvtnorm::dmvnorm(x,sigma=Sigma)),n.out,n.out)
  contour(xx, xx, like, drawlabels=F, nlevels=10, xlim=c(-3,3), ylim=c(-3,3), 
          xlab=expression(X[1]), ylab=expression(X[2]), main = paste(expression(rho), "=", i)) 
  X_sampler2_new = metropolized_composition_normal(10, n<-1000, rho<-i)
  points(X_sampler2_new[,1], X_sampler2_new[,2], xlab = expression(X[1]), ylab = expression(X[2]), pch=".", cex=2)
  plot(1:n, X_sampler2_new[,1], xlab="Iter", ylab=expression(X[1]), main="Metropolized Composition Sampler")
  plot(1:n, X_sampler2_new[,2], xlab="Iter", ylab=expression(X[2]), main="Metropolized Composition Sampler")
}

# sampler_3
layout(F_m2)
for (i in rho_vector_new) {
  n.out = 101
  xx <- seq(-3, 3, length=n.out)
  grid <- expand.grid(x=xx, y=xx)
  Sigma = matrix(c(1,i,i,1),nrow=2,ncol=2,byrow=TRUE)
  like <- matrix(apply(grid, 1, function(x) mvtnorm::dmvnorm(x,sigma=Sigma)),n.out,n.out)
  contour(xx, xx, like, drawlabels=F, nlevels=10, xlim=c(-3,3), ylim=c(-3,3), 
          xlab=expression(X[1]), ylab=expression(X[2]), main = paste(expression(rho), "=", i)) 
  X_sampler3_new = metropolis_hastings_normal(c(-10,10), n<-1000, rho<-i)
  points(X_sampler3_new[,1], X_sampler3_new[,2], xlab = expression(X[1]), ylab = expression(X[2]), pch=".", cex=2)
  plot(1:n, X_sampler3_new[,1], xlab="Iter", ylab=expression(X[1]), main="Random-Walk Metropolis-Hastings Sampler")
  plot(1:n, X_sampler3_new[,2], xlab="Iter", ylab=expression(X[2]), main="Random-Walk Metropolis-Hastings Sampler")
}

# sampler_4
layout(F_m2)
for (i in rho_vector_new) {
  n.out = 101
  xx <- seq(-3, 3, length=n.out)
  grid <- expand.grid(x=xx, y=xx)
  Sigma = matrix(c(1,i,i,1),nrow=2,ncol=2,byrow=TRUE)
  like <- matrix(apply(grid, 1, function(x) mvtnorm::dmvnorm(x,sigma=Sigma)),n.out,n.out)
  contour(xx, xx, like, drawlabels=F, nlevels=10, xlim=c(-3,3), ylim=c(-3,3), 
          xlab=expression(X[1]), ylab=expression(X[2]), main = paste(expression(rho), "=", i)) 
  X_sampler4_new = gibbs_bivariate_normal(c(-10,10), n<-1000, rho<-i)
  points(X_sampler4_new[,1], X_sampler4_new[,2], xlab = expression(X[1]), ylab = expression(X[2]), pch=".", cex=2)
  plot(1:n, X_sampler4_new[,1], xlab="Iter", ylab=expression(X[1]), main="Gibbs Sampler")
  plot(1:n, X_sampler4_new[,2], xlab="Iter", ylab=expression(X[2]), main="Gibbs Sampler")
}
