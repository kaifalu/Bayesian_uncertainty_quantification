## ----libraries, message=FALSE, warning=FALSE-----------------------------
library(reshape2)
library(plyr)
library(ggplot2)
library(xtable)
library(MASS)
library(mvtnorm)
library(BNPdensity)
library(mcmc)


data(logit) ## logit: 100*5 data array


## ------ Calculate log-posterior ----------------------------------
lupost_factory <- function(x, y, beta, sigma = 2) {
  eta <- as.numeric(x %*% beta)
  logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
  logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
  logl <- sum(logp[y == 1]) + sum(logq[y == 0])  # binomial(u), just keep y==1, p=u, q=1-u
  ## log: makes multiply become sum 
  return(logl + log(1/16) - sum(beta^2) / (2 * sigma^2)) # log-posterior
}

## ------ Calculate likelihood ----------------------------------
likelihood_factory <- function(x, y, beta) {
  eta <- as.numeric(x %*% beta)
  logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
  logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
  logl <- sum(logp[y == 1]) + sum(logq[y == 0])  # binomial(u), just keep y==1, p=u, q=1-u
  ## log: makes multiply become sum 
  return(logl)
}

## ----random walk metropolis hastings algorithm-------------------------
#sigma = 2
metropolis_hastings_normal = function(x,y,beta, n_points, scale) {
  #initialize
  S = diag(length(beta)) * scale^2  #
  current = beta
  coef = matrix(current, nrow=n_points, ncol=length(beta), byrow=TRUE)
  lupost = matrix(lupost_factory(x,y,beta), nrow=n_points, ncol=1, byrow=TRUE)
  likelihood = matrix(likelihood_factory(x,y,beta), nrow=n_points, ncol=1, byrow=TRUE)
  nn = 0
  #index = sample.int(length(y),blen)
  for (i in 2:n_points) {
    proposed = mvrnorm(1,current,2.4^2*S/length(beta))
    logr = lupost_factory(x,y,proposed,sigma=2)-lupost_factory(x,y,current,sigma=2)
    if (log(runif(1)) < logr) {current = proposed; nn = nn + 1}
    coef[i,] = current
    lupost[i] = lupost_factory(x,y,coef[i,],sigma=2)
    likelihood[i] = likelihood_factory(x,y,coef[i,])
    # if (i%%50 == 0) S = var(theta[1:i,])
  }
  sampler = list(coefficient=coef, log_posterior=lupost, acceptance=nn/n_points,
                 coef_means = colMeans(coef), # MCMC estimate
                 coef_variance = colMeans(coef^2)-colMeans(coef)^2,
                 likelihood_est = mean(likelihood),
                 likelihood_CI = c(sort(likelihood)[round(n_points*0.025)],sort(likelihood)[round(n_points*0.975)]))
  return(sampler)
}

# sigma comes from half-cauchy
metropolis_hastings_cauchy = function(x,y,beta, n_points, scale) {
  #initialize
  current = beta
  coef = matrix(current, nrow=n_points, ncol=length(beta), byrow=TRUE)
  lupost = matrix(lupost_factory(x,y,beta,2), nrow=n_points, ncol=1, byrow=TRUE)
  likelihood = matrix(likelihood_factory(x,y,beta), nrow=n_points, ncol=1, byrow=TRUE)
  nn = 0
  for (i in 2:n_points) {
    sigma = abs(rcauchy(1, 0, 1))
    S = diag(length(beta)) * scale^2
    proposed = mvrnorm(1,current,2.4^2*S/length(beta))
    logr = lupost_factory(x,y,proposed,sigma)-lupost_factory(x,y,current,sigma)
    if (log(runif(1)) < logr) {current = proposed; nn = nn + 1}
    coef[i,] = current
    lupost[i] = lupost_factory(x,y,coef[i,],sigma)
    likelihood[i] = likelihood_factory(x,y,coef[i,])
  }
  sampler = list(coefficient=coef, log_posterior=lupost, acceptance=nn/n_points,
                 coef_means = colMeans(coef), # MCMC estimate
                 coef_variance = colMeans(coef^2)-colMeans(coef)^2,
                 likelihood_est = mean(likelihood),
                 likelihood_CI = c(sort(likelihood)[round(n_points*0.025)],sort(likelihood)[round(n_points*0.975)]))
  return(sampler)
}

## -----------------------RWMH Importance Sampling-----------------------------
# sigma = 2
RWMH_importance_sampling_normal = function(x,y,beta, n_points) {
  #initialize
  current = beta
  S = diag(length(beta))
  coef = matrix(current, nrow=n_points, ncol=length(beta), byrow=TRUE)
  lupost = matrix(lupost_factory(x,y,beta,sigma=2), nrow=n_points, ncol=1, byrow=TRUE)
  likelihood = matrix(likelihood_factory(x,y,beta), nrow=n_points, ncol=1, byrow=TRUE)
  nn = 0
  for (i in 2:n_points) {
    proposed = mvrnorm(1,coef[which.max(lupost),],2.4^2*S/length(beta))
    logr = lupost_factory(x,y,proposed,sigma=2)-lupost_factory(x,y,current,sigma=2)
    if (log(runif(1)) < logr) {current = proposed; nn = nn + 1}
    coef[i,] = current
    lupost[i] = lupost_factory(x,y,coef[i,],sigma=2)
    likelihood[i] = likelihood_factory(x,y,coef[i,])
    if (i%%50 == 0) S = var(coef[1:i,])
  }
  sampler = list(coefficient=coef, log_posterior=lupost, acceptance=nn/n_points,
                 coef_means = colMeans(coef), # MCMC estimate
                 coef_variance = colMeans(coef^2)-colMeans(coef)^2,
                 likelihood_est = mean(likelihood),
                 likelihood_CI = c(sort(likelihood)[round(n_points*0.025)],sort(likelihood)[round(n_points*0.975)]))
  return(sampler)
}


# sigma comes from half cauchy
RWMH_importance_sampling_cauchy = function(x,y,beta, n_points) {
  #initialize
  current = beta
  S = diag(length(beta))
  coef = matrix(current, nrow=n_points, ncol=length(beta), byrow=TRUE)
  lupost = matrix(lupost_factory(x,y,beta,sigma=2), nrow=n_points, ncol=1, byrow=TRUE)
  likelihood = matrix(likelihood_factory(x,y,beta), nrow=n_points, ncol=1, byrow=TRUE)
  nn = 0
  for (i in 2:n_points) {
    sigma = abs(rcauchy(1, 0, 1))
    proposed = mvrnorm(1,coef[which.max(lupost),],2.4^2*S/length(beta))
    logr = lupost_factory(x,y,proposed,sigma)-lupost_factory(x,y,current,sigma)
    if (log(runif(1)) < logr) {current = proposed; nn = nn + 1}
    coef[i,] = current
    lupost[i] = lupost_factory(x,y,coef[i,],sigma)
    likelihood[i] = likelihood_factory(x,y,coef[i,])
    if (i%%50 == 0) S = var(coef[1:i,])
  }
  sampler = list(coefficient=coef, log_posterior=lupost, acceptance=nn/n_points,
                 coef_means = colMeans(coef), # MCMC estimate
                 coef_variance = colMeans(coef^2)-colMeans(coef)^2,
                 likelihood_est = mean(likelihood),
                 likelihood_CI = c(sort(likelihood)[round(n_points*0.025)],sort(likelihood)[round(n_points*0.975)]))
  return(sampler)
}

## ------------------ Model Index -----------------------------------------
## define the variable
model_index = matrix(c(0,0,0,0,0),nrow=2^4,ncol=5,byrow=TRUE)
model_index[,1] = 1
model_index[9:16,2] = 1
model_index[c(5:8,13:16),3] = 1
model_index[c(3:4,7:8,11:12,15:16),4] = 1
model_index[seq(2,2^4,2),5] = 1

## ---------- Example of full model (intercept, X1-X4 included) ----------------
i = 16
out <- glm(y~., data = logit[model_index[i,]==1], family = binomial, x = TRUE)
summary(out)
lupost <- lupost_factory(out$x, out$y, out$coefficients, sigma = 2)
print(lupost)

## Scenario 1: sigma = 2 for Problem 3.1

# coef_sampling
model_sampler = metropolis_hastings_normal(out$x, out$y, out$coefficients, n<-1000, scale<-0.2)

print(model_sampler$acceptance)
print(model_sampler$means)
print(model_sampler$standard_error)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$coefficient[,1], xlab="Iter", ylab=expression(beta[0]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,2], xlab="Iter", ylab=expression(beta[1]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,3], xlab="Iter", ylab=expression(beta[2]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,4], xlab="Iter", ylab=expression(beta[3]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,5], xlab="Iter", ylab=expression(beta[4]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$log_posterior, xlab="Iter", ylab='Log-posterior', main="Random Walk Metropolis-Hastings (RWMH)")

acf(model_sampler$coefficient[,1],100, plot=TRUE, main = expression(beta[0]))
acf(model_sampler$coefficient[,2],100, plot=TRUE, main = expression(beta[1]))
acf(model_sampler$coefficient[,3],100, plot=TRUE, main = expression(beta[2]))
acf(model_sampler$coefficient[,4],100, plot=TRUE, main = expression(beta[3]))
acf(model_sampler$coefficient[,5],100, plot=TRUE, main = expression(beta[4]))


## Scenario 2: sigma ~ Ca+(0,1) for Problem 3.2

# coef_sampling
model_sampler = metropolis_hastings_cauchy(out$x, out$y, out$coefficients, n<-1000, scale<-0.2)

print(model_sampler$acceptance)
print(model_sampler$means)
print(model_sampler$standard_error)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2,3,4,5,6),nrow=3,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$coefficient[,1], xlab="Iter", ylab=expression(beta[0]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,2], xlab="Iter", ylab=expression(beta[1]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,3], xlab="Iter", ylab=expression(beta[2]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,4], xlab="Iter", ylab=expression(beta[3]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$coefficient[,5], xlab="Iter", ylab=expression(beta[4]), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$log_posterior, ylim = c(-60,-47), xlab="Iter", ylab='Log-posterior', main="Random Walk Metropolis-Hastings (RWMH)")

acf(model_sampler$coefficient[,1],100, plot=TRUE, main = expression(beta[0]))
acf(model_sampler$coefficient[,2],100, plot=TRUE, main = expression(beta[1]))
acf(model_sampler$coefficient[,3],100, plot=TRUE, main = expression(beta[2]))
acf(model_sampler$coefficient[,4],100, plot=TRUE, main = expression(beta[3]))
acf(model_sampler$coefficient[,5],100, plot=TRUE, main = expression(beta[4]))


## ----------- Effect of scale of the proposal pdf on acceptance rate-----------
scale_vector = c(0.1,0.2,0.3,0.4,0.5,1,2,3)
for (i in 1:length(scale_vector)) {
  model_sampler = metropolis_hastings_normal(out$x, out$y, out$coefficients, n<-1000, scale<-scale_vector[i])
  print(data.frame(scale=scale_vector[i],
                   acceptance_rate=model_sampler$acceptance))
} 

## ------------------ MCMC of 16 Model/Variable Selection ----------------------
## RWMH sigma = 2
for (i in 1:length(model_index[,1])) {
  out <- glm(y~., data = logit[model_index[i,]==1], family = binomial, x = TRUE)
  model_sampler = metropolis_hastings_normal(out$x, out$y, out$coefficients, n<-1000, scale<-0.2)
  print(list(model_index=model_index[i,2:5],
             coef_means=model_sampler$coef_means,
             coef_variance=model_sampler$coef_variance,
             likelihood_est=model_sampler$likelihood_est,
             likelihood_CI=model_sampler$likelihood_CI))
}

## RWMH sigma comes from Ca+(0,1)
for (i in 1:length(model_index[,1])) {
  out <- glm(y~., data = logit[model_index[i,]==1], family = binomial, x = TRUE)
  model_sampler = metropolis_hastings_cauchy(out$x, out$y, out$coefficients, n<-1000, scale<-0.2)
  print(list(model_index=model_index[i,2:5],
             coef_means=model_sampler$coef_means,
             coef_variance=model_sampler$coef_variance,
             likelihood_est=model_sampler$likelihood_est,
             likelihood_CI=model_sampler$likelihood_CI))
}

## RWMH for importance sampling sigma = 2
for (i in 1:length(model_index[,1])) {
  out <- glm(y~., data = logit[model_index[i,]==1], family = binomial, x = TRUE)
  model_sampler = RWMH_importance_sampling_normal(out$x, out$y, out$coefficients, n<-1000)
  print(list(model_index=model_index[i,2:5],
             coef_means=model_sampler$coef_means,
             coef_variance=model_sampler$coef_variance,
             likelihood_est=model_sampler$likelihood_est,
             likelihood_CI=model_sampler$likelihood_CI))
}

## RWMH for importance sampling sigma comes from Ca+(0,1)
for (i in 1:length(model_index[,1])) {
  out <- glm(y~., data = logit[model_index[i,]==1], family = binomial, x = TRUE)
  model_sampler = RWMH_importance_sampling_cauchy(out$x, out$y, out$coefficients, n<-1000)
  print(list(model_index=model_index[i,2:5],
             coef_means=model_sampler$coef_means,
             coef_variance=model_sampler$coef_variance,
             likelihood_est=model_sampler$likelihood_est,
             likelihood_CI=model_sampler$likelihood_CI))
}

#sampler_1 = c(69.35,58.68,63.41,56.99,53.18,49.25,51.05,47.97,55.17,52.21,54.22,52.03,47.81,46.69,47.51,46.39)
#sampler_2 = c(69.38,59.50,63.72,57.92,54.51,50.16,52.02,49.50,55.65,52.85,55.08,52.06,49.12,47.57,48.85,46.66)
#sampler_3 = c(69.21,57.70,62.80,55.56,52.79,47.86,49.48,46.17,54.03,50.85,52.65,50.02,46.83,44.79,45.42,44.58)
#sampler_4 = c(69.24,57.70,63.39,55.56,52.14,48.00,51.38,46.96,55.74,52.44,53.84,50.02,48.63,44.79,45.42,43.83)
