## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("dplyr")
library("tidyr")
library("ggplot2")
library("pscl") 

## ----set_seed, echo=FALSE------------------------------------------------
set.seed(2)

## --------------------------Problem1-------------------------------

## ----data----------------------------------------------------------------
n = 10
y = 3
a = b = 1
d = data.frame(x = seq(0,1,by=0.01)) %>% 
  mutate(prior = dbeta(x,a,b),
         "normalized likelihood" = dbeta(x, y, n-y),
         posterior = dbeta(x, a+y, b+n-y))

m <- d %>%
  tidyr::gather(Distribution, density, -x)

ggplot(m, aes(x, density, group = Distribution, linetype = Distribution, color= Distribution)) +
  geom_line() +
  theme_bw() + 
  theme(legend.position="bottom") 

## ----estimates, dependson='data'-----------------------------------------
estimates = data.frame(median = qbeta(.5, a+y, b+n-y),
                       mode = (a+y-1)/(a+b+n-2)) #posterior median (Median) and posterior mode (MAP)

ggplot(d, aes(x, posterior, group = 1)) +
  geom_line() +
  geom_vline(data = estimates %>% gather(estimator, value), 
             aes(xintercept=value, color=estimator, linetype=estimator),
             show.legend = TRUE) +
  theme_bw() +
  theme(legend.position="bottom")

## ----Problem 1.1-----------------------------------------
theta_p1 = seq(0.01,0.99,by=0.01)
theta_p1_est <- c()
for (i in seq(1,99,by=1)) {
  f_p1 = function(theta, theta_1 = theta_p1[i]) {
    (theta-theta_1)*theta^(a+y-1)*(1-theta)^(b+n-y-1)
  }
  theta_p1_est = append(theta_p1_est, integrate(f_p1, theta_p1[i], 1)$value - integrate(f_p1, 0, theta_p1[i])$value)
}
theta_p1[which.min(theta_p1_est)] # the best solution of Bayes estimate
plot(theta_p1, theta_p1_est, main="Posterior Expected Loss Function", xlab=expression(theta), ylab="Posterior expected loss")

## ----problem 1.2-----------------------------------------
library("pscl") 
delta <- seq(1,8,by=1) 
interval <- 0.2 * 2^(-delta)
theta_list <- c()
for (i in delta) {theta_list = append(theta_list, mean(betaHPD(a+y,b+n-y,interval[i])))}
plot(delta, abs(theta_list - estimates$mode), main="Posterior Estimates Error", xlab="i", ylab="Absolute Error")


## --------Problem 1.3: Non-conjugate Prior--------------------------------------------
n = 10
y = 3
f = function(theta) {
  theta^y*(1-theta)^(n-y)*exp(theta)
}
curve(f, col="red", lwd=2, 
      main="Binomial, nonconjugate prior", ylab="Density (proportional)", xlab=expression(theta))
legend("topright", c("Posterior"), col=c("red"), lwd=2)

## ----integrate, dependson=c('data','plot_f'), echo=TRUE------------------
(i = integrate(f, 0, 1))

## ----nonconjugate_plot_normalized, dependson='integrate'-----------------
curve(f(x)/i$value, col='red', lwd=2, ylab="Density", xlab=expression(theta))
legend("topright", "Posterior", col='red', lwd=2)

## ----nonconjugate_grid, fig.width=9--------------------------------------
w = 0.001
theta = seq(w/2, 1-w/2, by=w)
d = f(theta)
d = d/sum(d)/w # probability or pdf (with w divided)
plot(theta, d, type="l", col="red", lwd=2,  
     main="Binomial, nonconjugate prior", ylab="Density", xlab=expression(theta))
legend("topright", "Posterior", col="red", lwd=2)

## ----estimates-----------------------------------------------------------
estimates_non = data.frame(median = theta[which(cumsum(d)*w>0.5)[1]-1],
                           mode = theta[which.max(d)])

## ----problem 1.3.1-----------------------------------------
theta_p2 = seq(0.001,0.999,by=0.001)
theta_p2_est <- c()
for (i in seq(1,999,by=1)) {
  f_p2 = function(theta, theta_2 = theta_p2[i]) {
    (theta-theta_2)*theta^y*(1-theta)^(n-y)*exp(theta)
  }
  theta_p2_est = append(theta_p2_est, integrate(f_p2, theta_p2[i], 1)$value - integrate(f_p2, 0, theta_p2[i])$value)
}
theta_p2[which.min(theta_p2_est)] # the best solution of Bayes estimate
plot(theta_p2, theta_p2_est, main="Posterior Expected Loss Function", xlab="theta_Bayes", ylab="Posterior expected loss")


## ----problem 1.3.2-----------------------------------------
delta_non <- seq(1,8,by=1) 
interval_non <- 0.1 * 2^(-delta_non)
theta_non_list <- c()
theta_max = (cumsum(d)*w)[which.max(d)]
for (i in delta_non) {
  theta_non_list = append(theta_non_list, mean(theta[c(which(cumsum(d)*w>theta_max-interval_non[i])[1]-1, which(cumsum(d)*w>theta_max+interval_non[i])[1])]))
}
plot(delta_non, abs(theta_non_list - estimates_non$mode), main="Posterior Estimates Error", xlab="i", ylab="Absolute Error")