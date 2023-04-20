## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("reshape2")
library("plyr")
library("ggplot2")

## ----set_seed------------------------------------------------------------
set.seed(0)

## ----normal approximation to normal posterior--------------------------------------
n=10
mut = 17
sigt = 1
y = rnorm(n, mut, sigt)

## ----Problem 3.1--------------------------------------
mu_seq <- seq(16, 18, by=0.1)
log_sigma_seq <- seq(-1, 1, by=0.1)
log_posterior <- matrix(, nrow = length(mu_seq), ncol = length(log_sigma_seq))
for (i in 1:length(mu_seq)) {
  for (j in 1:length(log_sigma_seq)) {
    log_posterior[i,j] = log(prod(dnorm(y, mean=mu_seq[i], sd=exp(log_sigma_seq[j]))))
  }
}

persp(mu_seq, log_sigma_seq, log_posterior, xlab=expression(mu), ylab=expression(tau), zlab="p(mu,tau|y)",
      main='log-posterior', col='pink', shade=.4, theta = 30, phi = 15, ticktype='detailed')

column_max=as.integer(which.max(log_posterior)/length(log_sigma_seq))
row_max=which.max(log_posterior) - column_max*length(mu_seq)
mut_hat=mu_seq[row_max]
taut_hat=log_sigma_seq[column_max+1]

## ----Problem 3.2&3.3--------------------------------------

# Approximation: "clear box"
mu_hat = mean(y)
s_square = sum((y-mean(y))^2)/(n-1)
tau_hat = log(sqrt((n-1)/n)*sqrt(s_square))
mu_var = s_square/n
tau_var = 1/(2*n)

#curve(dnorm(x,mu_hat,mu_var), xlim=c(16.5,18.5), lwd=2, xlab=expression(mu), ylab=expression(paste("p(", mu,"|y)")))
#curve(dnorm(x,tau_hat,tau_var), lwd=2, xlab=expression(tau), ylab=expression(paste("p(", tau,"|y)")))

mean_mu_tau_hat = c(mu_hat, tau_hat)
mat.data <- c(mu_var, 0, 0, tau_var)
clear_mu_tau_var = matrix(mat.data, nrow = 2)

# Approximation: "black box"
obj <- function(x) {-1*(-n*x[2] - 1/(2*exp(x[2])^2)*((n-1)*s_square + n*(mu_hat-x[1])^2))} # form of fun f
obj_origin <- function(x) {-n*x[2] - 1/(2*exp(x[2])^2)*((n-1)*s_square + n*(mu_hat-x[1])^2)}
init = c(16, 0.5)
out1 = optim(par=init, obj, lower=c(16,-1), upper=c(18,1), method="L-BFGS-B")

eps = 10^(-3)
e_vector = matrix(c(1,0,0,1), nrow = 2, ncol =2)
Hess_fxx <- function(g, x, matrix, delta=eps) {
  Hess_mat <- c()
  for (i in 1:2) {
    for (j in 1:2) {
      Hess_mat = append(Hess_mat, (g(x+delta*matrix[i,])+g(x-delta*matrix[j,])-2*g(x))/delta^2)
    }
  }
  Hess_mat
}

myJ = Hess_fxx(obj_origin,out1$par,e_vector)
myJ = matrix(myJ, nrow=2, ncol=2)
black_mu_tau_var = solve(myJ)

#curve(dnorm(x,out1$par[1],black_mu_tau_var[1]), xlim=c(15,20), lwd=2, xlab=expression(mu), ylab=expression(paste("p(", mu,"|y)")))
#curve(dnorm(x,out1$par[2],black_mu_tau_var[4]), xlim=c(-0.5,1), lwd=2, xlab=expression(tau), ylab=expression(paste("p(", tau,"|y)")))


curve(dnorm(x,mu_hat,mu_var), xlim=c(15,20), lwd=2, col="blue", xlab=expression(mu), ylab=expression(paste("p(", mu,"|y)")))
curve(dnorm(x,out1$par[1],black_mu_tau_var[1]), lwd=2, col="red", add = TRUE)
legend("topright", c("Clear-box Gaussian approximation","Black-box Gaussian approximation"), col=c("blue","red"), lwd=2)

curve(dnorm(x,tau_hat,tau_var), xlim=c(-0.5,1), lwd=2, col="blue", xlab=expression(tau), ylab=expression(paste("p(", tau,"|y)")))
curve(dnorm(x,out1$par[2],black_mu_tau_var[4]), lwd=2, col="red", add = TRUE)
legend("topright", c("Clear-box Gaussian approximation","Black-box Gaussian approximation"), col=c("blue","red"), lwd=2)
