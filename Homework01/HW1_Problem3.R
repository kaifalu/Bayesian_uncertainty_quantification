## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
#library(plyr)
library("dplyr")
library("tidyr")
library("ggplot2")
library("pscl") 
#library(xtable)
#library(Sleuth3)

## ----set_seed, echo=FALSE------------------------------------------------
set.seed(2)

## ----Problem 3.1: mixture of conjugate priors--------------------------------
p = 0.5
a = c(5,2)
b = c(2,5)

prior = function(theta,p,a,b) {
  p*dbeta(theta,a[1],b[1]) + (1-p)*dbeta(theta,a[2],b[2])
}

curve(prior(x,p,a,b), col="blue", lwd=2,
      main="Mixture of conjugate Beta priors", ylab="Density", xlab=expression(theta))
curve(p*dbeta(x,a[1],b[1]), col="red", lty=2, add=TRUE)
curve((1-p)*dbeta(x,a[2],b[2]), col="purple", lty=2, add=TRUE)
legend("topright", c("Mixture","Prior 1", "Prior 2"), col=c("blue","red", "purple"), lwd=2)

## ----Problem 3.2: Posterior----------------------------------------
p = 0.5
a = c(5,2)
b = c(2,5)
n = 10
y = 3
ppd = function(y,n,a,b) {
  exp(lchoose(n,y)+lbeta(a+y,b+n-y)-lbeta(a,b)) # l-series: default natural logarithm
}
posterior = function(theta,p,a,b,y,n) {
  q = p*ppd(y,n,a[1],b[1])
  p = q/(q+(1-p)*ppd(y,n,a[2],b[2]))
  p*dbeta(theta,a[1]+y,b[1]+n-y) + (1-p)*dbeta(theta,a[2]+y,b[2]+n-y)
}

curve(posterior(x,p,a,b,y,n), col="red", lwd=2,
      main="Binomial, mixture of betas", ylab="Density", xlab=expression(theta)) # posterior

curve(prior(x,p,a,b), col="blue", lwd=2, add=TRUE) 

curve(p*dbeta(x,a[1],b[1]), col="blue", lty=2, add=TRUE)
curve((1-p)*dbeta(x,a[2],b[2]), col="blue", lty=2, add=TRUE)

q = p*ppd(y,n,a[1],b[1])
q = q/(q+(1-p)*ppd(y,n,a[2],b[2]))

curve(q*dbeta(x,a[1]+y,b[1]+n-y), col="red", lty=2, add=TRUE)
curve((1-q)*dbeta(x,a[2]+y,b[2]+n-y), col="red", lty=2, add=TRUE)

legend("topright", c("Prior","Posterior"), col=c("blue","red"), lwd=2)