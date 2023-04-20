## ----libraries, message=FALSE, warning=FALSE-----------------------------
library("coda")

## ----set_seed------------------------------------------------------------
set.seed(2)

## ----random_walk---------------------------------------------------------
rar1 = function(n, mu, rho, sigma, theta0) {
  theta = rep(theta0, n+1)
  for (i in 1:n) theta[i+1] = rnorm(1, mu+rho*(theta[i]-mu), sigma)
  theta
}

n_vector = c(100,1000,10000)
rho_vector = c(0,0.25,0.5,0.75,0.9,0.95,0.99)

## method (i)
var_eff_sample_size = function(n, rho) {
  round(n^2*(1-rho)^2/(n*(1-rho^2)-2*rho*(1-rho^n)))
}

for (i in n_vector){
  d = ddply(data.frame(rho=c(0,0.25,0.5,0.75,0.9,0.95,0.99),n=c(i,i,i,i,i,i,i)), .(rho,n), function(effective_size) data.frame(effective_size=var_eff_sample_size(i,effective_size$rho)))
  print(ddply(d, .(rho,n), effective_size=round(d$effective_size)))
}



## method (ii)
mu = 0; sigma = 1
for (i in n_vector){
  d = ddply(data.frame(rho=c(0,0.25,0.5,0.75,0.9,0.95,0.99),n=c(i,i,i,i,i,i,i)), .(rho,n), function(x) data.frame(x=rar1(x$n,mu,x$rho,sigma,0)))
  dd = ddply(d, .(rho,n), summarize, acf_sum = sum(acf(x,i,plot=FALSE)$acf))
  print(data.frame(rho=c(0,0.25,0.5,0.75,0.9,0.95,0.99),n=c(i,i,i,i,i,i,i), effective_size=round(dd$n[1]/(1+2*dd$acf_sum))))
}


## method (iii)
mu = 0; sigma = 1
for (i in n_vector){
  d = ddply(data.frame(rho=c(0,0.25,0.5,0.75,0.9,0.95,0.99), n=c(i,i,i,i,i,i,i)), .(rho,n), function(x) data.frame(x=rar1(x$n,mu,x$rho,sigma,0)))
  print(ddply(d, .(rho,n), summarize, 
        effective_size = round(coda::effectiveSize(x))))
}
