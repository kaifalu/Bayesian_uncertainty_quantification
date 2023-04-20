## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("xtable")

## ----set_seed------------------------------------------------------------
set.seed(1)

## ----Problem 4.1----------------
# define f(x)

f = function(x,log = FALSE) {
  out = (1+x^2)^(-1.75)*as.numeric(x>0)
  if (log) return(log(out))
  return(out)
}

I_f = integrate(function(x) f(x), -Inf, Inf)

# define g2(x)
g2 = function(x, log = FALSE) {
    out = x^(-3)*as.numeric(x>1)
    if (log) return(log(out))
    return(out)
}

I_g2 = integrate(function(x) g2(x), -Inf, Inf)

# define g(x)
g = function(x, log = FALSE) {
  out = 0.5*as.numeric(x>0 & x<1) + 0.5*g2(x)/I_g2$value
  if (log) return(log(out))
  return(out)
}

sampler_1 = function(n_points, w=0.5) {
  sample = matrix(0, nrow=n_points, ncol=1, byrow=TRUE)
  for (i in 1:n_points) {
    if (runif(1)<w) {sample[i] = runif(1)}
    else {sample[i] = sqrt(1/(1-runif(1)))}
  }
  return(sample)
}

# rejection sampling
set.seed(1)

M = 2/I_f$value
n = 10^2
d = data.frame(x = sampler_1(n), 
               u = runif(n)) %>%
  mutate(u_scaled = u*M*g(x),
         accept   = u_scaled < f(x)/I_f$value)

gg = ggplot(d, aes(x=x,y=u_scaled,col=accept)) + 
  geom_point()

clrs = unique(ggplot_build(gg)$data[[1]]$colour)

gg + stat_function(fun=function(x) M*g(x), col=clrs[2]) +
  stat_function(fun=function(x) f(x)/I_f$value, col=clrs[1]) + 
  labs(x="sample",y=expression(paste("u M g(",x,")"))) +
  theme_bw()

cat("Observed acceptance rate was",mean(d$accept))



## ----Problem 4.2----------------

# define g2(x)
g2 = function(x,p,c, log = FALSE) {
  out = x^(-p)*(p-1)*as.numeric(x>c)/(c^(1-p))
  if (log) return(log(out))
  return(out)
}

g = function(x,p,c,w, log = FALSE) {
  out = w*(1/c)*as.numeric(x>0 & x<c) + (1-w)*g2(x,p,c)
  if (log) return(log(out))
  return(out)
}

g_cdf = function(x,p,c,w, log = FALSE) {
  out = w*(x/c)*as.numeric(x>0 & x<c) + (w+(1-w)*(1-(x/c)^(1-p)))*as.numeric(x>c)
  if (log) return(log(out))
  return(out)
}

sampler_2 = function(n_points, p, c, w) {
  sample = matrix(0, nrow=n_points, ncol=1, byrow=TRUE)
  for (i in 1:n_points) {
    if (runif(1)<w) {sample[i] = runif(1)*c}
    else {sample[i] = c*(1-runif(1))^(1/(1-p))}
  }
  return(sample)
}

# sampling
set.seed(2)
w = 0.5
c = 1
p = 3
n = 10^2
M = c/(w*I_f$value)

# plot cdf and ecdf
x = sampler_2(n,p,c,w)
x_cdf = g_cdf(x,p,c,w)
plot(ecdf(x),col="blue", main = 'cdf(x) and ecdf(x)')
lines(sort(x),sort(x_cdf), add = TRUE, col="black", lwd=5)
legend(8, 0.8, legend = c('ecdf', 'cdf'),fill = c('blue','black'))

# plot absolute error between cdf and ecdf
plot(x,abs(ecdf(x)(x)-x_cdf), ylab = 'Absolute error', main = 'The difference between cdf(x) and ecdf(x)')

# sampling from g
d = data.frame(x = x, 
               u = runif(n)) %>%
  mutate(u_scaled = u*M*g(x,p,c,w),
         accept   = u_scaled < f(x)/I_f$value)

gg = ggplot(d, aes(x=x,y=u_scaled,col=accept)) + 
  geom_point()

clrs = unique(ggplot_build(gg)$data[[1]]$colour)

gg + stat_function(fun=function(x) M*g(x,p,c,w), col=clrs[1]) +
  stat_function(fun=function(x) f(x)/I_f$value, col=clrs[2]) + 
  labs(x="sample",y=expression(paste("u M g(",x,")"))) +
  theme_bw()

cat("Observed acceptance rate was",mean(d$accept))




## ------Problem 4.3------------------------------

M = function(p,c,w) {
  m1 = c/(w*I_f$value)
  m2 = c^(1-p)*(p/(3.5-p))^(p/2)*(1+p/(3.5-p))^(-1.75)/(I_f$value*(1-w)*(p-1))
  m3 = c*(1+c^2)^(-1.75)/(I_f$value*(1-w)*(p-1))
  return(max(m1,m2,m3))
}

p_vector = seq(1.1,3.4,0.1)
c_vector = 1:10
w_vector = seq(0.1,0.9,0.1)
mm = array(rep(1, length(p_vector)*length(c_vector)*length(w_vector)), 
           dim=c(length(p_vector), length(c_vector), length(w_vector)))
for (i in 1:length(p_vector)) {
  for (j in 1:length(c_vector)) {
    for (k in 1:length(w_vector)) {
      mm[i,j,k] = M(p_vector[i],c_vector[j],w_vector[k])
    }
  }
}

## mm[13,1,8],M(2.3,1,0.8)
F_m1 <- matrix(c(1, 2, 3), nrow = 1, ncol = 3, byrow = TRUE)
layout(F_m1)
plot(p_vector,1/mm[,1,8],xlab='p',ylab='Efficiency',main='Rejection sampler') # effect of p
plot(c_vector,1/mm[13,,8],xlab='c',ylab='Efficiency',main='Rejection sampler') # effect of c
plot(w_vector,1/mm[12,1,],xlab='w',ylab='Efficiency',main='Rejection sampler') # effect of c




## ------Problem 4.4------------------------------

h = function(x, log = FALSE) {
  out = x*as.numeric(x>=10)
  if (log) return(log(out))
  return(out)
}

E_hf_x = integrate(function(x) f(x)*h(x), -Inf, Inf)

# sampling with the greatest efficiency
set.seed(3)
J = 10^5

w = 0.8
c = 1
p = 2.3
MM = M(p,c,w)

x = sampler_2(J,p,c,w)

F_m2 <- matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE)
layout(F_m2)

# sampling from proposal g
d = data.frame(x = x, 
               u = runif(J)) %>%
  mutate(u_scaled = u*MM*g(x,p,c,w),
         accept   = u_scaled < f(x)/I_f$value)

gg = ggplot(d, aes(x=x,y=u_scaled,col=accept)) + 
  geom_point()

clrs = unique(ggplot_build(gg)$data[[1]]$colour)

gg + stat_function(fun=function(x) MM*g(x,p,c,w), col=clrs[2]) +
  stat_function(fun=function(x) f(x)/I_f$value, col=clrs[1]) + 
  labs(x="sample",y=expression(paste("u M g(",x,")"))) +
  theme_bw()

cat("Observed acceptance rate was",mean(d$accept))

# calculate h(x)
hx = h(x[d$accept])
cummean = function(x) return(cumsum(x)/(1:length(x)))
plot(cmhx<-cummean(hx), type='l', ylab="Estimate", main="Monte Carlo estimate", xlab="Number of samples")

# variance or 95%CI
cumvar = function(x) {
  J = length(x)
  cumvar = numeric(J)
  for (j in 1:J) { cumvar[j] = var(x[1:j])/j }
  return(cumvar)
}

cvhx = cumvar(hx)
uci  = cmhx+qnorm(.975)*sqrt(cvhx)
lci  = cmhx-qnorm(.975)*sqrt(cvhx)
plot(cmhx, main="Monte Carlo estimate", type="n",
     ylim=range(uci,lci,na.rm=T), xlab="Number of samples", ylab="Estimate")
abline(h=E_hf_x$value, col="red",lwd=2)
lines(cmhx,lwd=2)
lines(uci, lty=2)
lines(lci, lty=2)
legend("bottomright", c("Truth","Estimate","95% CI"),
       lwd=c(2,2,1), lty=c(1,1,2), col=c("red","black","black"))


# plot relative error

set.seed(3)
J = 10^6

w = 0.8
c = 1
p = 2.3
MM = M(p,c,w)

x = sampler_2(J,p,c,w)

# sampling from proposal g
d = data.frame(x = x, 
               u = runif(J)) %>%
  mutate(u_scaled = u*MM*g(x,p,c,w),
         accept   = u_scaled < f(x)/I_f$value)

cat("Observed acceptance rate was",mean(d$accept))

# calculate h(x)
hx = h(x[d$accept])
cmhx = cummean(hx)
cvhx = cumvar(hx)
uci  = cmhx+qnorm(.975)*sqrt(cvhx)
lci  = cmhx-qnorm(.975)*sqrt(cvhx)

relative_error = (uci-lci)/(2*cmhx)
plot(relative_error, type='l', ylim = c(0,2), ylab="Relative error", main="Relative error of Monte Carlo estimate", xlab="Number of samples")
lines(c(rep(0.1, length(relative_error))),lty=2, col = 'blue')
lines(c(rep(0.01, length(relative_error))),lty=2, col = 'blue')
text(10,0.2,'0.1',col='blue')
text(10,0,'0.01',col = 'blue')

if (relative_error[length(relative_error)] < 0.1) {
  cat("Required MC sample size for relative error less than 10% is", 
      round(which.max(relative_error<0.1)/mean(d$accept)))
} else {
  cat("Required MC sample size for relative error less than 10% is", 
      "higher than 10^6")
}

if (relative_error[length(relative_error)] < 0.01) {
  cat("Required MC sample size for relative error less than 1% is", 
      round(which.max(relative_error<0.01)/mean(d$accept)))
} else {
  cat("Required MC sample size for relative error less than 1% is", 
      "higher than 10^6")
}



## -------------- Problem 4.5 ------------------

gj = function(x,j,c=10, log = FALSE) {
  out = x^(-j)*(j-1)*as.numeric(x>c)/(c^(1-j))
  if (log) return(log(out))
  return(out)
}

sampler_j = function(n_points, j, c=10) {
  sample = matrix(0, nrow=n_points, ncol=1, byrow=TRUE)
  for (i in 1:n_points) {
    sample[i] = c*(1-runif(1))^(1/(1-j))
  }
  return(sample)
}

# sample size
J = 10^5
E_hf_x = integrate(function(x) f(x)*h(x), -Inf, Inf) #truth

# sampling from proposal g (j=2)
set.seed(100000)
j=2
x = sampler_j(J,j)
weights = f(x)/(I_f$value*gj(x,j))

# calculate h(x)
hx = h(x)
cummean = function(x) return(cumsum(x*weights)/(1:length(x)))
cmhx<-cummean(hx)

# variance or 95%CI
cumvar = function(x) {
  J = length(x)
  cumvar = numeric(J)
  for (i in 1:J) { cumvar[i] = var(x[1:i]*weights[1:i])/i }
  return(cumvar)
}

cvhx = cumvar(hx)
uci  = cmhx+qnorm(.975)*sqrt(cvhx)
lci  = cmhx-qnorm(.975)*sqrt(cvhx)
plot(cmhx[100:20000], main="Monte Carlo estimate", type="n",
     ylim=range(uci[100:length(uci)],lci[100:length(lci)],E_hf_x$value,na.rm=T), xlab="Number of samples", ylab="Estimate")
abline(h=E_hf_x$value, col="red",lwd=2)
lines(cmhx[100:20000],lwd=2)
lines(uci[100:20000], lty=2)
lines(lci[100:20000], lty=2)
legend("bottomright", c("Truth","Estimate","95% CI"),
       lwd=c(2,2,1), lty=c(1,1,2), col=c("red","black","black"))


# sampling from proposal g (j=3)
set.seed(200000)
j=3
x = sampler_j(J,j)
weights = f(x)/(I_f$value*gj(x,j))

# calculate h(x)
hx = h(x)
cmhx<-cummean(hx)

# variance or 95%CI
cvhx = cumvar(hx)
uci  = cmhx+qnorm(.975)*sqrt(cvhx)
lci  = cmhx-qnorm(.975)*sqrt(cvhx)
plot(cmhx[100:20000], main="Monte Carlo estimate", type="n",
     ylim=range(uci[100:length(uci)],lci[100:length(lci)],E_hf_x$value, na.rm=T), xlab="Number of samples", ylab="Estimate")
abline(h=E_hf_x$value, col="red",lwd=2)
lines(cmhx[100:20000],lwd=2)
lines(uci[100:20000], lty=2)
lines(lci[100:20000], lty=2)
legend("bottomright", c("Truth","Estimate","95% CI"),
       lwd=c(2,2,1), lty=c(1,1,2), col=c("red","black","black"))


# sampling from a proposal g (j=4)
set.seed(300000)
j=4
x = sampler_j(J,j)
weights = f(x)/(I_f$value*gj(x,j))

# calculate h(x)
hx = h(x)
cmhx<-cummean(hx)

# variance or 95%CI
cvhx = cumvar(hx)
uci  = cmhx+qnorm(.975)*sqrt(cvhx)
lci  = cmhx-qnorm(.975)*sqrt(cvhx)
plot(cmhx[100:20000], main="Monte Carlo estimate", type="n",
     ylim=range(uci[100:length(uci)],lci[100:length(lci)],E_hf_x$value,na.rm=T), xlab="Number of samples", ylab="Estimate")
abline(h=E_hf_x$value, col="red",lwd=2)
lines(cmhx[100:20000],lwd=2)
lines(uci[100:20000], lty=2)
lines(lci[100:20000], lty=2)
legend("bottomright", c("Truth","Estimate","95% CI"),
       lwd=c(2,2,1), lty=c(1,1,2), col=c("red","black","black"))
