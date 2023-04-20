## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("coda")
library("ggplot2")
library("GGally")
library("dplyr")
library("rjags")


## ------------------------- joint posterior distribution ----------------------- ##
# Marginalized distribution (alpha, beta)|y and then posterior density (theta)
model_location = "
model {
## likelihood
for (i in 1:N) {
  y[i] ~ dnegbin(beta/(1+beta), alpha);
}

## hyperpriors

alpha <- u1/(u2^2);
beta <- (1-u1)/(u2^2);
u1 ~ dunif(0,1);
u2 ~ dunif(0,1);

# derived quantites
for (i in 1:N) {
  theta[i] ~ dgamma(alpha+y[i], beta+1)
}
}
"

# data
d = data.frame(
  vehicles = c(74,99,58,70,122,77,104,129,308,119)
)
dat = list(y = d$vehicles, N = nrow(d))

#sampling
m_location = jags.model(textConnection(model_location), dat)
r_location = coda.samples(m_location, c("alpha","beta","theta"), 1000)


# sampling statistics of joint posterior distribution
summary(r_location)
plot(r_location)
samples = as.matrix(r_location)


# plots of marginal distribution: (alpha, beta)|y
ggpairs(data.frame(log_alpha = log(as.numeric(samples[,"alpha"])), log_beta = log(as.numeric(samples[,"beta"]))),
        lower = list(continuous="density")) + theme_bw()


# scatter plots of alpha and beta
plot(samples[,"alpha"], samples[,"beta"], xlab = expression(alpha), ylab = expression(beta), pch=".", cex=2)
plot(log(samples[,"alpha"]/samples[,"beta"]), log(samples[,"alpha"]+samples[,"beta"]), xlab = expression(log(alpha~"/"~beta)), ylab = expression(log(alpha~"+"~beta)), pch=".", cex=2)


# plot of posterior density distribution: theta|y
posterior = samples[,c('theta[1]','theta[2]','theta[3]','theta[4]','theta[5]','theta[6]','theta[7]','theta[8]','theta[9]','theta[10]')] %>%
  as.data.frame() %>%
  tidyr::gather(parameter, value) %>%
  mutate(parameter = as.factor(parameter))

ggplot(posterior,
       aes(x = value,
           fill = parameter)) +
  geom_density(alpha = 0.5) +
  labs(x=expression('theta'), y='Density', title='Posterior Distribution') +
  theme_bw()