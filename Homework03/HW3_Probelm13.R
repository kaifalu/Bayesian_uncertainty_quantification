## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("StanHeaders")
library("ggplot2")
library("rstan")
library("GGally")
library("dplyr")


## ------------------------- joint posterior distribution ----------------------- ##
# Marginalized distribution (alpha, beta)|y and then posterior density (theta)
model_location = "
data {
  int<lower=0> N;
  int<lower=0> n[N];
  int<lower=0> y[N];
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}

model {
  target += -5*log(alpha+beta)/2;
  y ~ beta_binomial(n,alpha,beta);
}

generated quantities {
  real<lower=0,upper=1> theta[N];
  for (i in 1:N) theta[i] = beta_rng(alpha+y[i], beta+n[i]-y[i]);
}
"

# data
d = data.frame(
  bicycles = c(16,9,10,13,19,20,18,17,35,55,212),
  vehicles = c(74,99,58,70,122,77,104,129,308,119,1160)
)
dat = list(y = d$bicycles, n = d$vehicles, N = nrow(d))

#sampling
m_location = stan_model(model_code=model_location)
r_location = sampling(m_location, dat, c("alpha","beta","theta"))


# sampling statistics of joint posterior distribution
r_location


# plots of marginal distribution: (alpha, beta)|y
samples = extract(r_location, c("alpha","beta"))
ggpairs(data.frame(log_alpha = log(as.numeric(samples$alpha)), log_beta = log(as.numeric(samples$beta))),
        lower = list(continuous="density")) + theme_bw()


# plot of posterior credible interval: theta|y
plot(r_location, pars=c('theta[1]','theta[2]','theta[3]','theta[4]','theta[5]','theta[6]','theta[7]','theta[8]','theta[9]','theta[10]','theta[11]'))


# plot of posterior density distribution: theta|y
posterior = extract(r_location, c('theta[1]','theta[2]','theta[3]','theta[4]','theta[5]','theta[6]','theta[7]','theta[8]','theta[9]','theta[10]')) %>%
  as.data.frame() %>%
  tidyr::gather(parameter, value) %>%
  mutate(parameter = as.factor(parameter))

ggplot(posterior,
       aes(x = value,
           fill = parameter)) +
  geom_density(alpha = 0.5) +
  labs(x=expression('theta'), y='Density', title='Posterior Distribution') +
  theme_bw()


# compare the posterior distribution to the raw proportions
tmp = data.frame(summary(r_location)$summary[,c(1,4,8)])
new_d = mutate(d[1:10,],
               lcl = tmp[3:12,2],
               ucl = tmp[3:12,3],
               ml = tmp[3:12,1])
raw_proportions = new_d$bicycles/new_d$vehicles
ggplot(new_d,
       aes(x     = raw_proportions,
           xend  = raw_proportions,
           y     = lcl,
           yend  = ucl)) +
  geom_point(aes(y = ml), colour = 'black', size = 3) +
  geom_segment(lwd=1) +
  geom_segment(x = 0,
               y = 0,
               xend = 0.5,
               yend = 0.5,
               size = 0.5,
               col = 'grey',
               linetype = "dashed") +
  scale_color_brewer(palette='Set1') +
  labs(x='Raw Proportion', y='Posterior Estimate', title='95% Credible Intervals') +
  xlim(c(0,0.5)) + ylim(c(0,0.6)) + 
  theme_bw()


## the combined mean of all locations
tmp_mean = data.frame(summary(r_location)$summary[13,])
tmp_mean