# Zero inflated poisson model

# Author: Federico Boiocchi

# Y_i ~ ZIP (mu_i,pi_i)

# Y_i response variable of a Poisson regression
# mu_i mean of the Poisson
# pi_i probability of getting a false zero (a recorded zero that is not
# truly zero)

# Y_i is modeled as a Zero-inflated poisson if its PDF
# is a mixture of an atomic mass probability function located in
# zero with probability P(y_i=0) = pi_i + (1-pi_i)*exp(-mu_i)
# and with probability 1-P(y_i=0) a poisson pdf with mean mu_i

# mu_i is modeled with a poisson glm: mu_i=exp(eta_i)
# pi_i is modeled with a logistic glm: logit(pi_i)=(gamma_i)

# eta_i and gamma_i are two linear predictiors

# zero inflated poisson simulation

# Y_i ~ ZIP (mu_i,pi_i)

# mu_i : mean of the poisson distribution of Y_i
# pi_i : probability of having a false zero

# Y_i :  0 P(y_i=0) = pi_i+(1-pi_i)*exp(-mu_i)
#        y ~ Pois(mu_i) 1-P(y_i=0)=(1-mu_i)*dpois(y_i,lambda=mu_i)

# eta_i = b0+b1*x1+b2*x2
# gamma_i = a0+a1*x3+a2*x4

# Estimation simulation using glm() and zeroinfl()

set.seed(123)

n <- 10000

a0 <- -1
b1 <- 0.5
b2 <- 1.3
b0 <- 2
a1 <- -0.3
a2 <- 1.1

tot <- 4*n
X <- matrix(data = rnorm(tot), nrow = n, ncol = 4)
x1 <- X[, 1]
x2 <- X[, 2]
x3 <- X[, 3]
x4 <- X[, 4]

eta_i <- b0 + b1 * x1 + b2 * x2
gamma_i <- a0 + a1 * x3 + a2 * x4

mu_i <- exp(eta_i)
pi_i <- plogis(gamma_i)

py0 <- pi_i + (1 - pi_i) * exp(-mu_i)
pneq0 <- 1 - py0

w <- rbinom(n, 1, prob = 1 - py0)

yi <- w*rpois(n, lambda = mu_i)

hist(yi, prob=TRUE,breaks=100)

mod0 <- glm(yi ~ x1+x2,family = poisson(link="log"))
summ0 <- summary(mod0)
logLik(mod0)

mod1 <- zeroinfl(yi ~ x1+x2 | x3+x4,
                 dist=c("poisson"),
                 link=c("logit"))
summ1 <- summary(mod1)
summ1$loglik


# Analysis forest fires

# install.packages("pscl")
library(pscl)
library(tidyverse)

forest <- read.csv("C:\\Users\\andre\\Desktop\\statistica computazionale\\HW_stat_comp\\forestfires.csv")
str(forest)
forest <- forest[-c(1:8, 12)] # select useful variables (Temp, Wind, RH)

View(forest)
# istogramma risposta
g1 <- ggplot(forest, aes(x = area)) +
  geom_histogram(aes(y = after_stat(ncount)), fill = "green3", color = "white", binwidth = 10) +
  labs(
    title = "Andamento Area",
    x = "Area bruciata (ettari)",
    y = "densità di freq."
  )
g1

# the response variable area is clearly zero inflated
# but it is a continuos random variable.

# discretization of the response

forest$area <- round(forest$area)

mod_zip <- zeroinfl(
  forest$area ~ forest$temp + forest$RH + forest$wind |
    forest$temp + forest$RH + forest$wind,
  dist = c("poisson"),
  link = c("logit")
)

summ_zip <- summary(mod_zip)
summ_zip$loglik
summ_zip$converged

# classic poisson regression
mod_glm <- glm(forest$area ~ forest$temp + forest$RH + forest$wind, family = poisson(link = "log"))
summ_glm <- summary(mod_glm)
k <- 4
loglik_glm <- (summ_glm$aic - 2 * k) / (-2)
logLik(mod_glm)

