
setwd("c:/bayesx/r-test/linear")
library("VGAM")
library("foreign")
library("MASS")
library("gamlss.dist")

# general set up for univariate responses
n <- 2500
set.seed(777)

x1 <- round(runif(n, -pi, pi), 1)
x2 <- round(runif(n, -pi, pi), 1)

eta1 <- 0.4 + 0.5*x1 - 0.6*x2
eta2 <- -0.4 - 0.5*x1 + 0.6*x2


mydata <- data.frame(x1, x2, eta1)


# Gaussian responses

mydata$norm1 <- rnorm(n, eta1, sd=0.5)


# Bernoulli responses

mydata$bernlogit1 <- rbinom(n, prob=plogis(eta1), size=1)

mydata$bernprobit1 <- rbinom(n, prob=pnorm(eta1), size=1)

mydata$berncloglog1 <- rbinom(n, prob=pgumbel(eta1), size=1)


# Binomial responses

mydata$bin1 <- rbinom(n, prob=plogis(eta1), size=5)

mydata$binw <- rep(5, n)


# Poisson responses

mydata$pois1 <- rpois(n, lambda=exp(eta1))


# Cumulative probit with 3 categories

uhelp <- rnorm(n,eta1,sd=1) 
mydata$cump <- 1*(uhelp<=-0.6) + 2*(uhelp>-0.6 & uhelp <= 0.6) + 3*(uhelp > 0.6)


# Multinomial logit

pi1 <- exp(eta1)/(1+exp(eta1)+exp(eta2))
pi2 <- exp(eta2)/(1+exp(eta1)+exp(eta2))
pi3 <- 1-pi1-pi2
uhelp <- runif(n)


mydata$mlog <- 1*(uhelp <= pi1) + 2*(uhelp<=pi1+pi2 & uhelp > pi1) + 3*(uhelp > pi1+pi2)


# Negative Binomial

mydata$nbin <- rnegbin(n,mu=exp(eta1),theta =rep(0.6,n))


# ZIP

mydata$zip <- rZIP(n,mu=exp(eta1),sigma =rep(0.2,n))


# store data

mydata <- round(mydata, 3)
write.table(mydata, "data.raw", col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)
write.dta(mydata,"data.dta")



