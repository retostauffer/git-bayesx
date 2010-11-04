
setwd("c:/bayesx/r-test/linear")
library("VGAM")

# general set up for univariate responses
n <- 250
set.seed(777)

x1 <- round(runif(n, -pi, pi), 1)
x2 <- round(runif(n, -pi, pi), 1)

eta1 <- 0.4 + 0.5*x1 - 0.6*x2

mydata <- data.frame(x1, x2, eta)


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

# store data

mydata <- round(mydata, 3)
write.table(mydata, "data.raw", col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)


