
# setwd("c:/bayesx/r-test/pspline")
library("VGAM")

# general set up for univariate responses
n <- 1000
set.seed(777)

x1 <- round(runif(n, -pi, pi), 1)
x2 <- round(runif(n, -1, 1), 1)
x3 <- round(runif(n, -1, 1), 1)

z1 <- rbinom(n, prob=0.5, size=1)

eta1 <- sin(x1) + x2^2
eta2 <- sin(x1) + x2^2 + z1*x3

mydata <- data.frame(x1, x2, x3, z1, eta1, eta2)

# Gaussian responses

mydata$norm1 <- rnorm(n, eta1, sd=0.5)
mydata$norm2 <- rnorm(n, eta2, sd=0.5)

# Bernoulli responses

mydata$bernlogit1 <- rbinom(n, prob=plogis(eta1), size=1)
mydata$bernlogit2 <- rbinom(n, prob=plogis(eta2), size=1)

mydata$bernprobit1 <- rbinom(n, prob=pnorm(eta1), size=1)
mydata$bernprobit2 <- rbinom(n, prob=pnorm(eta2), size=1)

mydata$berncloglog1 <- rbinom(n, prob=pgumbel(eta1), size=1)
mydata$berncloglog2 <- rbinom(n, prob=pgumbel(eta2), size=1)

# Binomial responses

mydata$bin1 <- rbinom(n, prob=plogis(eta1), size=5)
mydata$bin2 <- rbinom(n, prob=plogis(eta2), size=5)

mydata$binw <- rep(5, n)

# Poisson responses

mydata$pois1 <- rpois(n, lambda=exp(eta1))
mydata$pois2 <- rpois(n, lambda=exp(eta2))

# store data

mydata <- round(mydata, 3)
write.table(mydata, "data.raw", col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)

