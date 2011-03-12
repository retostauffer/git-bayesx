
# setwd("c:/bayesx/r-test/gmrf")
library("VGAM")
library("BayesX")

# general set up for univariate responses

m <- read.bnd("kreisesim.bnd")
centroids <- get.centroids(m)
n <- nrow(centroids)
nrep <- 15
centroids[, 2] <- (centroids[, 2] - min(centroids[, 2])) / (max(centroids[, 2]) - min(centroids[, 2]))
centroids[, 3] <- (centroids[, 3] - min(centroids[, 3])) / (max(centroids[, 3]) - min(centroids[, 3]))
centroids <- centroids[rep(1:n, each=nrep),]
n <- nrep*n

set.seed(777)

z1 <- runif(n, -1, 1)

eta1 <- -2 + 6*centroids[,2]*centroids[,3]
eta2 <- -0.7 + 2*centroids[,2]*centroids[,3] + z1*1.5*sin(centroids[,2]*2*pi)

mydata <- data.frame(region=as.numeric(rownames(centroids)), eta1, eta2, z1)

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

