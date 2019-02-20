
setwd("c:/bayesx/trunk/testh")

set.seed(123)

n <- 100000

x1 <- runif(n, -pi, pi)
x2 <- runif(n, -pi, pi)
x3 <- runif(n, -pi, pi)
x4 <- runif(n, -pi, pi)
x5 <- runif(n, -pi, pi)

eta <- sin(x1) + cos(x2) + x3 - x4
p <- 1/(1+exp(-eta))
y <- rbinom(n, size=1, p)

data <- data.frame(x1, x2, x3, x4, x5, eta, p, y)

write.table(data, "testdata/large_binomial.raw", col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE)
