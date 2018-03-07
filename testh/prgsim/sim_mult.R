
setwd("c:/bayesx/trunk/testh/testdata")
set.seed(123)

n <- 1000

x <- runif(n, -pi, pi)
z <- runif(n, -1, 1)

etatilde <- z^2 - mean(z^2)
eta <- sin(x)*exp(etatilde)

y <- rnorm(n, mean=eta, sd=0.1)

data <- data.frame(x=x, z=z, y=y, y2=y)
data <- round(data, 3)
write.table(data, "gaussian_exp_mult.raw", row.names=FALSE, col.names=TRUE, sep=" ", quote=FALSE)



