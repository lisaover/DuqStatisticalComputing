M <- 5000
xc <- runif(1,0,1)

xValues <- NULL

for (i in 1:M) {

xstar <- runif(1,xc-0.2,xc+0.2)
while (xstar <= 0 | xstar >= 1) {
xstar <- runif(1,xc-0.2,xc+0.2)
}

xc <- xstar

xValues <- c(xValues, xc)

}

hist(xValues)

count2 <- 0
xc <- runif(1,0,1)

xValues2 <- NULL

for (j in 1:M) {

xstar <- runif(1,xc-0.2,xc+0.2)
while (xstar <= 0 | xstar >= 1) {
xstar <- runif(1,xc-0.2,xc+0.2)
}

alpha <- max(1/0.4, 1/(xstar+0.2), 1/(1-xstar+0.2))/max(1/0.4, 1/(xc+0.2), 1/(1-xc+0.2))

if (alpha >= 1) {
xc <- xstar
count2 <- count2 + 1
} else {
x_p <- runif(1,0,1)
if (x_p < alpha) {
xc <- xstar
count2 <- count2 + 1
}
}

xValues2 <- c(xValues2, xc)
}

hist(xValues2)


---------------------------------------------------------------

data <- read.table(file.choose(), header=F)
colnames(data) <- NULL
dataAvg <- sapply(data,mean)
dataVar <- sapply(data,var)

mu <- dataAvg
sigmaSq <- dataVar

N <- 10000

muValues <- NULL
sigmaSqValues <- NULL

muCount <- 0
sigmaSqCount <- 0

b <- 0.25
c <- 0.25

for (i in 1:N) {

muStar <- rnorm(1,mu,b^2)
sigmaSqStar <- rnorm(1,sigmaSq,c^2)

muExpComp <- 155*(dataAvg - muStar)^2 + muStar^2/10 - (155*(dataAvg - mu)^2 + mu^2/10)
muAlphaComp <- exp(-muExpComp/(2*sigmaSq))

if (muAlphaComp >= 1) {
mu <- muStar
muCount <- muCount + 1
} else {
mu_p <- runif(1,0,1)
if (mu_p < muAlphaComp) {
mu <- muStar
muCount <- muCount+1
}
}
muValues <- c(muValues, mu)

AStar <- (1/sigmaSqStar)^(155/2+4)
BStar <- (1/(2*sigmaSqStar))*(154*dataVar+155*(dataAvg-muStar)^2+muStar^2/10+2)
A <- (1/sigmaSq)^(155/2+4)
B <- (1/(2*sigmaSq))*(154*dataVar+155*(dataAvg-muStar)^2+muStar^2/10+2)
sigmaSqExpComp <- log(AStar)-BStar-log(A)+B
hastingsNum <- dnorm(sigmaSq,sigmaSqStar,c)/(1-pnorm(0,sigmaSqStar,c))
hastingsDenom <- dnorm(sigmaSqStar,sigmaSq,c)/(1-pnorm(0,sigmaSq,c))
hastings <- hastingsNum/hastingsDenom
sigmaSqAlphaComp <- exp(sigmaSqExpComp)*hastings

if (sigmaSqAlphaComp >= 1) {
sigmaSq <- sigmaSqStar
sigmaSqCount <- sigmaSqCount+1
} else {
sigmaSq_p <- runif(1,0,1)
if (sigmaSq_p < sigmaSqAlphaComp) {
sigmaSq <- sigmaSqStar
sigmaSqCount <- sigmaSqCount+1
}
}
sigmaSqValues <- c(sigmaSqValues, sigmaSq)
}

par(mfrow=c(2,2))
ts.plot(muValues,xlab="Iteration")
ts.plot(sigmaSqValues,xlab="Iteration")
hist(muValues,probability=T)
hist(sigmaSqValues,probability=T)
xx <- seq(0,2.5,length=250)
lines(xx, 42.45304^80/gamma(80)*(1/xx)^81*exp(-42.45304/xx))
lines(seq(0,60,length=250),1/dgamma(seq(0,60,length=250),43/2-1,21*dataVar))