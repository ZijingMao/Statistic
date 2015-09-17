## read in the data file
library(coda)
library(MCMCpack)
coal_ <- read.table("G:/Program Files (x86)/RStudio-0.98.1062/demo/coal_data.txt", header=TRUE)
attach(coal_, warn.conflicts=FALSE)
View(coal_)
## Plotting original data
par(mfrow=c(1,1))
plot(Year,Num,pch=19)

## k.post.dist
k.post.dist <- function(t1, t2, X) {
	n<-length(X)
	post.dist<-rep(0,n)
	for (i in 1:n)
		post.dist[i]<-exp(-(t1-t2)*i)*(t1/t2)^sum(X[1:i])
	post.dist<-post.dist/sum(post.dist)
	return(post.dist)
}

## Sampling from Posterior Distribution of coefficients
Y <- Num
n <- length(Y)
a1 <- a2 <- 5
c1 <- c2 <- 1
d2 <- 2
d1 <- d2 - 1

numsamp <- 5000
lambda1.samp <- rep(0,numsamp)
lambda2.samp <- rep(0,numsamp)
beta1.samp <- rep(0,numsamp)
beta2.samp <- rep(0,numsamp)
k.samp <- rep(0,numsamp)

beta1 <- rgamma(1, c1, d1)
beta2 <- rgamma(1, c2, d2)
k <- sample(1:n, 1)
for (s in 1:numsamp){
	lambda1.samp[s] <- lambda1 <- rgamma(1, a1 + sum(Y[1:k]), beta1 + k)
	lambda2.samp[s] <- lambda2 <- rgamma(1, a2 + sum(Y[(1+k):n]), beta2 + n - k)
	beta1.samp[s] <- beta1 <- rgamma(1, a1 + c1, lambda1 + d1)
	beta2.samp[s] <- beta2 <- rgamma(1, a2 + c2, lambda2 + d2)
	prob.k <- k.post.dist(lambda1, lambda2, Y)
	k.samp[s] <- k <- sample(1:n, 1, prob = prob.k)
}

hist(k.samp, breaks=50, xlim=c(30, 50), main="Decided k", col="gray")
abline(v=mean(k.samp),col=2)

