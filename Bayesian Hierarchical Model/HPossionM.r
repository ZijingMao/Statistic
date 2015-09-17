## load necessary libraries
library(coda)

## read in the data file
bike <- read.table("G:/Program Files (x86)/RStudio-0.98.1062/demo/bike.txt", header=TRUE)
View(bike)
attach(bike, warn.conflicts=FALSE)
y <- bikes+other
m <- length(y)

## ab.marg.llik:
ab.marg.llik <- function(alpha, beta, y){
	m <- length(y)
	const <- m*(alpha * log(beta) - lgamma(alpha))
	lprod <- sum(lgamma(alpha + y) - (alpha + y) * log(beta + 1))
	return(const + lprod)
}

## ab.lprior:
ab.lprior <- function(alpha, beta){
    return((alpha + beta)^(-1))
}

## set up the MCMC, allocate space, etc
S <- 100000
alpha <- beta <- rep(NA, S)
alpha[1] <- 1
beta[1] <- 1 
lpost <- ab.marg.llik(alpha[1], beta[1], y) + ab.lprior(alpha[1], beta[1])
L <- 1
U <- 2

## MCMC rounds for alpha and beta
for(s in 2:S) {

  ## propose new alpha and calculate new log posteiror
  aprime <- runif(1, L*alpha[s-1]/U, L*alpha[s-1]*U)
  lprime <- ab.marg.llik(aprime, beta[s-1], y) + ab.lprior(aprime, beta[s-1])

  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * alpha[s-1]/aprime) {  ## accept
    alpha[s] <- aprime
    lpost <- lprime
  } else alpha[s] <- alpha[s-1]  ## reject
  
  ## propose new BETA and calculate new log posteiror
  bprime <- runif(1, L*beta[s-1]/U, U*beta[s-1]*L)
  lprime <- ab.marg.llik(alpha[s], bprime, y) + ab.lprior(alpha[s], bprime)

  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * beta[s-1]/bprime) {  ## accept
    beta[s] <- bprime
    lpost <- lprime
  } else beta[s] <- beta[s-1]  ## reject

}

##
## trace plots of the sampled parameters, and ESS
##

## alpha
plot(alpha[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="alpha")
effectiveSize(alpha[-(1:10)])
## beta
plot(beta[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="beta")
effectiveSize(beta[-(1:10)])
## plot the joint samples
plot(alpha[seq(100,S,100)], beta[seq(100,S,100)],
     bty="n", xlab="alpha", ylab="beta")
## plot a useful transformation of the joint samples
plot(log(alpha/beta)[seq(100,S,100)], log(alpha + beta)[seq(100,S,100)], bty="n",
     xlab="log(alpha/beta)", ylab="log(alpha + beta)")
	 
## allocate space for the samples
theta <- matrix(NA, nrow=S, ncol=m)
## gather the samples by direct MC
for(j in 1:m) {
  theta[,j] <- rgamma(S, alpha + y[j], beta + 1)
}
## calculate the posterior means of the theta samples
theta.mean <- apply(theta[-(1:100),], 2, mean)

## theta
plot(theta[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="theta")
effectiveSize(theta[-(1:10)])

## plot the posterior means to assess shrinkage
plot(theta.mean, y, bty="n", ylab="theta-bar", xlab="y")
points(theta.mean[m], y[m], col="red", pch=19)
abline(0,1)

## plotting the difference between theta-bar
## and y/n as a function of n
plot(n, theta.mean - y, bty="n", xlab="n", ylab="theta-bar - y")
points(n[m], theta.mean[m] - y[m], col="red", pch=19)
abline(h=0)

## plotting the 95% CI of theta
ylim <- range(c(0, 320))
plot(y, theta.mean, bty="n", ylim=ylim, ylab="95% PI for theta", xlab="y")
abline(0,1)
for(j in 1:m) {
  temp <- theta[seq(1,S,100), j]
  segments((y)[j], (sort(temp)[c(25,975)])[1], 
			(y)[j], (sort(temp)[c(25,975)])[2])
}

## was beta reasonable for theta
theta.pooled <- apply(theta[-(1:100),], 1, mean)
dm <- density(theta.pooled)
dab <- density((alpha/beta)[-(1:100)])
ylim <- range(c(dm$y, dab$y))
plot(dm, lwd=2, bty="n", main="",
     xlab="theta")
lines(dab, col=2, lty=2, lwd=2)
legend("topright", c("bike", "a/b"), col=1:2,
       lwd=2, lty=1:2, bty="n")