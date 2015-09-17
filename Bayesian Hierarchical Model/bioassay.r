## load necessary libraries
library(arm)
library(coda)

## read in data
x = c(-0.86, -0.3, -0.05, 0.73)
n = c(5, 5, 5, 5)
y = c(0, 1, 3, 5)
data = cbind(x, n, y)

## ab.marg.llik:
ab.marg.llik <- function(alpha, beta, y, n, x){
	return( sum(dbinom (y, n, invlogit(alpha + beta*x), log=TRUE)) )
}

## ab.lprior:
ab.lprior <- function(alpha, beta){
    return(0)
}

## set up the MCMC, allocate space, etc
S <- 100000
alpha <- beta <- rep(NA, S)
alpha[1] <- 1
beta[1] <- 1 
lpost <- ab.marg.llik(alpha[1], beta[1], y, n, x) + ab.lprior(alpha[1], beta[1])
L <- 1
U <- 2

## MCMC rounds for alpha and beta
## Do as previous two problems
for(s in 2:S) {

  ## propose new alpha and calculate new log posteiror
  aprime <- runif(1, L*alpha[s-1]/U, L*alpha[s-1]*U)
  lprime <- ab.marg.llik(aprime, beta[s-1], y, n, x) + ab.lprior(aprime, beta[s-1])

  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * alpha[s-1]/aprime) {  ## accept
    alpha[s] <- aprime
    lpost <- lprime
  } else alpha[s] <- alpha[s-1]  ## reject
  
  ## propose new BETA and calculate new log posteiror
  bprime <- runif(1, L*beta[s-1]/U, U*beta[s-1]*L)
  lprime <- ab.marg.llik(alpha[s], bprime, y, n, x) + ab.lprior(alpha[s], bprime)

  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * beta[s-1]/bprime) {  ## accept
    beta[s] <- bprime
    lpost <- lprime
  } else beta[s] <- beta[s-1]  ## reject

}

## plot a useful transformation of the joint samples
plot(log(alpha/beta)[seq(100,S,100)], log(alpha + beta)[seq(100,S,100)], bty="n",
     xlab="log(alpha/beta)", ylab="log(alpha + beta)")

## ###################################################
## ab.marg.llik:
ab.marg.llik <- function(alpha, beta, sigma, y, n, x){
	theta <- rnorm( alpha + beta*x, sigma)
	return( sum(dbinom (y, n, invlogit(theta), log=TRUE)) )
}

## ab.lprior:
ab.lprior <- function(alpha, beta, sigma){
    return(0)
}

## set up the MCMC, allocate space, etc
S <- 100000
alpha <- beta <- sigma <- rep(NA, S)
alpha[1] <- 1
beta[1] <- 1 
sigma[1] <- 1
lpost <- ab.marg.llik(alpha[1], beta[1], sigma[1], y, n, x) + ab.lprior(alpha[1], beta[1])
L <- 1
U <- 2

## MCMC rounds for alpha and beta
for(s in 2:S) {

  ## propose new alpha and calculate new log posteiror
  aprime <- runif(1, L*alpha[s-1]/U, L*alpha[s-1]*U)
  lprime <- ab.marg.llik(aprime, beta[s-1], sigma[s-1], y, n, x) + ab.lprior(aprime, beta[s-1], sigma[s-1])

  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * alpha[s-1]/aprime) {  ## accept
    alpha[s] <- aprime
    lpost <- lprime
  } else alpha[s] <- alpha[s-1]  ## reject
  
  ## propose new beta and calculate new log posteiror
  bprime <- runif(1, L*beta[s-1]/U, U*beta[s-1]*L)
  lprime <- ab.marg.llik(alpha[s], bprime, sigma[s-1], y, n, x) + ab.lprior(alpha[s], bprime, sigma[s-1])

  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * beta[s-1]/bprime) {  ## accept
    beta[s] <- bprime
    lpost <- lprime
  } else beta[s] <- beta[s-1]  ## reject

  ## propose new sigma and calculate new log posteiror
  sprime <- runif(1, L*sigma[s-1]/U, U*sigma[s-1]*L)
  lprime <- ab.marg.llik(alpha[s], beta[s], sprime, y, n, x) + ab.lprior(alpha[s], beta[s], sprime)

  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * sigma[s-1]/sprime) {  ## accept
    sigma[s] <- sprime
    lpost <- lprime
  } else sigma[s] <- sigma[s-1]  ## reject
  
}

## plot a useful transformation of the joint samples
plot(log(alpha/beta)[seq(100,S,100)], log(alpha + beta)[seq(100,S,100)], bty="n",
     xlab="log(alpha/beta)", ylab="log(alpha + beta)"
## alpha
plot(alpha[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="alpha")
effectiveSize(alpha[-(1:10)])
## beta
plot(beta[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="beta")
effectiveSize(beta[-(1:10)])
## sigma
plot(sigma[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="sigma")
effectiveSize(sigma[-(1:10)])

## read in data and repeat
x = c(-0.86, -0.3, -0.05, 0.73)
n = c(5000, 5000, 5000, 5000)
y = c(500, 1000, 3000, 4500)
data = cbind(x, n, y)