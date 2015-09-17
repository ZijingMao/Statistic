## read data
y <- c(10, 3, 2, 5, 14)
## t.likelihood:
t.likelihood <- function(theta, y){
	m <- length(y)
	return(1/theta^m)
}
## t.prior:
t.prior <- function(x0, alpha, theta){
    return(alpha*x0^alpha/theta^(alpha+1))
}
## set up the MCMC, allocate space, etc
S <- 100100
x0 <- 5
alpha <- 3
theta <- rep(NA, S)
theta[1] <- 20
post <- t.likelihood(theta[1], y) * t.prior(x0, alpha, theta[1])
L <- 1
U <- 100
## MCMC rounds for theta
for(s in 2:S) {
	## propose new theta and calculate new posterior
	theta_prime <- runif(1, max(alpha, L*theta[s-1]/U), L*theta[s-1]*U)
	prime <- t.likelihood(theta_prime, y) * t.prior(x0, alpha, theta_prime)

	## accept or reject the proposal
	if(runif(1) < (prime/post) * (theta[s-1]/theta_prime)) {  ## accept
		theta[s] <- theta_prime
		post <- prime
	} else theta[s] <- theta[s-1]  ## reject
}
t <- theta[seq(101,S,100)]
par(mfrow=c(2,1))
plot(t, type="l", bty="n", xlab="trace", ylab="theta")
hist(t, xlab="theta", ylab="n", breaks=seq(3,20,0.1), cex=2)