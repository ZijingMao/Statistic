## read in the data file
library(coda)
library(MCMCpack)
library(ggplot2)
BeerData <- read.table("G:/Program Files (x86)/RStudio-0.98.1062/demo/BeerData.txt", header=TRUE)
attach(BeerData, warn.conflicts=FALSE)
View(BeerData)

## Plotting Raw Data of Scores and Drinks (Types)
df <- data.frame(Obs, Budlight, BlueMoon, Guiness, DosEquis, Shiner)
ggplot(df, aes(Obs, y = Scores, color = variable)) + 
     geom_point(aes(y = Budlight, col = "Budlight")) + 
     geom_point(aes(y = BlueMoon, col = "BlueMoon")) +
     geom_point(aes(y = Guiness, col = "Guiness")) +
     geom_point(aes(y = DosEquis, col = "DosEquis")) +
     geom_point(aes(y = Shiner, col = "Shiner"))

#######################################
## Process in a hierarchical model
#######################################
##
## Gibbs Sampling from Condition Distribution
##
numsamp <- 1000
J <- 5 	# How many Types
N <- 30	# How many participates

y <- matrix(NA, nrow=N, ncol=J)
theta.samp <- matrix(NA,nrow=numsamp,ncol=J)
sigma <- rep(NA, J)
mu.samp <- rep(NA,numsamp)
tau_2.samp <- rep(NA,numsamp)

## store data into y
y <- cbind(Budlight, BlueMoon, Guiness, DosEquis, Shiner)
pred <- y[37,]
y <- y[1:36,]
## initialize value
sigma <- rep(sqrt(var(c(Budlight, BlueMoon, Guiness, DosEquis, Shiner))/J), J)
y.bar <- apply(y, 2, mean)
mu <- sum(y.bar)/J
tau_2 <- rep(0.1, J)
for (s in 1:numsamp){
	V <- 1/(1/(sigma^2) + 1/(tau_2))
	theta.hat <- (y.bar/(sigma^2) + mu/(tau_2)) * V
	theta.samp[s,] <- theta <- rnorm(5, theta.hat, sqrt(V))
	
	theta.bar <- mean(theta.samp[s,])
	mu.samp[s] <- mu <- rnorm(1, theta.bar, sqrt(tau_2/J))
	
	tau_2.samp[s] <- tau_2 <- 
		rinvgamma(1, J/2-1, sum((mu-theta.samp[s,])^2/2))
}

##
## trace plots of the sampled parameters
##
par(mfrow=c(1,2))
## mu
plot(mu.samp[seq(1,numsamp,10)], type="l", bty="n", xlab="s", ylab="mu")
## tau
tau <- sqrt(tau_2.samp[seq(1,numsamp,10)])
plot(tau, type="l", bty="n", xlab="s", ylab="tau")
## plot the joint samples
par(mfrow=c(1,1))
plot(log(mu.samp/tau), log(mu.samp + tau), bty="n",
     xlab="log(mu/tau)", ylab="log(mu + tau)")

##	 
## plotting the 95% CI of theta
##
ylim <- range(c(0, 5))
plot(apply(theta.samp, 2, mean), bty="n", ylim=ylim, 
	ylab="95% PI for theta", xlab="Drink Types", pch=19, col="blue", xaxt = "n")
axis(1, at=1:5, labels=c("Budlight", "Blue Moon", "Guiness", "Dos Equis", "Shiner"))
abline(mean(apply(theta.samp, 2, mean)), 0, col="red")
for(j in 1:J) {
  temp <- theta.samp[, j]
  segments(j, (sort(temp)[c(25,975)])[1], 
			j, (sort(temp)[c(25,975)])[2])
}
legend(x=4,y=5,c("Mean = 2.744"),lty = 1, col="red")

##	 
## compare a new data to our 95% CI of theta
##
par(new=TRUE)
ylim <- range(c(0, 5))
plot(pred, pch=17, col="green", bty="n", ylim=ylim, 
	ylab="", xlab="", main=NULL, xaxt = "n", yaxt="n")

	
#######################################
## Process in a 
## multivariate Gaussian model
#######################################
##
## Gibbs Sampling from Condition Distribution
##	 
n <- nrow(y)
k <- ncol(y)
bar.x <- colMeans(y)
S <- rowSums(apply(y, 1, function(x) (x - bar.x) %*% t(x - bar.x)))
S <- matrix(S, nrow=k, ncol=k)

Sigma <- array(0, dim=c(k, k, numsamp))
Mu <- matrix(0, ncol=k, nrow=numsamp)
for(i in 1:numsamp){
	Sigma[,,i] <- riwish(n, S)
	Mu[i,] <- mvrnorm(1, bar.x, Sigma[,,i]/n)
}

##
## trace plots of the sampled parameters
##
## Mu
df.Mu <- data.frame(1:100,Mu[seq(1,numsamp,10),1],Mu[seq(1,numsamp,10),2],
	Mu[seq(1,numsamp,10),3],Mu[seq(1,numsamp,10),4],Mu[seq(1,numsamp,10),5])
ggplot(df.Mu, aes(X1.100, y = Trace, color = variable)) + 
  geom_line(aes(y=Mu[seq(1,numsamp,10),1], col = "Budlight")) +  # first layer
  geom_line(aes(y=Mu[seq(1,numsamp,10),2], col = "BlueMoon")) +  # second layer
  geom_line(aes(y=Mu[seq(1,numsamp,10),3], col = "Guiness")) +
  geom_line(aes(y=Mu[seq(1,numsamp,10),4], col = "DosEquis")) +
  geom_line(aes(y=Mu[seq(1,numsamp,10),5], col = "Shiner")) +
  ylab("Mu") + xlab("trace")
	
##	 
## plotting the 95% CI of mu for each individual model
##
CI.Mu <- matrix(NA, nrow=2, ncol=J)
for(j in 1:J) {
	CI.Mu[,j] <- quantile(Mu[,j], c(0.025, 0.975))
}
ylim <- range(c(0, 5))
Mu.hat = apply(Mu, 2, mean)
plot(Mu.hat, bty="n", ylim=ylim, 
	ylab="95% PI for theta", xlab="Drink Types", pch=19, col="blue", xaxt = "n")
axis(1, at=1:5, labels=c("Budlight", "Blue Moon", "Guiness", "Dos Equis", "Shiner"))
abline(mean(Mu.hat), 0, col="red")
for(j in 1:J) {
	segments(j, CI.Mu[1,j], 
			j, CI.Mu[2, j])
}
legend(x=4,y=5,c("Mean = 2.742"),lty = 1, col="red")

##
## compare the same new data to our 95% CI of mu 
## for each individual model
##
points(pred, pch=17, cex=1, col=rgb(0,0.8,0,1))
	

#######################################
## Process models individually
## As a comparison
#######################################
##	 
## for each drink we build a single model, 
## assuming the joint prior proportion to 1/sigma^2
##
CI <- matrix(NA, nrow=2, ncol=J)
mu.old.hat <- rep(NA, J)
for(j in 1:J) {
	S = sum((y[, j] - y.bar[j])^2)
	sigma2.old = S/rchisq(1000, N - 1)
	mu.old = rnorm(1000, y.bar[j], sqrt(sigma2.old)/sqrt(N))
	CI[,j] <- quantile(mu.old, c(0.025, 0.975))
	mu.old.hat[j] <- mean(mu.old)
}

##	 
## plotting the 95% CI of mu for each individual model
##
ylim <- range(c(0, 5))
plot(mu.old.hat, bty="n", ylim=ylim, 
	ylab="95% PI for theta", xlab="Drink Types", pch=19, col="blue", xaxt = "n")
axis(1, at=1:5, labels=c("Budlight", "Blue Moon", "Guiness", "Dos Equis", "Shiner"))
abline(mean(mu.old.hat), 0, col="red")
for(j in 1:J) {
  segments(j, CI[1,j], 
			j, CI[2, j])
}
legend(x=4,y=5,c("Mean = 2.738"),lty = 1, col="red")

##
## compare the same new data to our 95% CI of mu 
## for each individual model
##
par(new=TRUE)
ylim <- range(c(0, 5))
plot(pred, pch=17, col="green", bty="n", ylim=ylim, 
	ylab="", xlab="", main=NULL, xaxt = "n", yaxt="n")

#######################################
## Evaluation of our three built models
#######################################
##	 
## combine three results 
## to show a better comparison representation
##	
xlim <- range(c(1, 5.4))
ylim <- range(c(0, 5))
## plot the hierarchical model
plot(1:5, apply(theta.samp, 2, mean), bty="n", xlim=xlim, ylim=ylim, 
	ylab="95% PI for theta", xlab="Drink Types", pch=19, cex=1.5, col=rgb(0.5,0.2,0.1,1), xaxt = "n")
axis(1, at=1:5, labels=c("Budlight", "Blue Moon", "Guiness", "Dos Equis", "Shiner"))
abline(mean(apply(theta.samp, 2, mean)), 0, col="red", lwd=1)
for(j in 1:J) {
  temp <- theta.samp[, j]
  segments(j, (sort(temp)[c(25,975)])[1], 
			j, (sort(temp)[c(25,975)])[2], lwd=2)
}
## plot the each individual models
points(1.2:5.2, mu.old.hat, pch=19, cex=1.5, col=rgb(0.2,0.1,0.5,1))
for(j in 1:J) {
  segments(j+0.2, CI[1,j], 
			j+0.2, CI[2, j], lwd=2)
}
abline(mean(mu.old.hat), 0, col="blue", lwd=1)
## plot the each individual models
points(1.1:5.1, Mu.hat, pch=19, cex=1.5, col=rgb(0.1,0.5,0.1,1))
for(j in 1:J) {
  segments(j+0.1, CI.Mu[1,j], 
			j+0.1, CI.Mu[2, j], lwd=2)
}
abline(mean(Mu.hat), 0, col="green", lwd=1)
points(1.1:5.1, pred, pch=17, cex=1.5, col=rgb(0,0.8,0,1))
## legend
legend(x=2.5,y=5.2,c("MVN Mean = 2.742", "HM Mean = 2.744", "IM Mean = 2.738"),lty = c(1, 1, 1), col=c("green", "red", "blue"), lwd = c(2, 2, 2))
