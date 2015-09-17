y <-  0
a1 <- a2 <- 50
b1 <- b2 <- 10
K <- 1000
THETA <- array(NA,dim=c(5,K,2))
for(n in 1:5) {
	theta2 <- runif(1, min=40, max=60)
	for(i in 1:K) {
		mu1 <- (a1+b1^2*(y-theta2))/(1+b1^2)
		sigmasq1 <-  b1^2/(1+b1^2)
		theta1 <- rnorm(1, mean=mu1, sd=sqrt(sigmasq1))
		mu2 <- (a2+b2^2*(y-theta1))/(1+b2^2)
		sigmasq2 <-  b2^2/(1+b2^2)
		theta2 <- rnorm(1, mean=mu2, sd=sqrt(sigmasq2))
		THETA[n,i,] <- c(theta1, theta2)
	}
}

pdf(file="set2014_6problem2c7.pdf",height=4,width=12,points=12)
par(mfrow=c(1,3))
acf(THETA[1,,1],main=expression(theta[1]))
acf(THETA[1,,2],main=expression(theta[2]))
acf(THETA[1,,1]+THETA[1,,2],main=expression(mu))
dev.off()

pdf(file="set2014_6problem2c8.pdf",height=4,width=12,points=12)
par(mfrow=c(1,3))
plot(THETA[1,,1], type="l", lty=1, xlab="Iteration", ylab=expression(theta[1]), ylim=c(-100,0))
for(i in 2:5) lines(THETA[i,,1], lty=i)
plot(THETA[1,,2], type="l", lty=1, xlab="Iteration", ylab=expression(theta[2]), ylim=c(0,100))
for(i in 2:5) lines(THETA[i,,2], lty=i)
plot(THETA[1,,1]+THETA[1,,2], type="l", lty=1, xlab="Iteration", ylab=expression(mu), ylim=c(-5,5) )
for(i in 2:5) lines(THETA[i,,1]+THETA[i,,2], lty=i)
dev.off()

mean(THETA[1,,1])
## -65.29 (K=100, b=100)
## -32.66 (K=1000, b=100)
## -17.62 (K=100, b=10)
##    1.77 (K=1000, b=10)
sd(THETA[1,,1])
##   8.00 (K=100, b=100)
## 12.05 (K=1000, b=100)
## 13.22 (K=100, b=10)
##   9.80 (K=1000, b=10)
mean(THETA[1,,2])
## 65.35 (K=100, b=100)
## 32.67 (K=1000, b=100)
## 13.21 (K=100, b=10)
##  -1.28 (K=1000, b=10)
sd(THETA[1,,2])
##  8.11 (K=100, b =100)
##12.02 (K=1000,b=100)
##13.21 (K=100, b=10)
##  9.73 (K=1000, b=10)
mean(THETA[1,,1]+THETA[1,,2])
## 0.06 (K=100, b=100)
## 0.02 (K=1000, b=100)
## 0.14 (K=100, b=10)
## 0.49 (K=1000, b=10)
sd(THETA[1,,1]+THETA[1,,2])
## 1.02 (K=100, 100)
## 1.02 (K=1000, b=100)
## 0.87 (K=100, b=10)
## 1.02 (K=1000, b=10)

## Code for class

y <-  0
a1 <- a2 <- 50
b1 <- b2 <- 10
K <- 1000
theta2 <- 1000
THETA <- array(NA,dim=c(1,K,2))
for(i in 1:K) {
	mu1 <- (a1+b1^2*(y-theta2))/(1+b1^2)
	sigmasq1 <-  b1^2/(1+b1^2)
	theta1 <- rnorm(1, mean=mu1, sd=sqrt(sigmasq1))
	mu2 <- (a2+b2^2*(y-theta1))/(1+b2^2)
	sigmasq2 <-  b2^2/(1+b2^2)
	theta2 <- rnorm(1, mean=mu2, sd=sqrt(sigmasq2))
	THETA[1,i,] <- c(theta1, theta2)
}

par(mfrow=c(2,2))
acf(THETA[1,,1]+THETA[1,,2],lag.max=500,main=expression(mu))
plot(THETA[1,,1], type="l", lty=1, xlab="Iteration", ylab=expression(theta[1]), ylim=c(-100,100))
plot(THETA[1,,2], type="l", lty=1, xlab="Iteration", ylab=expression(theta[2]), ylim=c(-100,100))
plot(THETA[1,,1]+THETA[1,,2], type="l", lty=1, xlab="Iteration", ylab=expression(mu), ylim=c(-5,5) )
cat(" Mean of theta1: ", mean(THETA[1,,1]), "\n",
	"SD of theta1: ", sd(THETA[1,,1]), "\n",
	"Mean of theta2: ", mean(THETA[1,,2]), "\n",
	"SD of theta2: ", sd(THETA[1,,2]), "\n",
	"Mean of theta1+theta2: ", mean(THETA[1,,1]+THETA[1,,2]), "\n",
	"SD of theta1+theta2: ", sd(THETA[1,,1]+THETA[1,,2]), "\n")
