## read in the data file
library(coda)
salary <- read.table("G:/Program Files (x86)/RStudio-0.98.1062/demo/salary.txt", header=TRUE)
View(salary)
attach(salary, warn.conflicts=FALSE)
## Plotting Outcome (Salary) as a function of Covariates (Years)
plot(Salary,Years,pch=19)

## MLE fit of Regression Model
model <- lm(Salary~Years)
summary(model)

par(mfrow=c(1,1))
plot(model$fitted.values,model$residuals,pch=19)
abline(h=0,col=2)

## Sampling from Posterior Distribution of coefficients beta and variance sigsq
beta.hat <- model$coef
n <- length(Salary)
p <- length(beta.hat)
s2 <- (n-p)*summary(model)$sigma^2
V.beta <- summary(model)$cov.unscaled

numsamp <- 1000
beta.samp <- matrix(NA,nrow=numsamp,ncol=p)
sigsq.samp <- rep(NA,numsamp)
for (i in 1:numsamp){
	temp <- rgamma(1,shape=(n-p)/2,rate=s2/2)
	cursigsq <- 1/temp
	curvarbeta <- cursigsq*V.beta
	curvarbeta.chol <- t(chol(curvarbeta))
	z <- rnorm(p,0,1)
	curbeta <- beta.hat+curvarbeta.chol%*%z
	sigsq.samp[i] <- cursigsq
	beta.samp[i,] <- curbeta
}

## Plotting Posterior Samples for coefficients beta and variance sigsq
par(mfrow=c(2,2))
hist(beta.samp[,1],main="Intercept",col="gray")
abline(v=beta.hat[1],col="red")
hist(beta.samp[,2],main="Years",col="gray")
abline(v=beta.hat[2],col="red")
hist(sigsq.samp,main="Sigsq",col="gray")
abline(v=summary(model)$sigma^2,col="red")

## Posterior Predictive Sampling for new data with year = 4
Xstar <- c(1,4)  # new tree with year = 4

Ystar.samp <- rep(NA,numsamp)
for (i in 1:numsamp){
 Ystar.hat <- sum(beta.samp[i,]*Xstar) 
 Ystar.samp[i] <- rnorm(1,mean=Ystar.hat,sd=sqrt(sigsq.samp[i]))
}   

## Plotting Predictive Samples and Data on same scale
par(mfrow=c(2,1))
xmin <- min(Salary,Ystar.samp)
xmax <- max(Salary,Ystar.samp)
hist(Salary,main="Salary in Observed Data",xlim=c(xmin,xmax),breaks=10,col="gray")
abline(v=mean(Salary),col=2)
hist(Ystar.samp,main="Predicted Salary with Years = 4",xlim=c(xmin,xmax),col="gray")
abline(v=mean(Ystar.samp),col=2)

## Lasso Part
library(lars)
model.lasso <- lars(as.matrix(Years), Salary, type="lasso")
lambda.lasso <- c(model.lasso$lambda,0)
beta <- coef(model.lasso)