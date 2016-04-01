#Basic example for Spectrode: Plot a mixture
#Change this to your working directory, where all code resides
setwd('C:/Git/Spectrode R/Code')

#parameters
num_clus <- 6
t <- seq(1,num_clus,length=num_clus)
w <- rep(1,num_clus)/num_clus
gamma <- 1/20

#data
n <- 1e4
p <- floor(gamma*n)
X <- matrix(rnorm(n*p),n) #random Gaussian
Sigma <-sample(t,p,replace=TRUE,prob=w) #covariance
X <- X%*%diag(Sigma^(1/2))
S <- t(X)%*%X/n
lambda <- eigen(S)$values
h<-hist(lambda, breaks = floor(5*sqrt(p)))

#Test Spectrode
library("deSolve")
source("spectrode.R")
s <- spectrode(t,w,gamma,ep = 1e-6)
m <- max(h$counts)/max(s$density)
matplot(s$grid,s$density*m, add=T,type="l",lwd=2, main="FPA")

#save fig
dev.print(pdf,"../Examples/Example 1 - Mixture Density/mixture_example.pdf")

