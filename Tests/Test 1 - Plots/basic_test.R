#Basic tests for Spectrode and its dependencies.
#Change this to your working directory, where all code resides
setwd('C:/Git/Spectrode R/Code')

#parameters
t <- c(1,2,10)
w <- c(1,1,1)/3
gamma <- 1/2

#data
n <- 1e3
p <- floor(gamma*n)
X <- matrix(rnorm(n*p),n) #random Gaussian
Sigma <-sample(t,p,replace=TRUE,prob=w) #covariance
X <- X%*%diag(Sigma^(1/2))
S <- t(X)%*%X/n
lambda <- eigen(S)$values
h<-hist(lambda, breaks = floor(5*sqrt(p)))

#test Fixed-Point method
source("compute_esd_fp.R")
s <- compute_esd_fp(t,w,gamma)
m <- max(h$counts)/max(s$density)
matplot(s$grid,s$density*m, add=T,type="l",lwd=2, main="FPA")
#save fig
dev.print(pdf,"../Tests/FPA_example.pdf")

#Test Inverse ST
source("evaluate_inverse_ST.R")
grid <- seq(-1.1, 0.1, length=1e3)
ist <- evaluate_inverse_ST(t,w,gamma,grid)
plot(grid, ist$grid)
#save fig
dev.print(pdf,"../Tests/inv_ST_example.pdf")


#Test Spectrode
library("deSolve")
source("spectrode.R")
h <- hist(lambda, breaks = floor(5*sqrt(p)))
s <- spectrode(t,w,gamma,ep = 1e-4)
m <- max(h$counts)/max(s$density)
matplot(s$grid,s$density*m, add=T,type="l",lwd=2, main="FPA")
#save fig
dev.print(pdf,"../Tests/spectrode_example.pdf")

