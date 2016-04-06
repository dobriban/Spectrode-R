#Basic example for Spectrode: Plot a mixture
#Change this to your working directory, where all code resides
#setwd('C:/Git/Spectrode-R/Code')
setwd('~/Git/Spectrode-R/Code')
library("deSolve")
source("spectrode.R")

num_clus <- 2
t <- c(1,6)
eps <- 0.05
w <- c(1-eps,eps)
gamma <- 1
s <- spectrode(t,w,gamma)

matplot(s$grid,s$density, type="l",lwd=2)
dev.print(pdf,"../Examples/Example 2 - STAT 325 Lec 4/twopoint_mixture_example.pdf")

## vary epsi
eps_g <- c(0.1,0.25)
for (i in (1:length(eps_g))) {
w <- c(1-eps_g[i],eps_g[i])
s <- spectrode(t,w,gamma)
matplot(s$grid,s$density, add = T, col = i+1, type="l",lwd=2)
}
legend("topright", legend = c(eps,eps_g), col=1:3, pch=1)
dev.print(pdf,"../Examples/Example 2 - STAT 325 Lec 4/twopoint_mixture_example_2.pdf")

