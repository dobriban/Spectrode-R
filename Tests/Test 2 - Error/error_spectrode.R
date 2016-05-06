#From Jeha Yang
#This code is supposed to give an error on a Mac, but not on a Windows
#setwd('C:/Git/Spectrode-R/Code')
setwd('~/Git/Spectrode-R/Code')

library(RMTstat)
library(deSolve)
source("spectrode.R")

gamma = 1  #gamma=1 does not work
supp.min = (1-sqrt(gamma))^2 ; supp.max = (1+sqrt(gamma))^2
n.bin = 1e2*5 # n.bin = 3 also resulted in a similar bug as below
t = seq(supp.min,supp.max,length.out=n.bin+1)
w = (pmp(t,svr=1/gamma) - pmp(c(0,t[1:n.bin]),svr=1/gamma))
# Approx. generalized MP F_(gamma,F_gamma) 
approx.gMP = spectrode(t,w,gamma)
approx.gMP$density
# > approx.gMP$density
#  [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
# [27] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
# [53] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
# [79] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0
#[105]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#[131]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#[157]  0  0  0  0  0  0  0  0  0  0  0  0  0

matplot(approx.gMP$grid,approx.gMP$density, type="l",lwd=2)

