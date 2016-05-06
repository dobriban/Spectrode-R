#Test the deSolve package
library("deSolve")
fun <-function(time,x,param) list(x)
v0 <- ode(y=0.1,times=seq(0,1,by=0.01),func=fun)
head(v0)
plot(v0)
