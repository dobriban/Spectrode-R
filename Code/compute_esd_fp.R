compute_esd_fp <- function(t, w, gamma, epsilon = 1e-4, grid = seq(0.1,20,length=1e3)) {
  #Compute limit spectrum of covariance matrices with Fixed Point method
  
  #Inputs
  #t - population eigenvalues, non-negative vector of length p
  # w - mixture weights >=0 default: w <- uniform
  #gamma - aspect ratio p/n
  #epsilon - size of imaginary part of grid
  #grid - real grid where MP transform should be computed
  
  #outputs list containing
  # density - approximation to the density on the grid
  # m - numerical approximation to Stieltjes transform on real line, complex vector of size Nx1
  # v - numerical approximation to dual Stieltjes transform on real line, complex vector of size Nx1
  # numIter - number of iterations taken by the algorithm for each element of the grid
  # stepSize - last stepsize taken by the algorithm for each element of the grid
  # grid - useful if have a default

multiplier_num_iter <- 1e2
maxIter <- multiplier_num_iter/epsilon
tol <- epsilon
grid_imag <- grid + 1i*epsilon^2
L <- length(grid)
v <- matrix(0,L,1)
v_d <- matrix(0,L,6)
numIter <- matrix(0,L,1)
lastStepSize <- matrix(0,L,1)

for (i in 1:L) {
  #define function
  z <- grid_imag[i]
  fun <- function(x) - 1/(z-gamma*sum(w*t/(1+t*x)))
  
  #starting point
  if (i==1) {
    v1 <- -1/grid_imag[i]
  } else {
    v1 <- v[i-1]
  }
  v1 <- c(v1,fun(v1[1]))
  j <- 2
  
  while ((abs(v1[j]-v1[j-1])>tol)&&(j<maxIter)) {
    v1 <- c(v1,fun(v1[j]))
    j <- j+1
  }
  v[i] <- v1[j]
  numIter[i] <- j
  lastStepSize[i] <- abs(v1[j]-v1[j-1])
}

m <- 1/gamma*v- (1-1/gamma)/grid_imag
density= 1/pi*Im(m)
#[density,m,v,numIter,lastStepSize]
results <- list(
  "density" = density, "m" = m, "v" = v, "numIter" = numIter,
  "lastStepSize" = lastStepSize, "grid" = grid 
)
return(results)
}