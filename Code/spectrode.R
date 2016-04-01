#' Spectrode: Compute the limit empirical spectral distribution of large sample covariance matrices
spectrode <- function(t = c(1,2,3), w = c(1,1,1)/3, gamma = 1/2,ep = 1e-4) {
  # Inputs
  # t - population eigenvalues >=0
  # gamma - aspect ratio p/n
  # w - mixture weights >=0; default: w = uniform
  # ep - accuracy parameter: square root of the imaginary part of the grid where MP is solved
  
  # Outputs: List containing: 
  # grid - real grid where Stieltjes transform is evaluated, real vector of size Nx1
  # density - approximation to the density on the grid
  # m - numerical approximation to Stieltjes transform on real line,
  #       complex vector of size Nx1
  # v - numerical approximation to dual Stieltjes transform on real line,
  #       complex vector of size Nx1
  # K_hat - numerical approximation to number of disjoint clusters of support
  # l_hat - lower endpoints of support intervals; real vector of size K_hat;
  # u_hat - upper endpoints of support intervals; real vector of size K_hat;
  
  #source two necessary files
  source("compute_esd_fp.R")
  source("evaluate_inverse_ST.R")
  
  #Error checking
  if (any(t < 0)) {
    stop("Positive eigenvalues required")
  }
  if (any(w < 0)) {
    stop("Non-negative weights required")
  } 
  if (length(w)!=length(t)) {
    stop("Weight vector must have same weight as vector of eigenvalues")
  }
  if (gamma < 0) {
    stop("Positive aspect ratio required")
  }
  if (ep < 0) {
    stop("Positive accuracy required")
  }

# define auxiliary variables
M <- floor(sqrt(1/ep))+3
T <- unique(t[t>0])
B <- sort(-1/T)
num_unique <- length(B)
z_l_endpoints <- matrix(0,num_unique+2,1)
z_u_endpoints <- matrix(0,num_unique+2,1)

#Part I: outside the support
z_grid <- matrix(0,num_unique+2,M) #grids in z-space
v_grid <- matrix(0,num_unique+2,M) #grids in v-space
num_grid_points <- matrix(0,num_unique+2,1) #number of points in each grid segment

#go interval by interval;
num_clus <- 1
z_l_endpoints[1] <- 0
if (gamma <1) { #if need to look at lower half
  v_L <- min(B) -1 
  v <- seq(v_L+ep, min(B)-ep, length=M)
  ist <- evaluate_inverse_ST(t,w,gamma,v)
  #if f is decreasing on the current interval, need to widen it
  g_grid <- ist$grid[ist$grid<Inf] 
  while (ist$maxf==g_grid[1])  { #need to deal with Inf's
    v_L <- 10*(v_L - min(B)) + min(B)
    v<- seq(v_L+ep, min(B)-ep, length=M)
    ist <- evaluate_inverse_ST(t,w,gamma,v)
    g_grid <- ist$grid[ist$grid<Inf]
  }
  z_u_endpoints[num_clus] <- ist$maxf
  num_grid_points[num_clus] <- ist$ind_max #number of grid points in first interval
  z_grid[num_clus,1:ist$ind_max] <- ist$grid[1:ist$ind_max]
  v_grid[num_clus,1:ist$ind_max] <- v[1:ist$ind_max]
} else {
  num_grid_points[num_clus] <- 0
}

#each interval between point masses -1/t
for (i in (1:(num_unique-1))) {
  v_L <- B[i] #lower endpoint -1/t
  v_U <- B[i+1] #upper endpoint -1/t
  if (length(v_L)==0) {
    v_L <-2*v_U
  }
  #Find the edges of the spectrum
  v <- seq(v_L+ep, v_U-ep,  length=M) #adaptive grid-forming
  ist <- evaluate_inverse_ST(t,w,gamma,v)

  if (ist$decreasing==0) {
    num_clus <- num_clus + 1
    z_l_endpoints[num_clus] <- ist$local_min
    z_u_endpoints[num_clus] <- ist$local_max
    num_grid_points[num_clus] <- ist$local_max_ind - ist$local_min_ind + 1 #num grid points where increasing
    c <- num_grid_points[num_clus]
    z_grid[num_clus,1:c] <- ist$grid[ist$local_min_ind:ist$local_max_ind]
    v_grid[num_clus,1:c] <- v[ist$local_min_ind:ist$local_max_ind]
  }
}

#last interval between -1/t & 0
v_L <- B[num_unique]
v_U <- 0
v <- seq(v_L+ep, v_U-ep,  length=M) #adaptive grid-forming
ist <- evaluate_inverse_ST(t,w,gamma, v) #grid,minf,ind_min
num_clus <- num_clus + 1
z_l_endpoints[num_clus]  <- ist$minf
z_u_endpoints[num_clus]  <- Inf
c <- M-ist$ind_min+1
num_grid_points[num_clus] <- c #num grid points where increasing
z_grid[num_clus,1:c] <- ist$grid[ist$ind_min:M]
v_grid[num_clus,1:c] <- v[ist$ind_min:M]

#the positive line: between 0 and inf
num_clus <- num_clus+1
z_l_endpoints[num_clus] <- -Inf
if (gamma >1) { #if upper half maximum may exceed 0
  v_U <- 1
  v<- seq(ep, v_U-ep,  length=M)
  ist <- evaluate_inverse_ST(t,w,gamma, v) #grid, maxz,ind_max
  # while (ist$ind_max==length(v)) {  #no need to deal with Inf's
  #   v_U <- 10*v_U
  #   v <- seq(ep, v_U-ep,  length=M)
  #   ist <- evaluate_inverse_ST(t,w,gamma, v)
  # }
  z_u_endpoints[num_clus] <- ist$maxz
  num_grid_points[num_clus] <- ist$ind_max #number of grid points in last interval
  z_grid[num_clus,1:ist$ind_max] <- ist$grid[1:ist$ind_max]
  v_grid[num_clus,1:ist$ind_max] <- v[1:ist$ind_max]
  
} else { #if gamma<1, this interval gives me the whole region m>0;
  #know that the function is strictly increasing in this case
  v_U <- M
  v <- seq(ep, v_U-ep,  length=M)
  ist <- evaluate_inverse_ST(t,w,gamma, v) #grid
  z_u_endpoints[num_clus] <- 0
  num_grid_points[num_clus] <- M #number of grid points in last interval
  z_grid[num_clus,1:M] <- ist$grid[1:M]
  v_grid[num_clus,1:M] <- v[1:M]
}

## Part II: inside the support
#Now stitch together all intervals where Stieltjes transform is real, and
#fill in the intervals where the Stieltjes transform is complex-valued
# First Stitch together the real ST
# pre-allocate grids
grid_length <- (num_clus-2)*M + sum(num_grid_points)
v <- matrix(0,grid_length,1)
grid <- matrix(0,grid_length,1)
K_hat <- num_clus-2
l_hat <- matrix(0,K_hat,1)
u_hat <- matrix(0,K_hat,1)


grid[1:num_grid_points[num_clus]] <- z_grid[num_clus, 1:num_grid_points[num_clus]]
v[1:num_grid_points[num_clus]] <- v_grid[num_clus, 1:num_grid_points[num_clus]]
#store the increasing intervals
ind <- num_grid_points[num_clus] #current index

#if gamma>1 need to set the u-endpoint to higher than 0
if (gamma>1) {
  z_u_endpoints[1] <- max(grid[1:num_grid_points[num_clus]])
}

#if gamma <1 need to add in the first interval
if (gamma < 1) {
  grid[(ind + 1):(ind+num_grid_points[1])] <- z_grid[1, 1:num_grid_points[1]]
  v[(ind + 1):(ind+num_grid_points[1])] <- v_grid[1, 1:num_grid_points[1]]
  ind <- ind + num_grid_points[1]
}

#function for ODE solving
MP_diff <- function(time,v_arr,param) { #transform to complex then to real
  v <- v_arr[1]+1i*v_arr[2]
  d <- (1/v^2 - gamma* sum( w*(t^(-1) + v)^(-2)))^(-1)
  list(c(Re(d),Im(d)))
}

#from now on, can handle all points in a uniform way
#have num_clus-2 pairs of intervals left: look at the support first, then
#outside
for (i in 2:(num_clus-1)) {
  #within the support
  #define the parameters of an interval in the support
  endpoint_1 <-  z_u_endpoints[i-1] #lower endpoint in f-space
  endpoint_2 <-  z_l_endpoints[i] #upper endpoint in f-space
  l_hat[i-1] <- endpoint_1
  u_hat[i-1] <-  endpoint_2
  
  #the grid within the i-th support interval
  grid_current <-  seq(endpoint_1,endpoint_2, length=M)

  #find starting point for ODE using iterative method:
  c0 <-  compute_esd_fp(t,w,gamma,ep, grid_current[1])
  v_start <- c0$v
  #find dual Stieltjes transform using ode
  ode_out <- ode(y=c(Re(v_start),Im(v_start)),times=grid_current,func=MP_diff,rtol=ep,atol=ep)
  v0 <- ode_out[,2]+1i*ode_out[,3]
  
  #set the output
  grid[(ind+1):(ind + M)] <- grid_current
  v[(ind+1):(ind + M)] <- v0
  ind <- ind + M
  
  #add in the corresponding interval outside the suppost
  grid[(ind + 1):(ind+num_grid_points[i])] <-  z_grid[i, 1:num_grid_points[i]]
  v[(ind + 1):(ind+num_grid_points[i])] <- v_grid[i, 1:num_grid_points[i]]
  ind <- ind + num_grid_points[i]
}

#retain the relevant part of the grid
good_ind <- (grid>0)&(grid<(1+sqrt(gamma))^2*max(t)*1.1)
grid <- grid[good_ind]
v <- v[good_ind]

#compute the usual ST from the dual ST
m <- 1/gamma*v- (1-1/gamma)/grid
density <- abs(1/pi*Im(m))

results <- list("grid" = grid, "density" = density, "m" = m, "v" = v,
  "K_hat" = K_hat,  "l_hat" = l_hat, "u_hat" = u_hat )
return(results)
}
