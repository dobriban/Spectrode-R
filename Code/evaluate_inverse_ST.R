evaluate_inverse_ST <- function(t,w,gamma,m) {
  #evaluate the inverse Stieltjes transform with spectrum t on grid m
  #Inputs
  # t - population eigenvalues >=0
  # gamma - aspect ratio p/n
  # m - real-valued grid to evaluate inverse ST
  # w - mixture weights >=0 default: w <- uniform
  #
  # Outputs list containing
  # f - inverse ST
  # minf - its minimum
  # ind_min - where it's attained
  # maxf, ind_max - same for max
  # decreasing - is the inverse decreaing on this interval
  # localmin, localmax - if f is not incresaing, the indices where the
  #   two inflexion points happen . Here we assume that m is an interval
  #   between two spikes in the population spectrum, so that it is
  #   known that the f goes down, up, down
  
  p <- length(t)
  num <- length(m)
  inv_ST <- function(m) - 1/m + gamma* sum(w*t/(1 + t*m))
  
  f <- matrix(0,num,1)
  decreasing <- 1 #start out decreasing
  local_min <- NaN
  local_max <- NaN
  local_min_ind <- NaN
  local_max_ind <- NaN
  
  for (i in 1:num) {
    f[i] <- inv_ST(m[i])
    
    #look for first inflexion point
    if ((i>1)&&(decreasing==1)) {
      if (f[i]>f[i-1]) {
        decreasing <- 0
        local_min <- f[i-1]
        local_min_ind <- i-1
      }
    }
    #if found first inflexion point, look for second inflexion point
    if ((i>1)&&(decreasing==0)&&(is.na(local_max))) {
      if (f[i]<f[i-1]) {
        local_max <- f[i-1]
        local_max_ind  <- i-1
      } else {
        if (i==num) { #need to give special treatment to the special case where f is increasing all the way
          local_max <- f[i]
          local_max_ind  <- i
        }
      }
    }
  }
  
  minf <- min(f[f>-Inf])
  ind_min <- which(f==minf)
  maxf <- max(f[f<Inf])
  ind_max <- which(f==maxf)
  
  results <- list(
    "grid" = f, "minf" = minf, "ind_min" = ind_min, "maxf" = maxf, "ind_max" = ind_max,
    "decreasing" = decreasing,  "local_min"=local_min, "local_max"=local_max,
    "local_min_ind"=local_min_ind, "local_max_ind"=local_max_ind
  )
  return(results)
}

