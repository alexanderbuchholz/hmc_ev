# function to run hmc for eigenvalues
#library(mixAK)
f.HMC_simulation_eigenvalues <- function(big_M, epsilon, L, v, s, lambda, intervall_ev=NULL){
  d <- s # set dimension
  # preallocate
  x = array(NA, dim=c(d,big_M)) # output 1: eigenvalues
  xx = array(NA, dim=c(c(d,d),big_M)) #  output 2: Wishart dist
  xxx = array(NA, dim=c(c(d,d),big_M)) #  output 3: rotation
  vec.bounce_counter <- rep(0,big_M) # output bounce counter
  vec.reject_counter <- rep(0,big_M) # output bounce counter
  
  inter.f.hmc = f.HMC_eigenvalues( epsilon, L, v, s,  lambda, intervall_ev)
  vec.reject_counter[1] = inter.f.hmc[[2]]
  vec.bounce_counter[1] = inter.f.hmc[[3]]
  
  x[,1] = inter.f.hmc[[1]]
  rotate <- f.random_rotation_incl_flip(s)
  xx[,,1] = rotate%*%diag(x[,1])%*%t(rotate)
  xxx[,,1] = rotate
  
  counter = 0 
  # loop over number of desired iterations
  for(i in 2:big_M){

    inter.f.hmc = f.HMC_eigenvalues( epsilon, L, v, s, x[,i-1], intervall_ev)
    inter = inter.f.hmc[[1]]
    
    vec.reject_counter[i] = inter.f.hmc[[2]]
    vec.bounce_counter[i] = inter.f.hmc[[3]]
    
    if(is.vector(inter)){
      if(vec.reject_counter[i]==0){
        # proposal has been accepted, generate random rotations
        x[,i] = inter
        rotate <- f.random_rotation_incl_flip(s)
        xx[,,i] = rotate%*%diag(x[,i])%*%t(rotate)
        xxx[,,i] = rotate
      }
      else{
        # proposal has been rejected, keep former values
        x[,i] = x[,i-1]
        xx[,,i] = xx[,,i-1]
        xxx[,,i] = xxx[,,i-1]
      }
    }
    else { 
      print("Error, no vector where needed!")
      x[,i] = x[,i-1]
      xx[,,i] = xx[,,i-1]
      xxx[,,i] = xxx[,,i-1]
      counter = counter + 1
    }
    
  }
  output <- list(x, xx, xxx, v, s, counter, vec.reject_counter, vec.bounce_counter,  epsilon, L)
  return(output)
}