# function to run hmc 
f.HMC_simulation <- function(big_M, epsilon, L, v, s, Sigma, Sigma_inv, l){
  d <- dim(l)
  # preallocate
  x = array(NA, dim=c(d,big_M)) # output 1: Cholesky decomp
  xx = array(NA, dim=c(d,big_M)) #  output 2: Covar
  vec.bounce_counter <- rep(0,big_M) # output bounce counter
  vec.reject_counter <- rep(0,big_M) # output bounce counter
  
  inter.f.hmc = f.HMC( epsilon, L, v, s, Sigma, Sigma_inv, l)
  vec.reject_counter[1] = inter.f.hmc[[2]]
  vec.bounce_counter[1] = inter.f.hmc[[3]]
    
  x[,,1] = inter.f.hmc[[1]]
  xx[,,1] = x[,,1]%*%t(x[,,1])
  
  counter = 0 
  for(i in 2:big_M){
    is.natural <- function(x){ x>0 && identical(round(x), x)    }
    if(is.natural((i/500))){print(paste("Iteration", i))}
    
    inter.f.hmc = f.HMC( epsilon, L, v, s, Sigma, Sigma_inv, x[,,i-1])
    inter = inter.f.hmc[[1]]
    
    vec.reject_counter[i] = inter.f.hmc[[2]]
    vec.bounce_counter[i] = inter.f.hmc[[3]]
    
    if(is.matrix(inter)){
      x[,,i] = inter}
    else { 
      x[,,i] = x[,,i-1]
      counter = counter + 1
    }
    xx[,,i] = x[,,i]%*%t(x[,,i])
  }
  output <- list(x, xx, Sigma, v, s, counter, vec.reject_counter, vec.bounce_counter, epsilon, L)
  return(output)
}