# function that returns the negative log potential for a given set of eigenvalues
f.potential_eigenvalue_density <- function(lambda,v,s.s){
  mat.lambda <- diag(lambda)
  sum_lambda <- 0
  
  for(it.i in 2:s.s){
    sum_lambda <- sum_lambda+sum(log(abs(lambda[it.i]-lambda[1:(it.i-1)])))
  }
  
  alpha = (v-s.s-1)/2
  beta = (v*s.s/2)*log(2)+f.mvgamma(s.s,v/2,T)
  return(+beta -alpha*log(det(mat.lambda)) +1/2*tr(mat.lambda)-sum_lambda)
}
