# Function that returns the the gradient of the negativ log potential of the eigenvalue density
# Input:
# lambda: vector with eigenvalues at which the density is to be evaluated
# v : int degrees of freedom
# s.s : int dimension of the target Wishart distribution, normally s.s = length(lambda)

# we need lamda to be ordered!
f.gradient_potential_eigenvalue_density <- function(lambda,v,s.s){
  mat.lambda <- diag(lambda)
  alpha = (v-s.s-1)/2
  # preallocate gradient
  gradient_out <- 0 
  # loop over entries for gradient
  for(i.int in 1:s.s){
    gradient_out[i.int] <- sum(1/(lambda[i.int]-lambda[-i.int]))
  }
  gradient_out <- gradient_out -1/2+alpha*(1/lambda)
  return(-gradient_out) # make it negative since we want the - log density
}

# test if gradient works 
#library(numDeriv)
#grad(f.potential_eigenvalue_density, lambda, method = "Richardson", side = NULL, method.args = list(),v,s)
#f.gradient_potential_eigenvalue_density(lambda,v,s)
#f.potential_eigenvalue_density(lambda,v,s)
