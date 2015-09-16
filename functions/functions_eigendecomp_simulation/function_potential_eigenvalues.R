# library(mixAK)
# 
# 
# rRotationMatrix(1,2)
# 
# P <- rRotationMatrix(n=1, dim=5)
# print(P)
# round(P %*% t(P), 10)
# round(t(P) %*% P, 10)
# det(P)
# 
# n <- 10
# P <- rRotationMatrix(n=n, dim=5)
# for (i in 1:3){
#   cat(paste("*** i=", i, "\n", sep=""))
#   print(P[[i]])
#   print(round(P[[i]] %*% t(P[[i]]), 10))
#   print(round(t(P[[i]]) %*% P[[i]], 10))
#   print(det(P[[i]]))
# }


# A function that returns the potential of the eigenvalue density
#library(matrixcalc)
#library(psych)
#setwd("~/DATA/masterarbeit/simulations/functions")
#setwd("/home/alex/Arbeitsfläche/exchange_MA/help_files")
#source(file= "function_MASTER.R", local=T)
# case where the covariance matrix is equal to one
# v = 8 
# s = 5
# lambda <- 10+10*runif(s)
# lambda <- 1:3

# we need lamda to be ordered!


f.potential_eigenvalue_density <- function(lambda,v,s.s){
  mat.lambda <- diag(lambda)
  # lambda = eigenvalues a vector
  # need v = degrees of freedom
  # need s= number of eigenvalues
  
  #beta = log(2^(v*s/2)*det(Sigma)^(v/2)* f.mvgamma(s,v/2,FALSE))
  #prod_lambda <- 1
  # get the product of lambda's entries
  #for(it.i in 2:s.s){
  #  prod_lambda <- prod_lambda*prod(abs(lambda[it.i]-lambda[1:(it.i-1)]))
  #}
  # change for summation
  sum_lambda <- 0
  for(it.i in 2:s.s){
    sum_lambda <- sum_lambda+sum(log(abs(lambda[it.i]-lambda[1:(it.i-1)])))
  }
  
  alpha = (v-s.s-1)/2
  #beta = log(2^(v*s.s/2)* f.mvgamma(s.s,v/2,FALSE))
  #return(-log(exp(-beta)*det(mat.lambda)^alpha*exp(-1/2*tr(mat.lambda))*prod_lambda))
  beta = (v*s.s/2)*log(2)+f.mvgamma(s.s,v/2,T)
  
  #return(+beta -alpha*log(det(mat.lambda)) +1/2*tr(mat.lambda)-log(prod_lambda))
  return(+beta -alpha*log(det(mat.lambda)) +1/2*tr(mat.lambda)-sum_lambda)
}

#f.potential_eigenvalue_density(lambda,v,s)
