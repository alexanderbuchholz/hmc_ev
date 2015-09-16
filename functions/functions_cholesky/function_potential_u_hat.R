# potential u hat

# Version 29.5.2015
# correction 1.6.2015
library(matrixcalc)
library(psych)
f.u_hat = function(v, s, Sigma, Sigma_inv, l){
  l <- vech2full(l)
  l[!lower.tri(l, diag=T)] <- 0 ### ATTENTION ! 
  
  alpha = (v-s-1)/2
  beta = log(2^(v*s/2)*det(Sigma)^(v/2)* f.mvgamma(s,v/2,FALSE))
  m = dim(l)[1]
  
  # Change was effected here
  T_mm <- f.vectorize_tranpose(m)
  I_mm <- diag(1, m*m)
  I_m <- diag(1, m)
  #derive_ll <- (diag(1,m) %x% l)%*% T_mm + (l %x% diag(1,m))
  
  derive_ll <- (I_mm + T_mm) %*% (l %x% I_m)
  
  derive_ll <- elimination.matrix(m) %*% derive_ll %*% t(elimination.matrix(m))
  out <- -log(det(derive_ll ))-alpha*log(det(l%*%t(l)))+beta + 0.5*tr(Sigma_inv %*% l%*%t(l))
  return(out)
}

# test
#v = 5 
#s = 2
#Sigma = diag(c(1,1))
#Sigma_inv = Sigma
#l = matrix(c(1,1,0,1),2,2)
#f.u_hat(v, s, Sigma, Sigma_inv, l)