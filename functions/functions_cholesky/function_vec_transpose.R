# function that creates the permuation matrix for vector operations
# 11-5-2015

# Version 29.5.2015

f.vectorize_tranpose <- function(m){
  # Attention! Works only for symmetric matrices!
  n = m # still need solution for non quadratic matrics!
  d = m*n
  i = c(1:d)
  
  Tmn = matrix(0,d,d)
  rI = 1+m*(i-1)-(m*n-1)*floor((i-1)/n)
  diag(Tmn) = 1
  Tmn = Tmn[rI,i]
  return(Tmn)
  }