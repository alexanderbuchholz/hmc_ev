# mvgamma function

# Version 29.5.2015

library(psych)


## Function for MV Gamma
f.mvgamma <- function(p, x, log = TRUE)
{
  if(p>1)
  {
    out <- (p-1)/2*log(pi) + f.mvgamma(p-1, x, log = TRUE) + lgamma(x+(1-p)/2) 
  }
  else if (p == 1)
  {
    out <- lgamma(x) 
  }
  else stop("The integer p must exceed 0.")
  
  if(log == TRUE) # return the log form of multivariate gamma function
  {
    return(out)
  }
  else if(log  == FALSE) # return the multivariate gamma function
  {
    return(exp(out))
  }
  else stop("Wrong argument input!") 
}
