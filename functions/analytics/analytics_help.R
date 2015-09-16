# analytics
# version 31.5.2015
f.analytics_mean <- function(xx,v){
  big_M <- dim(xx)[3]
  i.dimension1 <- dim(xx)[1]
  test <- matrix(0,i.dimension1, i.dimension1)
  for(i in 1:big_M){
    test <- test + xx[,,i]
  }
  return(test/(big_M*v)  )
}

f.analytics_coordinates <- function(xx){
  big_M <- dim(xx)[3]
  y = rep(0,big_M)
  z = rep(0,big_M)
  w = rep(0,big_M)
  for(i in 1:big_M){
    y[i] = xx[1,1,i]
    z[i] = xx[2,1,i]
    w[i] = xx[2,2,i]
  }
  return(as.data.frame(cbind(y,z,w)))
}
