p_size <- 2
deg_free <- 4
N_runs <- 1000
D_constraint=diag(2,p_size)

test_results <- f.constraint_wish(deg_free,p_size,N_runs, D_constraint)
a <- apply(test_results[[1]],3,function(x){eigen(x)$value})

hist(a)

install.packages("tlnise")
library(tlnise)
library(mvtnorm)
x <- rmvtnorm(10,c(0,0),diag(1,2))  ## Second level
y <- rnorm(10, x)  ## First level means
out <- tlnise(Y = y, V = rep(1, 10), w = rep(1, 10), seed = 1234)
out$B0
