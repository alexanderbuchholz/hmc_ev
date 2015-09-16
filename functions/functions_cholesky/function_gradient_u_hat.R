# function gradient U hat

# Version 9.6.2015
library(OpenMx)

f.gradient_u_hat <- function(v,s, Sigma_inv, l){
  l = vech2full(l)
  l[!lower.tri(l, diag=T)] <- 0 ### ATTENTION ! 
  
  alpha = (v-s-1)/2
  m = dim(l)[1]
  T_mm <- f.vectorize_tranpose(m)
  I_m <- diag(1,m)
  I_mm <- diag(1,m^2)
  elimination_m <- elimination.matrix(m)
  t_elimination_m <- t(elimination_m)
  
  # PART 1
  derive_ll <- (I_mm + T_mm) %*% (l %x% I_m)
  derive_ll <- derive_ll + 0.00000001*I_mm # regularization 
  derive_ll_redu <- elimination_m %*% derive_ll %*% t_elimination_m
  
  vec_inv_derive_ll <- t(as.vector(t(solve(derive_ll_redu))))
  large_kroneckers_1 <- (elimination_m %x% elimination_m) + (elimination_m %x%  (elimination_m %*% T_mm))
  large_kroneckers_2 <- (I_m %x% T_mm %x% I_m ) %*% (I_mm %x% as.vector(I_m))
  part1 <- - vec_inv_derive_ll %*% large_kroneckers_1 %*% large_kroneckers_2 %*% t_elimination_m
  
  #large_kroneckers <- (I_mm %x% I_mm) %*% (I_m %x% as.vector(I_m) %x% I_m) + (I_mm %x% T_mm) %*% (I_m %x% as.vector(I_m) %x% I_m) %*% elimination_m
  
  # PART 2
  #derive_log_hl <- -alpha *as.vector(solve(l%*%t(l)))%*% derive_ll %*% t_elimination_m
  derive_log_hl <- -alpha *as.vector(solve(l%*%t(l)+ 0.00000001*I_m))%*% derive_ll %*% t_elimination_m
  # PART 3
  derive_tr_hl <- 0.5 *as.vector(Sigma_inv) %*% derive_ll %*% t_elimination_m
  
  half_vec_output <- part1 + derive_log_hl + derive_tr_hl
  #vec_output <- vech2full(half_vec_output)
  vec_output <- half_vec_output
  return(vec_output)
}
