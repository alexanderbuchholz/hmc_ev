## function to sample form cwish according to everson 2000

library(MASS)
f.constraint_chisq <- function(deg_free, d_constraint){
  return(qchisq((runif(1)*pchisq(d_constraint,deg_free)),deg_free))
}

#f.constraint_chisq(1, 1)
f.constraint_wish <- function(deg_free,p_size,N_runs, D_constraint){
  break_counter = 0
  arry_results = array(0,dim=c(p_size,p_size,N_runs))
  successful_tries = 0
  i = 1
  # while loop in order to complete tries
  while(successful_tries < (N_runs) ){
    
    # print something that lets the use know where we are
    if((break_counter %% 1000) ==0){
      print("Already too many breaks!")
      print(break_counter)
    }
    # set criteria for too many iterations
    if(break_counter>1000000) break
    
    # loop over matrices
    for(t in 1:p_size){
      # initialize
      if(t==1){
        u11 <- f.constraint_chisq(deg_free, D_constraint[1,1])
        T_t <- as.matrix(sqrt(u11))
      }
      # loop over other t 
      else{
        # necessary for breaking condition
        Delta_t_1 <- D_constraint[1:(t-1),1:(t-1)]-(T_t[1:(t-1),1:(t-1)]%*%t(T_t[1:(t-1),1:(t-1)]))
        # Specify necessary start values
        U_tt = f.constraint_chisq(deg_free, D_constraint[t,t])
        X_star_t = rchisq(1,df = (deg_free-t+1))
        Z_star_t = mvrnorm(1,  mu= rep(0,(t-1)), Sigma = diag(1,(t-1)))
        
        # calculate the values for j
        
        B_t <- (Z_star_t^2)/(t(Z_star_t)%*%Z_star_t+X_star_t )
        Z_t <- sign(Z_star_t)*(sqrt(B_t*U_tt))
        X_t <- U_tt-(t(Z_t)%*%Z_t)
        lambda_t = D_constraint[t,t] - U_tt- (t(Z_t)%*%t(T_t)%*%solve(Delta_t_1)%*%T_t%*%Z_t)
        # check if we have to break
        if( (lambda_t<0) | (X_t<0) ) break_counter = break_counter +1 
        if((lambda_t<0) | (X_t<0)) break
        # else construct new matrix
        
        T_t <- cbind(rbind(T_t, t(Z_t)), c(rep(0,(t-1)),sqrt(X_t)))
        
        if(t==p_size){
          # save results
          successful_tries = successful_tries + 1
          print(successful_tries)
          print("Another successfull try!")
          arry_results[,,successful_tries] <- T_t%*%t(T_t)
        }
        
      }
      
      
    }
    
  }
  return(list(arry_results, successful_tries, break_counter))
}


