# Hamiltonian function
# Version 30.5.2015
# set.seed(2)

library(base)
f.HMC = function (epsilon, L, v, s, Sigma, Sigma_inv, current_q){
  
  i.bounce_counter = 0
  i.reject_counter = 0
  
  current_q <- vech(current_q)
  q = vech(current_q)
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * f.gradient_u_hat(v, s, Sigma_inv, q) / 2
  # Alternate full steps for position and momentum
  for(i in 1:L){
    # Make a full step for the position
    q = q + epsilon * p
    ## Check if constraint violated:
      if(sum(diag(vech2full(q))<0) == 0){
        # Make a full step for the momentum, except at end of trajectory
        p = p - epsilon * f.gradient_u_hat(v, s,  Sigma_inv, q)
        #print(vech2full(p))
      }
    
      if(sum(diag(vech2full(q))<0) != 0){ # bounce
        I_d <- diag(1, dim(vech2full(q)))
        diag(I_d) = diag(I_d)*(diag(vech2full(q))<0) # change here!! 
        I_d_flat <- as.vector(vech(I_d))
        n <- I_d_flat
        p_flat <- as.vector(p)
        reflect_vec <- (t(p_flat)%*%n)*n
        p = (p - 2 * reflect_vec)
        
        i.bounce_counter = i.bounce_counter + 1
        
      }
    
  }
  q = q + epsilon * p
  # Make a half step for momev, p, Sigma, Sigma_invntum at the end.
  p = p - epsilon * f.gradient_u_hat(v, s, Sigma_inv, q) / 2
  
  if(sum(diag(vech2full(q))<0) != 0){
    q = current_q
   # print("Error! Despite bounces bad proposal!")
  }
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = f.u_hat(v, s, Sigma, Sigma_inv, current_q)
  current_K = sum(t(current_p)%*% current_p) / 2
  proposed_U = f.u_hat(v, s, Sigma, Sigma_inv, q)
  proposed_K = sum(t(p)%*% p) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  { out <- vech2full(q)
    out[!lower.tri(vech2full(q), diag=T)] <- 0
    output <- out # accept
    if(current_U==proposed_U){i.reject_counter <- i.reject_counter + 1}
  }
  else{
    out_alt <- vech2full(current_q)
    out_alt[!lower.tri(vech2full(current_q), diag=T)] <- 0
    output <- out_alt # reject
    i.reject_counter <- i.reject_counter + 1
    #print(output)
  }
  return(list(output, i.reject_counter, i.bounce_counter))
}