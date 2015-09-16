# Hamiltonian function for eigenvalue decomp
# Version 6.7.2015
# set.seed(2)

library(base)
f.HMC_eigenvalues = function (epsilon, L, v, s,  current_q, intervall_ev=NULL){
  
  i.bounce_counter = 0
  i.reject_counter = 0
  
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * f.gradient_potential_eigenvalue_density(q, v, s) / 2
  # Alternate full steps for position and momentum
  for(i in 1:L){
    # Make a full step for the position
    q = q + epsilon * p
    
    ## Check if constraint violated, if not, normal procedure
    if(sum(q<0) == 0){
      # Make a full step for the momentum, except at end of trajectory
      p = p - epsilon * f.gradient_potential_eigenvalue_density(q, v, s)
      #print(vech2full(p))
    }
    
    # include eigenvalue constraint
    if(!is.null(intervall_ev)){
      lower_ev = intervall_ev[1]
      upper_ev = intervall_ev[2]
      # Bounce procedure with ev constraints
      if((sum(q<lower_ev) != 0)|(sum(q>upper_ev) != 0)){ # bounce
        
        # routine for bouncing of the particle
        I_d <- rep(1, length(q))
        I_d = I_d*((q<lower_ev)| (q>upper_ev))
        n <- I_d
        p_flat <- as.vector(p)
        reflect_vec <- p_flat*n # select the components that have to be reflected
        p = (p - 2 * reflect_vec)
        i.bounce_counter = i.bounce_counter + 1
        
      }
    } 
    # Bounce procedure with out ev constraints, only one lower bound
    else if(sum(q<0) != 0){ # bounce
      
      I_d <- rep(1, length(q))
      I_d = I_d*(q<0) 
      n <- I_d
      p_flat <- as.vector(p)
      reflect_vec <- (t(p_flat)%*%n)*n
      p = (p - 2 * reflect_vec)
      i.bounce_counter = i.bounce_counter + 1
      
    }
    
  }
  q = q + epsilon * p
  # Make a half step for momev, p, Sigma, Sigma_invntum at the end.
  p = p - epsilon * f.gradient_potential_eigenvalue_density(q, v, s) / 2
  
  # check if despite bounces there is a bad proposal
  if(sum(q<0) != 0){
    q = current_q
     #print("Error! Despite bounces bad proposal!")
  }
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  s.s <- s
  current_U = f.potential_eigenvalue_density(current_q, v, s.s)
  # error if there are numerical problems
  if(!is.finite(current_U)){print("Error! gradient not finite!")}
  current_K = sum(t(current_p)%*% current_p) / 2
  proposed_U = f.potential_eigenvalue_density(q, v, s)
  proposed_K = sum(t(p)%*% p) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){ 
    output <- q # accept
    # check if the bounce procedure failed
    if(current_U==proposed_U){i.reject_counter <- i.reject_counter + 1}
  }
  else{
    output <- current_q # reject
    i.reject_counter <- i.reject_counter + 1
  }
  return(list(output, i.reject_counter, i.bounce_counter))
}
