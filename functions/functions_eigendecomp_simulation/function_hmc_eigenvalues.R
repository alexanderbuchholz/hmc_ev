# Hamiltonian function for eigenvalue decomp
# Version 6.7.2015
# set.seed(2)

library(base)
f.HMC_eigenvalues = function (epsilon, L, v, s,  current_q, intervall_ev=NULL){
  
  i.bounce_counter = 0
  i.reject_counter = 0
  
  q = current_q
  #print(q)
  #p = matrix(rnorm(length(q),0,1),2,2) # independent standard normal variates
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * f.gradient_potential_eigenvalue_density(q, v, s) / 2
  # Alternate full steps for position and momentum
  for(i in 1:L){
    # Make a full step for the position
    q = q + epsilon * p
    
    ## Check if constraint violated:
    if(sum(q<0) == 0){
      # Make a full step for the momentum, except at end of trajectory
      p = p - epsilon * f.gradient_potential_eigenvalue_density(q, v, s)
      #print(vech2full(p))
    }
    
    if(!is.null(intervall_ev)){
      lower_ev = intervall_ev[1]
      upper_ev = intervall_ev[2]
      
      if((sum(q<lower_ev) != 0)|(sum(q>upper_ev) != 0)){ # bounce
        
        
        I_d <- rep(1, length(q))
        I_d = I_d*((q<lower_ev)| (q>upper_ev))# change here!! 
        I_d_norm = 1
        I_d_flat <- I_d
        n <- I_d_flat/I_d_norm
        p_flat <- as.vector(p)
        #reflect_vec <- (t(p_flat)%*%n)*n
        reflect_vec <- p_flat*n
        p = (p - 2 * reflect_vec)
        #print("Bounce!")
        #print(vech2full(p))
        #print(vech2full(q))
        #print(vech2full(q + epsilon * p))
        i.bounce_counter = i.bounce_counter + 1
        
      }
    } 
    
    else if(sum(q<0) != 0){ # bounce
      
      I_d <- rep(1, length(q))
      I_d = I_d*(q<0) # change here!! 
      I_d_norm = 1
      I_d_flat <- I_d
      n <- I_d_flat/I_d_norm
      p_flat <- as.vector(p)
      reflect_vec <- (t(p_flat)%*%n)*n
      p = (p - 2 * reflect_vec)
      #print("Bounce!")
      #print(vech2full(p))
      #print(vech2full(q))
      #print(vech2full(q + epsilon * p))
      i.bounce_counter = i.bounce_counter + 1
      
    }
    
  }
  q = q + epsilon * p
  # Make a half step for momev, p, Sigma, Sigma_invntum at the end.
  p = p - epsilon * f.gradient_potential_eigenvalue_density(q, v, s) / 2
  
  if(sum(q<0) != 0){
    #print(vech2full(q))
    q = current_q
     #print("Error! Despite bounces bad proposal!")
  }
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  s.s <- s
  current_U = f.potential_eigenvalue_density(current_q, v, s.s)
  if(!is.finite(current_U)){print("Error! gradient not finite!")}
  current_K = sum(t(current_p)%*% current_p) / 2
  proposed_U = f.potential_eigenvalue_density(q, v, s)
  proposed_K = sum(t(p)%*% p) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  { output <- q # accept
    #print(output)
  }
  else{
    out_alt <- current_q
    output <- out_alt # reject
    i.reject_counter <- i.reject_counter + 1
    #print(output)
  }
  #print("Number of bounces")
  #print(i.bounce_counter)
  #print("Number of rejections")
  #print(i.reject_counter)
  return(list(output, i.reject_counter, i.bounce_counter))
}

#current_q <- f.HMC_eigenvalues (epsilon, 100, v, s,  current_q)[[1]]
#print(current_q)
