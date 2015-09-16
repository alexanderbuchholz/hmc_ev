# Adaptive HMC function
# 22-6-2015
library(nloptr)
#############################################################################
f.k_dist <- function(gamma_i, gamma_j, Sigma){
  return(exp(-0.5*as.matrix(gamma_i)%*%solve(Sigma)%*%t(as.matrix(gamma_j))))
}
#############################################################################
f.k_gamma <- function(new_gamma, old_gamma, Sigma){
  nr.rep <- dim(old_gamma)[1]
  vec.k <- rep(0,nr.rep)
  for(i.k in 1:nr.rep){
    vec.k[i.k] <- f.k_dist(t(new_gamma), old_gamma[i.k,], Sigma)
  }
  return(vec.k)
}
#############################################################################

f.K_matrix <- function(grid_epsilon_L, Sigma){
  i <- dim(grid_epsilon_L)[1]
  K <- matrix(0,i,i)
  for(i.col in 1:i){
    for(i.row in 1:i){
      K[i.row, i.col] <- f.k_dist(grid_epsilon_L[i.row,],grid_epsilon_L[i.col,] , Sigma)
    }
  }
  return(K)
  
}
#############################################################
f.mu_i <- function(new_gamma, grid_epsilon_L, Sigma, sigma_eta, r_i, s_adaptive){
  k <- f.k_gamma(new_gamma, grid_epsilon_L, Sigma)
  K <- f.K_matrix(grid_epsilon_L, Sigma)
  return(t(k)%*% solve(K+diag(1,dim(K))*sigma_eta)%*%r_i*s_adaptive)
}
#############################################################
f.sigma_i <- function(new_gamma, grid_epsilon_L, Sigma, sigma_eta){
  k_small <- f.k_dist(t(new_gamma), t(new_gamma), Sigma)
  k <- f.k_gamma(new_gamma, grid_epsilon_L, Sigma)
  K <- f.K_matrix(grid_epsilon_L, Sigma)
  return(k_small-t(k)%*% solve(K+diag(1,dim(K))*sigma_eta)%*%k)
}
#############################################################
f.UCB <- function(new_gamma, grid_epsilon_L, Sigma, sigma_eta, r_i, k_hyper, s_adaptive){
  i <- dim(grid_epsilon_L)[1]
  delta <- 0.1
  d <- 2
  beta_i <- 2*log((((i+1)^((d/2)+2))*pi^2)/(3*delta))
  p_i <- max(i-k_hyper+1,1)^(-0.5)
  return(-(f.mu_i(new_gamma, grid_epsilon_L, Sigma, sigma_eta, r_i, s_adaptive)+p_i*beta_i^0.5*f.sigma_i(new_gamma, grid_epsilon_L, Sigma, sigma_eta)))
}
###########################################################
f.bayesian_optimization <- function(new_gamma_start, Gamma_space, k_hyper, alpha_hyper, sigma_eta_hyper, r_i, epsilon_L_i, method , s_adaptive){
  
  # define the Gamma_space 
  upper_epsilon <- Gamma_space[2,1]
  lower_epsilon <- Gamma_space[1,1]
  
  upper_L <- Gamma_space[2,2]
  lower_L <- Gamma_space[1,2]
  
  Sigma <- as.matrix(diag(alpha_hyper^2*c((upper_epsilon-lower_epsilon)^2, (upper_L-lower_L)^2)))
  ###########################################
  if(method=="L-BFGS-B"){
    return(optim(par=new_gamma_start, fn=f.UCB, gr = NULL, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive, method= "L-BFGS-B", lower=Gamma_space[1,],upper=Gamma_space[2,], control=list()))
  }
  else if(method=="direct"){
    return(direct(fn=f.UCB,   lower=Gamma_space[1,],upper=Gamma_space[2,],   control=list(maxeval = 100), grid_epsilon_L = epsilon_L_i, Sigma=Sigma, sigma_eta=sigma_eta_hyper, r_i =r_i, k_hyper=k_hyper, s_adaptive=s_adaptive))
  }
  
  else if(method=="SANN"){
    return(optim(new_gamma_start+0.00001, fn=f.UCB, gr = NULL, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive, method= "SANN", control=list(maxit=500)))
  }
}
