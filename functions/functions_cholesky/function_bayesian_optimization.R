# Adaptive HMC function
# 22-6-2015
#library(optimx)
#library(minqa)
library(nloptr)
#############################################################################
f.k_dist <- function(gamma_i, gamma_j, Sigma){
  return(exp(-0.5*as.matrix(gamma_i)%*%solve(Sigma)%*%t(as.matrix(gamma_j))))
}
# test 
#gamma_i = grid_epsilon_L[1,]
#gamma_j = grid_epsilon_L[4,]
#f.k_dist(gamma_i, gamma_j, Sigma)
#############################################################################
f.k_gamma <- function(new_gamma, old_gamma, Sigma){
  nr.rep <- dim(old_gamma)[1]
  vec.k <- rep(0,nr.rep)
  for(i.k in 1:nr.rep){
    vec.k[i.k] <- f.k_dist(t(new_gamma), old_gamma[i.k,], Sigma)
  }
  return(vec.k)
}
# test
# f.k_gamma(c(0.015,2), grid_epsilon_L, Sigma)
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
# test
# f.K_matrix(grid_epsilon_L, Sigma)
#############################################################
f.mu_i <- function(new_gamma, grid_epsilon_L, Sigma, sigma_eta, r_i, s_adaptive){
  k <- f.k_gamma(new_gamma, grid_epsilon_L, Sigma)
  K <- f.K_matrix(grid_epsilon_L, Sigma)
  return(t(k)%*% solve(K+diag(1,dim(K))*sigma_eta)%*%r_i*s_adaptive)
}
# test
# f.mu_i(t(new_gamma_start), grid_epsilon_L, Sigma, sigma_eta_hyper, r_i, 1)

#############################################################
f.sigma_i <- function(new_gamma, grid_epsilon_L, Sigma, sigma_eta){
  k_small <- f.k_dist(t(new_gamma), t(new_gamma), Sigma)
  k <- f.k_gamma(new_gamma, grid_epsilon_L, Sigma)
  K <- f.K_matrix(grid_epsilon_L, Sigma)
  return(k_small-t(k)%*% solve(K+diag(1,dim(K))*sigma_eta)%*%k)
}
# test
# f.sigma_i(t(new_gamma_start), grid_epsilon_L, Sigma, sigma_eta_hyper)
#############################################################
f.UCB <- function(new_gamma, grid_epsilon_L, Sigma, sigma_eta, r_i, k_hyper, s_adaptive){
  i <- dim(grid_epsilon_L)[1]
  delta <- 0.1
  d <- 2
  beta_i <- 2*log((((i+1)^((d/2)+2))*pi^2)/(3*delta))
  p_i <- max(i-k_hyper+1,1)^(-0.5)
  return(-(f.mu_i(new_gamma, grid_epsilon_L, Sigma, sigma_eta, r_i, s_adaptive)+p_i*beta_i^0.5*f.sigma_i(new_gamma, grid_epsilon_L, Sigma, sigma_eta)))
}
# test
# f.UCB(t(new_gamma_start), epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, 1)

# defining input parameters to test
#l.epsilon <- (1:5)/10#c(0.01,0.02,0.03)
#l.L <- 1 #(1:10)*5 #c(1,2,3)
#epsilon_L_i <- expand.grid(l.epsilon,l.L)
#grid_epsilon_l <- epsilon_L_i

#alpha_hyper <- 0.2
#s_adaptive <- 1

#new_gamma_start <- tail(epsilon_L_i, n=1)
#set.seed(11111)
#r_i <- runif( dim(epsilon_L_i)[1])+1
#k_hyper <- 2
#sigma_eta_hyper <- 1

# f.UCB(new_gamma_start, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, 1)
#Gamma_space <- matrix(c(0.1, 1, 1, 20),2,2)


f.bayesian_optimization <- function(new_gamma_start, Gamma_space, k_hyper, alpha_hyper, sigma_eta_hyper, r_i, epsilon_L_i, method , s_adaptive){
  
  # define the Gamma_space 
  upper_epsilon <- Gamma_space[2,1]
  lower_epsilon <- Gamma_space[1,1]
  
  upper_L <- Gamma_space[2,2]
  lower_L <- Gamma_space[1,2]
  
  Sigma <- as.matrix(diag(alpha_hyper^2*c((upper_epsilon-lower_epsilon)^2, (upper_L-lower_L)^2)))
  ###########################################
  if(method=="L-BFGS-B"){
    #print(f.UCB(new_gamma_start+0.00001, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive))
    return(optim(par=new_gamma_start, fn=f.UCB, gr = NULL, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive, method= "L-BFGS-B", lower=Gamma_space[1,],upper=Gamma_space[2,], control=list()))
  }
  else if(method=="direct"){
    #print(f.UCB(new_gamma_start+0.00001, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive))
    #Sigma1 <- Sigma
    #return(bobyqa(par=new_gamma_start, fn=f.UCB, lower=Gamma_space[1,],upper=Gamma_space[2,], epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive))
    #new_gamma_start <- as.numeric(as.vector(new_gamma_start))+0.1
    #optimx(new_gamma_start, fn=f.UCB,   lower=Gamma_space[1,],upper=Gamma_space[2,], control=list(all.methods=T), epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive)
    #optim(new_gamma_start, fn=f.UCB,  method="SANN",epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive  )
    # f.UCB(new_gamma_start, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive)
    return(direct(fn=f.UCB,   lower=Gamma_space[1,],upper=Gamma_space[2,],   control=list(maxeval = 100), grid_epsilon_L = epsilon_L_i, Sigma=Sigma, sigma_eta=sigma_eta_hyper, r_i =r_i, k_hyper=k_hyper, s_adaptive=s_adaptive))
    
    #direct(f.function, lower=-100, upper=100)
    #    bobyqa(x0=new_gamma_start, fn=f.UCB,   lower=Gamma_space[1,], upper=Gamma_space[2,],  grid_epsilon_L = epsilon_L_i, Sigma=Sigma, sigma_et=sigma_eta_hyper, r_i =r_i, k_hyper=k_hyper, s_adaptive=s_adaptive)
  }
  
  else if(method=="SANN"){
    return(optim(new_gamma_start+0.00001, fn=f.UCB, gr = NULL, epsilon_L_i, Sigma, sigma_eta_hyper, r_i, k_hyper, s_adaptive, method= "SANN", control=list(maxit=500)))
  }
}

# test
#f.bayesian_optimization(Gamma_space, k_hyper, alpha_hyper, sigma_eta_hyper, r_i, epsilon_L_i, "SANN")
