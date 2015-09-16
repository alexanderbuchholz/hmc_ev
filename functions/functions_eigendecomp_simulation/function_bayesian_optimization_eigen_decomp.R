# adaptive hmc for eigen decom
# this function runs the bayesian optimization of the hyper-parameters
# function bayesian optimization

f.adaptive_hmc_eigen <- function(lambda_start, epsilon_start, L_start, v, s, big_M, number_runs, Gamma_space, k_hyper, alpha_hyper, sigma_eta_hyper, alpha2_hyper, intervall_ev=NULL){
  # Define starting values
  lambda <- lambda_start
  list.results <- list()
  L = L_start
  epsilon = epsilon_start
  r_i <- NULL
  epsilon_L_i <- as.data.frame(t(c(epsilon, L)))
  new_gamma_start <- epsilon_L_i
  s_adaptive = 1
  
  # Iterate over 
  for(i in 1:number_runs){
    # run simulation and store results
    test_output <- f.HMC_simulation_eigenvalues(big_M, epsilon, L, v, s, lambda, intervall_ev)
    print("Rejections")
    print(test_output[[7]])
    print("Bounces")
    print(test_output[[8]])
    list.results <- c(list.results, test_output)
    
    # prepare data for maximization
    int.value_objective_funtion <- f.objective_function_eigen(test_output)
    print(paste("Value objective function", int.value_objective_funtion))
    
    if((int.value_objective_funtion > max(r_i)) & (int.value_objective_funtion>0)){ s_adaptive = alpha2_hyper/int.value_objective_funtion}
    # augment data
    r_i <- rbind(r_i, int.value_objective_funtion)
    
    # choose if we optimize or not
    u = runif(1)
    i <- dim(r_i)[1]
    p_i <- max(i-k_hyper+1,1)^(-0.5)
    
    if(u<p_i){
      print("Start optimization")
      next_gamma_opt <- f.bayesian_optimization(new_gamma_start, Gamma_space, k_hyper, alpha_hyper, sigma_eta_hyper, r_i, epsilon_L_i, "direct", s_adaptive)
      print("Optimization finished")
    }  
    
    # filter out negative values of function
    if(sum((next_gamma_opt$par<0))){print("Negative result from Optimization!")}
    
    next_gamma <- (next_gamma_opt$par)*(next_gamma_opt$par>0)+ (Gamma_space[1,]+apply(Gamma_space, 2, mean)*2*runif(2))*(next_gamma_opt$par<0)
    next_gamma[2]  <- ceiling(next_gamma[2])
    new_gamma_start <- next_gamma
    epsilon_L_i <- rbind(epsilon_L_i, next_gamma)
    print("Next set of hyperparameters")
    print(next_gamma)      
    epsilon <- next_gamma[1]
    L <- next_gamma[2]
    
    
    
  }
  
  # run algorithm a last time 
  test_output <- f.HMC_simulation_eigenvalues(big_M, epsilon, L, v, s, lambda, intervall_ev)
  print("Rejections")
  print(test_output[[7]])
  print("Bounces")
  print(test_output[[8]])
  list.results <- c(list.results, test_output)
  
  # return value of the objective function
  int.value_objective_funtion <- f.objective_function_eigen(test_output)
  print(paste("Value objective function", int.value_objective_funtion))
  
  # combine results in data frame
  r_i <- rbind(r_i, int.value_objective_funtion)
  
  out_stat_eps_L_r <- cbind(epsilon_L_i, r_i)
  print(out_stat_eps_L_r)
  # save results to local file
  save(out_stat_eps_L_r, file=paste("stats_adaptive_hmc_bo_info",  "_v",v,"_s",s,".Rda", sep=""))
  return(list.results)
}