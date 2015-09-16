f.objective_function_eigen <- function(output_list_hmc){
  mc_results <- as.array(output_list_hmc[[2]]) # getting hmc positions 
  L <- output_list_hmc[[1]][[10]] # getting L of simulation
  length_list <- dim(mc_results)[3] # number of simulations
  distance_jumps <- mc_results[,,2:length_list]-mc_results[,,1:(length_list-1)] # calculate jumps
  average_norm <- mean(apply(distance_jumps, 3, norm, "F")) # calculate the frobenius norm of jumps and their mean
  return(average_norm/sqrt(L)) # return the output
}

f.objective_function_cholesky <- function(output_list_hmc){
  mc_results <- as.array(output_list_hmc[[1]][[2]]) # getting hmc positions 
  L <- output_list_hmc[[1]][[10]] # getting L of simulation
  length_list <- dim(mc_results)[3] # number of simulations
  distance_jumps <- mc_results[,,2:length_list]-mc_results[,,1:(length_list-1)] # calculate jumps
  average_norm <- mean(apply(distance_jumps, 3, norm, "F")) # calculate the frobenius norm of jumps and their mean
  return(average_norm/sqrt(L)) # return the output
}
