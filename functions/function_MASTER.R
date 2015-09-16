# functions masterfile: contains packages to load 

# Version 29.5.2015
library(matrixcalc)
library(psych)
library(ggplot2)


# Load functions for the 
source("functions_cholesky/function_gradient_u_hat.R")
source("functions_cholesky/function_mvgamma.R")
source("functions_cholesky/function_potential_u_hat.R")
source("functions_cholesky/function_vec_transpose.R")
source("functions_cholesky/function_hmc.R")
source("functions_cholesky/function_hmc_simulation.R")
source("functions_cholesky/function_tune_hmc.R")
source("functions_cholesky/function_bayesian_optimization.R")
source("functions_cholesky/function_objective_maximization_adaptive_hmc.R")
source("functions_cholesky/function_adaptive_hmc.R")


source("analytics/analytics_help_plot_GP.R")
source("analytics/analytics_help_plot.R")
source("analytics/analytics_help.R")

## function for eigen decomp
source("functions_eigendecomp_simulation/function_gradient_potential_eigenvalues.R", chdir=T)
source("functions_eigendecomp_simulation/function_hmc_eigenvalues.R", chdir=T)
source("functions_eigendecomp_simulation/function_hmc_simulation_eigenvalues.R", chdir=T)
source("functions_eigendecomp_simulation/function_potential_eigenvalues.R", chdir=T)
source("functions_eigendecomp_simulation/function_random_rotation_matrix.R", chdir=T)
source("functions_eigendecomp_simulation/function_bayesian_optimization_eigen_decomp.R", chdir=T) 
source("functions_eigendecomp_simulation/function_grid_hmc_eigen.R", chdir=T) 


