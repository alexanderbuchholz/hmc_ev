# simulation Bayesian optimization eigen decomp
# 16-9-2015
# version 3


setwd("/home/alex/Dropbox/Master_Statistik/Masterarbeit/code_simulations/hmc_ev/functions")
#setwd("~/DATA/masterarbeit/simulations_cleaned/functions")
source(file= "function_MASTER.R", local=T)

working_directory = "/home/alex/Dropbox/Master_Statistik/Masterarbeit/code_simulations/test_results"
#working_directory = "~/DATA/masterarbeit/simulations_cleaned/simulations_HMC_eigen_BO/simulation_results/"
setwd(working_directory)

library(Matrix)

v = 30 # degrees of freedom
s = 10



diagonal <- diag(1,s,s)
low_tri <- lower.tri(matrix(0,s,s))*0
up_tri <- upper.tri(matrix(0,s,s))*0

Sigma = diagonal+low_tri+up_tri
############################################

lambda_start <- runif(s)

#####################################################



big_M = 30
epsilon_start = 0.4#0.4
L_start = 20#20
alpha_hyper <- 0.2
alpha2_hyper <- 4
k_hyper <- 100#100

Gamma_space <- matrix(c(0.001, 2, 1, 100), 2,2)


sigma_eta_hyper <- 0.1
number_runs <- 5
intervall_ev <- c(0,1)

output_adaptive_hmc <- f.adaptive_hmc_eigen(lambda_start, epsilon_start, L_start, v, s, big_M, number_runs, Gamma_space, k_hyper, alpha_hyper, sigma_eta_hyper, alpha2_hyper, intervall_ev)

save(output_adaptive_hmc, file=paste("results_adaptive_hmc_eigendecomp_", number_runs,"_big_M", big_M,"_",s,"_", ".Rda", sep=""))

### Rerun with good parameters
# load stats file
load(file=paste("stats_adaptive_hmc_bo_info",  "_v",v,"_s",s,".Rda", sep=""))

# obtain adequate hyper-parameters
max_index <- which.max(out_stat_eps_L_r$r_i)
epsilon <- as.numeric(out_stat_eps_L_r[max_index,1])
L <- as.numeric(out_stat_eps_L_r[max_index,2])

big_M = 10000

hmc_final_simulation <- f.HMC_simulation_eigenvalues(big_M, epsilon, L, v, s, lambda_start, intervall_ev)

save(hmc_final_simulation, file=paste("results_final_hmc_eigendecomp_", "big_M", big_M,"_",s,"_epsilon",round(epsilon,3),"_L_",L , ".Rda", sep=""))




