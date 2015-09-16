# simulation adaptive HMC cholesky decomp
# 16-9-15
# version 2
setwd("/home/alex/Dropbox/Master_Statistik/Masterarbeit/code_simulations/hmc_ev/functions")
#setwd("~/DATA/masterarbeit/simulations_cleaned/functions")
source(file= "function_MASTER.R", local=T)

working_directory = "/home/alex/Dropbox/Master_Statistik/Masterarbeit/code_simulations/test_results"
#working_directory = "~/DATA/masterarbeit/simulations_cleaned/simulations_HMC_eigen_BO/simulation_results/"
setwd(working_directory)


v = 20 # degrees of freedom
s = 2


diagonal <- diag(1,s,s)
low_tri <- lower.tri(matrix(0,s,s))*0
up_tri <- upper.tri(matrix(0,s,s))*0

Sigma = diagonal+low_tri+up_tri
Sigma_inv = solve(Sigma) #matrix(c(1,0,0,1),2,2)
l = t(chol(Sigma))

big_M = 20#20
epsilon_start = 0.4#0.4
L_start = 20#20
alpha_hyper <- 0.2
alpha2_hyper <- 4
k_hyper <- 100#100

Gamma_space <- matrix(c(0.01, 1, 5, 50), 2,2)


sigma_eta_hyper <- 0.01
number_runs <- 5

test_adaptive_hmc <- f.adaptive_hmc(Sigma, epsilon_start, L_start, v, s, big_M, number_runs, Gamma_space, k_hyper, alpha_hyper, sigma_eta_hyper, alpha2_hyper)

save(test_adaptive_hmc, file=paste("run_adaptive_hmc_runs_", number_runs,"_big_M", big_M,"_",s,"_", ".Rda", sep=""))

### Rerun with good parameters
# load stats file
load(file=paste("test_adaptive_hmc_bo_info",  "_v",v,"_s",s,".Rda", sep=""))

n_max <- which.max(out_stat_eps_L_r[,3])
epsilon <- as.numeric(out_stat_eps_L_r[n_max,1])
L <- as.numeric(out_stat_eps_L_r[n_max,2])

big_M = 10000

print("start final simulation")
hmc_final_simulation <- f.HMC_simulation(big_M, epsilon, L, v, s, Sigma, Sigma_inv, l)

save(hmc_final_simulation, file=paste("run_final_hmc_", "big_M", big_M,"_",s,"_epsilon",round(epsilon,3),"_L_",L , ".Rda", sep=""))
print("final simulation saved -- terminate")




