# This is a sampler for constraint Wishart distributions based on Hamiltonian Monte Carlo.
## For any questions or comments please refer to alexander.buchholz@ensae.fr
### We sample from Wishart distributions based on their Cholesky and eigenvalue decomposition.
### Arising hyper-parameters are adaptively tuned via Bayesian optimization based on Gaussian process regression

The functions needed for the simulation are in the folder functions and its subfolders respectively. 
In order to run the algorithm one has to exectute the scripts in the folder "scripts_simulations". 
There are three different scripts. One script for simulation of Wishart distributed matrices based on its Cholesky decomposition, one script based on the eigenvalue decomposition and one script based on the eigenvalue decomposition for sampling from Wishart distributions with eigenvalue constraints. 

In the scripts the parameters of the algorithm can be specified like the dimension of the matrix, the grid on which the hyper-parameters should be optimized, the resulting sample size and other parameters for the optimization procedure based on Gaussian process regression. 
