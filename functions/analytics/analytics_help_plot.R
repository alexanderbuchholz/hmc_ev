library(ggplot2)
library(xtable)
#setwd("~/DATA/masterarbeit/simulations/functions")
#setwd("/home/alex/Arbeitsfläche/exchange_MA/help_files")
#source(file= "function_MASTER.R", local=T)

#setwd("H:/masterarbeit/simulations/figures")
#l.input_all <- list.results
#acceptance_rate <- rep(0,length(l.input_all))
#big_M <- 10000

#for(i in 1:length(l.input_all)){
#  acceptance_rate[i] <- 1 - sum(unlist(l.input_all[[i]][7]))/big_M
#}
#hmc_final_simulation
#which.max(acceptance_rate)
#l.input_one_combination <- hmc_final_simulation
#l.input_one_combination <- l.input_all[[2]]

# obtain information of results
f.plot_results_simple_sigma <- function(l.input_one_combination, file_string){
  big_M <- length(l.input_one_combination[8][[1]]) # big_M <- 1000
  v <- l.input_one_combination[4][[1]] #  v <- 20
  s <- l.input_one_combination[5][[1]] # s <- 2
  epsilon <- round(l.input_one_combination[9][[1]],3) # 
  L <- l.input_one_combination[10][[1]] # 
  bounces <- l.input_one_combination[8][[1]] # s <- 2
  
  Sigma <- as.array(l.input_one_combination[3])[[1]] # Sigma <- matrix(c(2,1,1,2),2,2) 
  if(length(dim(Sigma))>2){Sigma <- matrix(c(1,0,0,1),2,2)}
  mc_results <- as.array(l.input_one_combination[2][1])[[1]]
  mc_results <- mc_results[,,1001:5000]
  
  sink(paste("output_latex_mean", "epsilon", gsub(".","_",as.character(epsilon), fixed=T), "L ", L, file_string,".txt", sep=""))
  print(xtable(f.analytics_mean(mc_results,v)))
  sink()
  cord_hmc <- f.analytics_coordinates(mc_results)
  
  yy <- rWishart(8000,v,Sigma)
  cord_wish <- f.analytics_coordinates(yy)
  
  pdf(paste("density_save_y_z_", "_epsilon_", gsub(".","_",as.character(epsilon), fixed=T), "_L_", L,file_string,".pdf", sep=""), width=15, height=8)
  g.left <- ggplot(cord_hmc, aes(x=y, y=z)) +  theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white")) +geom_point()+ geom_density2d(color="gray", size=1.5)  +   ggtitle(paste("Simulation results HMC \nand realized density"))  # Use hollow circles
  g.middle <- ggplot(cord_hmc, aes(x=y, y=z)) + theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white"))+geom_point()  + geom_density2d(data = cord_wish, aes(x= y, y= z), color="gray", size=1.5) + ggtitle(paste("Simulation results HMC \nand target density"))  # Use hollow circles
  g.right <- ggplot(cord_hmc[1:500,], aes(x=y, y=z)) +  theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white"))+ geom_path()+  ggtitle(paste("Simulation results HMC \n- jumps in target space"))  # Use hollow circles
  grid.arrange(g.left,g.middle, g.right, widths = c(1/3, 1/3,1/3)) 
  dev.off()
  
  #dev.new()
  pdf(paste("density_save_y_w_", "_epsilon_", gsub(".","_",as.character(epsilon), fixed=T), "_L_", L,file_string,".pdf", sep=""), width=15, height=8)
  g.left <- ggplot(cord_hmc, aes(x=y, y=w)) +  theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white")) +geom_point()+ geom_density2d(color="gray", size=1.5)  +   ggtitle(paste("Simulation results HMC \nand realized density"))  # Use hollow circles
  g.middle <- ggplot(cord_hmc, aes(x=y, y=w)) + theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white"))+geom_point()  + geom_density2d(data = cord_wish, aes(x= y, y= w), color="gray", size=1.5) + ggtitle(paste("Simulation results HMC \nand target density"))  # Use hollow circles
  g.right <- ggplot(cord_hmc[1:500,], aes(x=y, y=w)) +  theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white"))+ geom_path()+  ggtitle(paste("Simulation results HMC \n- jumps in target space"))  # Use hollow circles
  grid.arrange(g.left,g.middle, g.right, widths = c(1/3, 1/3,1/3)) 
  dev.off()
  
  pdf(paste("density_save_z_w_", "_epsilon_", gsub(".","_",as.character(epsilon), fixed=T), "_L_", L,file_string,".pdf", sep=""), width=15, height=8)
  g.left <- ggplot(cord_hmc, aes(x=z, y=w)) +  theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white")) +geom_point()+ geom_density2d(color="gray", size=1.5)  +   ggtitle(paste("Simulation results HMC \nand realized density"))  # Use hollow circles
  g.middle <- ggplot(cord_hmc, aes(x=z, y=w)) + theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white"))+geom_point()  + geom_density2d(data = cord_wish, aes(x= z, y= w), color="gray", size=1.5) + ggtitle(paste("Simulation results HMC \nand target density"))  # Use hollow circles
  g.right <- ggplot(cord_hmc[1:500,], aes(x=z, y=w)) +  theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white"))+ geom_path()+  ggtitle(paste("Simulation results HMC \n- jumps in target space"))  # Use hollow circles
  grid.arrange(g.left,g.middle, g.right, widths = c(1/3, 1/3,1/3)) 
  dev.off()
  
  
  pdf(file=paste("autocorr_first_coord", "epsilon", gsub(".","_",as.character(epsilon), fixed=T), "L_", L, file_string,".pdf", sep=""))
  acf(cord_hmc$y, main=paste("Autocorrelation first coordinate"))
  dev.off()
  
  pdf(file=paste("autocorr_second_coord", "epsilon", gsub(".","_",as.character(epsilon), fixed=T), "L_", L, file_string,".pdf", sep=""))
  acf(cord_hmc$w, main=paste("Autocorrelation second coordinate"))
  dev.off()
  
  pdf(file=paste("autocorr_third_coord", "epsilon", gsub(".","_",as.character(epsilon), fixed=T), "L_", L,file_string, ".pdf", sep=""))
  acf(cord_hmc$z, main=paste("Autocorrelation third coordinate"))
  dev.off()
  
  dev.new()
  #pdf(file=paste("density_first_coord_", "epsilon", epsilon, "L_", L,file_string, ".pdf", sep=""))
  m <- ggplot(cord_hmc, aes(x=y)) + labs(title="Histogram and density \nfirst coordinate")+ theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white")) + geom_histogram(aes(y=..density..),fill = "gray", binwidth=5)+ geom_density()+stat_function(fun=dchisq, args=list(df=v), colour="black", linetype="dashed")
  m
  ggsave(file=paste("density_first_coord_", "epsilon", gsub(".","_",as.character(epsilon), fixed=T), "L_", L,file_string, ".pdf", sep=""))
  dev.off()
  
  dev.new()
  #pdf(file=paste("density_second_coord_", "epsilon", epsilon, "L_", L,file_string, ".pdf", sep=""))
  m <- ggplot(cord_hmc, aes(x=w)) + labs(title="Histogram and density \nsecond coordinate")+ theme(plot.title =  element_text(size = rel(1.5)), panel.background = element_rect(fill = "white")) + geom_histogram(aes(y=..density..),fill = "gray", binwidth=5)+ geom_density()+stat_function(fun=dchisq, args=list(df=v), colour="black", linetype="dashed")
  m
  ggsave(file=paste("density_second_coord_", "epsilon", gsub(".","_",as.character(epsilon), fixed=T), "L_", L,file_string, ".pdf", sep=""))
  dev.off()
}

#for(i in 1:10){
#  print(i)
  #f.plot_results_simple_sigma (l.input_all[[i]])
#}

#f.plot_results_simple_sigma(hmc_final_simulation)
