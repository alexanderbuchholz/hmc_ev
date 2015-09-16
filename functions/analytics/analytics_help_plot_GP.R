# plot.gaussian processes results


#setwd("~/DATA/masterarbeit/simulations/functions")
#setwd("/home/alex/Arbeitsfläche/exchange_MA/help_files")
#source(file= "function_MASTER.R", local=T)

library(reshape2)
library(ggplot2)
library(gridExtra)
library(akima)



f.plot_bo_convergence <- function(out_stat_eps_L_r, file_string){
  out_stat_eps_L_r$i <- 1:51
  p <- ggplot(out_stat_eps_L_r, aes(x=i, y=V2))
  p + geom_line()
  
  
  g.top <- ggplot(out_stat_eps_L_r, aes(x = i, y = V1)) +  geom_line() +  theme_bw() +  labs(y = "Epsilon - step size")+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank()) + ggtitle("Bayesian optimization via Gaussian process regression")
  g.middle <- ggplot(out_stat_eps_L_r, aes(x = i, y = V2)) +  geom_line() +  theme_bw() +  labs(y = "L - number of steps")+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())
  g.bottom <- ggplot(out_stat_eps_L_r, aes(x = i, y = r_i)) +  geom_line() +  theme_bw() +  labs(x = "Iteration i", y = "r_i - squared jump dist")
  
  #theme(      axis.text.x=element_blank(),      axis.ticks=element_blank(),    axis.title.x=element_blank(),      )
  ## plot graphs and set relative heights
  pdf(file=paste(file_string, "_BO_lines.pdf", sep=""))
  grid.arrange(g.top,g.middle, g.bottom, heights = c(1/3, 1/3, 1/3)) 
  dev.off()
  #g <- arrangeGrob(g.top,g.middle, g.bottom, heights = c(1/3, 1/3, 1/3)) 
  #ggsave(file=paste(file_string, "_BO_lines.pdf", sep=""),g)
  
  
  
  
  # interpolation
  fld <- with(out_stat_eps_L_r, interp(x = V1, y = V2, z = r_i, duplicate="mean"))
  
  df <- melt(fld$z, na.rm = TRUE)
  names(df) <- c("x", "y", "z")
  df$x <- fld$x[df$x]
  df$y <- fld$y[df$y]
  
  ggplot(data = df, aes(x = x, y = y, z = z)) +
    geom_tile(aes(fill = z)) +
    stat_contour() +  
    ggtitle("ESJD") +
    xlab("Epsilon - stepsize") +
    ylab("L - number of steps") +
    scale_fill_continuous(name = "Exp sq. \njum. dist.",
                          low = "white", high = "black") +
    theme(plot.title = element_text(size = 25, face = "bold"),
          legend.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.title.x = element_text(size = 20, vjust = 0.2),
          axis.title.y = element_text(size = 20, vjust = 1),
          legend.text = element_text(size = 10), panel.background = element_rect(fill = "white"))
  ggsave(file=paste(file_string, "_BO_surface.pdf", sep=""))
}