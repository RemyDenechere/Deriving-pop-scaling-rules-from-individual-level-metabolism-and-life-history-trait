library(ggplot2)
# setwd("C:/Users/rden/OneDrive - Danmarks Tekniske Universitet/Ph.D/Project_1/Deriving-pop-scaling-rules-from-individual-level-metabolism-and-life-history-trait/shinny app")


rmax = function(A_0, b, M_0scl = "scaling"){
  n = 3/4
  a = 0.42
  epsR = 0.03
  M_0 = 0.01
  
  M = exp(seq(log(M_0), log(1000000), length.out=1000))
  
  if (M_0scl == "constant"){
    ratio = exp(10.299)* M
  } else if(M_0scl == "scaling") {
    ratio =  1612106 # based on elasmobtanches
  }
  
  rmax = A_0 * M^b * (1 - n) * M^(n-1) * ((1-a)*log(ratio) + log(epsR) )
  
  res = cbind.data.frame(rmax, M)
  return(res)  
}

plot_rmax = function(res){

  ggplot(res, aes(x = M, y = rmax))+
    geom_line(size=1)+ 
    xlab("Adult mass (g)")+
    ylab("r_max (1/yr)")+
    scale_x_log10(limits = c(0.001, 1000000))+
    scale_y_log10(limits = c(10^(-2), 10^2))+ 
    geom_vline(xintercept = 0.01, linetype = "dotted", 
               color = "blue", size=1)+
    annotate("text", 0.3 , 50, label = "Offspring size", color= "blue")
}