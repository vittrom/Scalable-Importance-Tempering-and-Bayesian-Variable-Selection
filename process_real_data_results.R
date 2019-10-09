rm(list=ls())
setwd("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code\\Results\\Simulation_Study\\dld\\")

plot_and_save <- function(data, title, vers, width=8, height=6.5){
  require(ggplot2)
  
  ggplot(data, aes(x=x, y=y)) + 
    geom_point() + labs(x="", y="", title = title) +
    geom_abline(slope=1, intercept=0) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste(title,"_", vers, ".pdf", sep=""), width = width, height = width, units = "cm")
  
}

scenario= "dld"
path = paste("D:/Studies/PhD/Year1/STAT548 - Qualifying course/Paper 1 - TGS/Code/Results/Simulation_Study/", scenario, "/", sep="")
files <- list.files(path, pattern="*.RData")

res_parsed = list()
count=1
for(res in files){
  full_path = paste(path, res, sep="")
  load(full_path)
  res_parsed[count] = list(res_par)
  count = count + 1
}

for(res in res_parsed){
  gs = res$pip_GS
  hbs = res$pip_HBS
  wtgs = res$pip_wTGS
  
  rows = nrow(gs); cols=ncol(gs)
  results_gs = matrix(NA, 1,2)
  results_hbs = matrix(NA, 1, 2)
  results_wtgs = matrix(NA, 1, 2)
  
  for(i in 1:rows){
    for(j in 1:cols){
      x = rep(gs[i, j], rows - 1)
      y = gs[-i, j]
      results_gs = rbind(results_gs, cbind(x, y))
      
      x = rep(hbs[i, j], rows - 1)
      y = hbs[-i, j]
      results_hbs = rbind(results_hbs, cbind(x, y))
      
      x = rep(wtgs[i, j], rows - 1)
      y = wtgs[-i, j]
      results_wtgs = rbind(results_wtgs, cbind(x, y))
    }
  }
  results_gs = as.data.frame(results_gs[-1,])
  results_hbs = as.data.frame(results_hbs[-1,])
  results_wtgs = as.data.frame(results_wtgs[-1,])
  
  plot_and_save(results_gs, "GS", res$c_vers)
  plot_and_save(results_hbs, "HBS", res$c_vers)
  plot_and_save(results_wtgs, "wTGS", res$c_vers)
  
}
