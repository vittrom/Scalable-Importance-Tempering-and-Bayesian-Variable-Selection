setwd("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code\\Results\\Simulation_Study\\")
library(tidyverse)

myFun <- function(files, path) {
  load(paste(path, files, sep=""))
  res = as.data.frame(res_par)
  res
}

scenario= "scenario_1"
path = paste("D:/Studies/PhD/Year1/STAT548 - Qualifying course/Paper 1 - TGS/Code/Results/Simulation_Study/", scenario, "/", sep="")

files <- list.files(path, pattern="*.RData")

results = as_tibble(do.call(rbind.data.frame, lapply(files, myFun, path))) %>% 
  arrange(p, n)
results
