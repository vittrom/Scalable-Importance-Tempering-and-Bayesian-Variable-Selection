rm(list=ls())
setwd("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code")

source("functions_for_BVS.R")
library(foreach)
library(doParallel)

precompute_mult <- function(dataset, c_version="n", path = ".\\Data\\"){
  dt = read.table(paste(path, dataset, sep=""), header = TRUE)
  y = dt[,1]
  X = dt[,-1]
  n = nrow(X)
  p = ncol(X)
  
  prior_p_incl = 5/p
  c=n

  #Standardize
  X<-t(t(X)-colMeans(X))
  y<-y-mean(y)
  
  ### precompute matrices needed for samplers and create list to output ###
  XtX<-t(X)%*%X
  ytX<-t(y)%*%X
  yty<-sum(y^2)
  
  hyper_par<-list(
    n=n,y=y,X=X,XtX=XtX, ytX=ytX, yty=yty,
    prior_p_incl=prior_p_incl, c=c)
  return(hyper_par)
}

reps<-20
datasets <- c("dld.txt")#, "tgfb", "tgbf172")
T <- c(500)#, 1000, 30000)
burn_in = 0.1
c_version = c("n", "p2")
thin=2
full_cond=full_cond
n_HBS=1000
n_GS=1000

#Precompute Xtx etc


#Make combinations
combs = list()
for(i in 1:length(datasets)){
  dataset_name = datasets[i]
  hyper_par = precompute_mult(dataset_name)
  p = ncol(hyper_par$X)
  hyper_par_1 = hyper_par
  hyper_par_1$c = p^2
  combs[[length(combs) + 1]] = list(dataset=dataset_name, hyper_par=hyper_par, T=T[i], burn_in=ceiling(T[i] * burn_in), 
                                    reps=reps, thin=thin, full_cond=full_cond, n_HBS=n_HBS,n_GS=n_GS, p=p, c_vers="n")
  combs[[length(combs) + 1]] = list(dataset=dataset_name, hyper_par=hyper_par_1, T=T[i], burn_in=ceiling(T[i] * burn_in), 
                                    reps=reps, thin=thin, full_cond=full_cond, n_HBS=n_HBS,n_GS=n_GS, p=p, c_vers="p2")
}

#Settings for parallel run
cores = detectCores()
cl = makeCluster(2)
registerDoParallel(cl)

#Get results parallel run
res = foreach(comb=combs) %dopar% {
  source("functions_for_BVS.R")
  
  reps = comb$reps
  T= comb$T
  burn_in=comb$burn_in
  hyper_par = comb$hyper_par
  p=comb$p
  thin=comb$thin
  full_cond=comb$full_cond
  T_GS = comb$n_GS
  T_HBS = comb$n_HBS
  
  pip_HBS <- matrix(NA, reps, p)
  pip_GS<-matrix(NA, reps, p)
  pip_wTGS<-matrix(NA, reps, p)
  
  for(z in 1:reps){
    output_GS<-GS(p=p,hyper_par=hyper_par,T=T_GS,burn_in=burn_in,thin=thin)
    pip_GS[z,] <- output_GS$est_inclusion_probs
    
    output_wTGS<-wTGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in)
    pip_wTGS[z,] <- output_wTGS$est_inclusion_probs
    
    output_HBS<-HBS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in, full_cond=full_cond)
    pip_HBS[z,] <- output_HBS$est_inclusion_probs
  }
  
  res_par = list(dataset=comb$dataset, c_vers=comb$c_vers, pip_GS=pip_GS, pip_wTGS=pip_wTGS, pip_HBS=pip_HBS)
  save(res_par, file=paste(paste(".\\Results\\Simulation_Study\\dld\\res", comb$c_vers, sep = "_"), ".RData", sep=""))
  
}

#save(res, file = ".\\Results\\results_real_data.RData")