setwd("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code")
source("functions_for_BVS.R")

library(foreach)
library(doParallel)
library(tidyverse)
set.seed(12345)

#### Simulation study comparison TGS vs GS ####
####
T<-50000 #n.iterations
burn_in<-5000
thin_GS<-2
n_scenarios<- 3
p_s<- c(100, 200, 1000)
n_s <- c(50, 200, 500)
SNR_s <- c(0.5, 1, 2, 3)
c<- 10^3
reps<-50

#Implement in parallel
cores = detectCores()
cl = makeCluster(cores - 3)
registerDoParallel(cl)

#Make combination of inputs
combs = list()
for(i in 1:n_scenarios){
  for(j  in 1:length(p_s)){
    for(k in SNR_s){
      combs[[length(combs) + 1]] = list(T=T, burn_in=burn_in, thin_GS=thin_GS, c=c, n=n_s[j], SNR=k,
                                        p=p_s[j], prior_p_incl=5/p_s[j], scenario=i, reps=reps)
    }
  }
}

#Source functions needed for parallel implementation and run the loop
res = foreach(comb = combs) %dopar% {
  source("functions_for_BVS.R", local = TRUE)
  hyper_par <- simulate_data(n=comb$n, p=comb$p, c=comb$c, SNR=comb$SNR, scenario=comb$scenario, prior_p_incl=comb$prior_p_incl)
  
  scenario=comb$scenario
  SNR=comb$SNR
  reps = comb$reps
  p=comb$p
  T=comb$T
  burn_in=comb$burn_in
  thin=comb$thin
  n=comb$n
  
  pip_TGS <- matrix(NA, reps, p)
  pip_GS<-matrix(NA, reps, p)
  pip_wTGS<-matrix(NA, reps, p)
  
  for(z in 1:reps){
    output_GS<-GS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in,thin=thin)
    pip_GS[z,] <- output_GS$est_inclusion_probs
    
    output_TGS<-TGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in)
    pip_TGS[z,] <- output_TGS$est_inclusion_probs
    
    output_wTGS<-wTGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in)
    pip_wTGS[z,] <- output_wTGS$est_inclusion_probs
  }
  
  #For table 1
  pip_var_GS <- median(apply(pip_GS, 2, var))
  pip_var_TGS <- median(apply(pip_TGS, 2, var))
  pip_var_wTGS<- median(apply(pip_wTGS, 2, var))
  eff_TGS_1 = pip_var_GS/pip_var_TGS
  eff_wTGS_1 = pip_var_GS/pip_var_wTGS
  
  #For table 2
  #Find columns to include
  cols_GS = which(apply(pip_GS >= 0.05, 2, sum) > 0)
  cols_TGS = which(apply(pip_TGS >= 0.05, 2, sum) > 0)
  cols_wTGS =  which(apply(pip_wTGS >= 0.05, 2, sum) > 0)
  if(length(cols_GS) == 0 || length(cols_TGS) == 0 || length(cols_wTGS) == 0){
    eff_TGS_2 = NA
    eff_wTGS_2 = NA
    list(scenario=scenario, SNR=SNR, p=p, n=n, eff_TGS_1=eff_TGS_1, eff_wTGS_1=eff_wTGS_1, eff_TGS_2=eff_TGS_2, eff_wTGS_2=eff_wTGS_2)
    next
  }
  cols_to_take = unique(c(cols_GS, cols_TGS, cols_wTGS))
  
  pip_var_GS <- exp(mean(log(apply(pip_GS[,cols_to_take], 2, var))))
  pip_var_TGS <- exp(mean(log(apply(pip_TGS[, cols_to_take], 2, var))))
  pip_var_wTGS<- exp(mean(log(apply(pip_wTGS[, cols_to_take], 2, var))))
  
  if(pip_var_GS == 0){
    eff_TGS_2 = NA
    eff_wTGS_2 = NA
  }else{
    eff_TGS_2 = pip_var_GS/pip_var_TGS
    eff_wTGS_2 = pip_var_GS/pip_var_wTGS
  }
  ## Output
  res_par = list(scenario=scenario, SNR=SNR, p=p, n=n, eff_TGS_1=eff_TGS_1, eff_wTGS_1=eff_wTGS_1, eff_TGS_2=eff_TGS_2, eff_wTGS_2=eff_wTGS_2)
  save(res_par, file=paste(paste(".\\Results\\Simulation_Study\\res", scenario, SNR, p, n, sep = "_"), ".RData", sep=""))
}

#Extract results
results_sim = as_tibble(t(sapply(res, unlist)))
save(results_sim, file=".\\Results\\results_simulation.RData")
