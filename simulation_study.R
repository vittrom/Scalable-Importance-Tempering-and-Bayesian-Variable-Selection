source("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code\\functions_for_BVS.R")

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

results_sim <- matrix(NA, n_scenarios*length(p_s), 2 * length(SNR_s))
results_sim_2 <- matrix(NA, n_scenarios*length(p_s), 2 * length(SNR_s))

for(i in 1:n_scenarios){
  print(paste("Scenario:", i))
  for(j in 1:length(p_s)){
    print(paste("Dimension:", j))
    prior_p_incl<-5/p_s[j] #prior prob of inclusion
    for(k in SNR_s){
      print(paste("SNR:", k))
      hyper_par<-simulate_data(n=n_s[j],p=p_s[j],c=c,SNR=k,scenario=i) # creates list to be passed to samplers
      pip_TGS <- matrix(NA, reps, p_s[j])
      pip_GS<-matrix(NA, reps, p_s[j])
      pip_wTGS<-matrix(NA, reps, p_s[j])
      for(z in 1:reps){
        print(paste("Iter:", z))
        output_GS<-GS(p=p_s[j],hyper_par=hyper_par,T=T,burn_in=burn_in,thin = thin_GS)
        pip_GS[z,] <- output_GS$est_inclusion_probs
        
        output_TGS<-TGS(p=p_s[j],hyper_par=hyper_par,T=T,burn_in=burn_in)
        pip_TGS[z,] <- output_TGS$est_inclusion_probs
        
        output_wTGS<-wTGS(p=p_s[j],hyper_par=hyper_par,T=T,burn_in=burn_in)
        pip_wTGS[z,] <- output_wTGS$est_inclusion_probs
      }
      
      #For table 1
      pip_var_GS <- median(apply(pip_GS, 2, var))
      pip_var_TGS <- median(apply(pip_TGS, 2, var))
      pip_var_wTGS<- median(apply(pip_wTGS, 2, var))
      results_sim[length(p_s) * (i - 1) +  j, which(SNR_s == k)] <- pip_var_GS/pip_var_TGS
      results_sim[length(p_s) * (i - 1) +  j, which(SNR_s == k) + length(SNR_s)] <- pip_var_GS/pip_var_wTGS
      
      #For table 2
      #Find columns to include
      cols_GS = which(apply(pip_GS >= 0.05, 2, sum) > 0)
      cols_TGS = which(apply(pip_TGS >= 0.05, 2, sum) > 0)
      cols_wTGS =  which(apply(pip_wTGS >= 0.05, 2, sum) > 0)
      if(length(cols_GS) == 0){
        next
      }
      cols_to_take = unique(c(cols_GS, cols_TGS, cols_wTGS))
      
      pip_var_GS <- mean(apply(pip_GS[,cols_to_take], 2, var))
      pip_var_TGS <- mean(apply(pip_TGS[, cols_to_take], 2, var))
      pip_var_wTGS<- mean(apply(pip_wTGS[, cols_to_take], 2, var))
      
      if(pip_var_GS == 0){
        next
      }else{
        results_sim_2[length(p_s) * (i - 1) +  j, which(SNR_s == k)] <- pip_var_GS/pip_var_TGS
        results_sim_2[length(p_s) * (i - 1) +  j, which(SNR_s == k) + length(SNR_s)] <- pip_var_GS/pip_var_wTGS
      }
      
    }
  }
}