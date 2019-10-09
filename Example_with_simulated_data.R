rm(list=ls())
setwd("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code")
source("functions_for_BVS.R")

set.seed(12345)

###### SPECIFYING THE MODEL AND SIMULATING DATA #######
####
n<-100 # n.observations
p<-100 # n.covariates
c<-10^3 #prior hyperparameter for the covariance matrix
prior_p_incl<-5/p #prior prob of inclusion
hyper_par<-simulate_data(n=n,p=p,c=c,SNR=3,scenario=1, prior_p_incl=prior_p_incl) # creates list to be passed to samplers

###### RUNNING THE SAMPLERS #######
####
T<-50000 #n.iterations
reps<-2
burn_in<-5000
thin_GS<-2
output_GS<-GS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in,thin = thin_GS, vars_selected = c(3))
output_TGS<-TGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in, vars_selected = c(3))
output_wTGS<-wTGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in, vars_selected = c(3))
output_wTGS_true<-wTGS(p=p,hyper_par=hyper_par,T=reps*T,burn_in=burn_in, vars_selected = c(3))

##### PLOTTING THE OUTPUT ##########
####

plot_results_bvs(output_GS, output_TGS, output_wTGS, output_wTGS_true$est_inclusion_probs, 
                 variables = c(1,2,17), save = TRUE, height = 4.5)

