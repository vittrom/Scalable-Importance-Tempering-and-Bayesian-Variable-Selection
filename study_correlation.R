setwd("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code")

source("tgs_normal.R")
library(matrixcalc)
#Correlation study for more than 2 dimensions
create_moments <- function(scenario, rho, d){
  mu = rep(0, d)
  sigma = diag(d)
  if(scenario == 1){
    for(i in 1:(d/2)){
      sigma[2*i - 1, 2*i] = rho
      sigma[2*i, 2*i - 1] = rho
    }
  }
  if(scenario == 2){
    sigma[which(sigma== 0)] = rho
  }
  if(scenario == 3){
    sigma[which(sigma== 0)] = -rho/(d-1)
  }
  return(list(mu=mu, sigma=sigma))
}

compute_var_est <- function(x, weights=NULL, method="NULL", corr_lag=200){
  
  h = x #- mean(results$x[,1])
  corr = acf(h, lag.max=corr_lag, plot = FALSE, type = "correlation")
  if(method=="TGS"){
    f_x =  (h * weights)
    corr = acf(f_x, lag.max = corr_lag, plot=FALSE, type="correlation")
  }
  var_h = mean(h^2)
  if(method=="TGS"){
    var_h = mean(f_x/sum(weights))^2 #sum((f_x^2 * results$weights^2)/length(f_x)
  }
  var_est = var_h * (1 + 2 * sum(corr$acf))
  return(var_est)
}

#Scenario 1: pairwise correlation
#Scenario 2: all positive correlation
#Scenario 3: all negative correlation

set.seed(12345)
d=2; burn_in=0
rho = c(0, 0.9, 0.99, 0.999, 0.9999)
reps = length(rho); n_scenarios = 2#3

tempering = 1-rho^2
var_TGS = matrix(NA, nrow = reps, ncol = n_scenarios)
var_GS =  c() #matrix(NA, nrow = reps, ncol = n_scenarios)
var_TGS_t =  matrix(NA, nrow = reps, ncol = n_scenarios)
lag_corr = 100; shape=5
iter = 1000
start_x = rep(3, d)

var_w_tgs = matrix(NA, nrow=reps, ncol=n_scenarios)
var_w_tgs_t = matrix(NA, nrow=reps, ncol=n_scenarios)

repeated = 100
for(j in 2:n_scenarios){
  for(dim in 1:reps){
    extra_pars = create_moments(scenario=j, rho[dim], d)
    
    # res_TGS = TGS(start_x = start_x, d=d, T=iter, burn_in=burn_in, tempering=tempering[dim], cond_distr = cond_distr_full_normal,
    #           single_cond = sample_g_cond_normal, extra_pars = extra_pars, version = "mixed")
    # var_TGS[dim, j] = mean(apply(res_TGS$x, 2, compute_var_est, res_TGS$weights, "TGS", lag_corr))
    
    # for(k in 1:repeated){
      res_GS = GS(start_x = start_x, d=d, T=iter, burn_in=burn_in, cond_moments = cond_moments_normal,
                  extra_pars = extra_pars, type = "deterministic")
      h_gs = iter * apply(res_GS$x, 2, var)
      var_GS = c(var_GS, h_gs)
      #var_GS[dim, j] = mean(apply(res_GS$x, 2, compute_var_est, NULL, "NULL",lag_corr))
    #}
    # res_GS = GS(start_x = start_x, d=d, T=iter, burn_in=burn_in, cond_moments = cond_moments_normal,
    #             extra_pars = extra_pars, type = "deterministic")
    # var_GS[dim, j] = mean(apply(res_GS$x, 2, compute_var_est, NULL, "NULL",lag_corr)) 

    # res_TGS_t = TGS(start_x = start_x, d=d, T=iter, burn_in=burn_in, tempering=tempering[dim], cond_distr = cond_distr_full_t,
    #                 single_cond = sample_g_cond_t, extra_pars = extra_pars, version = "t", shape=shape)
    #var_TGS_t[dim, j] = compute_var_est(res_TGS_t, method="TGS", corr_lag = lag_corr)
    #var_w_tgs_t[dim, j] = var(res_TGS_t$weights)
  }
}
