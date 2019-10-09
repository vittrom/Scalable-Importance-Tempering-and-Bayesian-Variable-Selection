rm(list= ls())
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

run_var_est <- function(n, extra_pars, tempering, cond_distr_full){
  require(MASS)
  mu = extra_pars$mu
  sigma = extra_pars$sigma
  x = mvrnorm(n, mu, sigma)
  p_x = apply(x, 1, cond_distr_full, mu, sigma, tempering)
  weights = 1/apply(p_x, 1, mean)
  
  return(mean(x[,1]^2 * weights))
}

compute_var_est <- function(x, extra_pars=NULL, tempering=NULL, cond_distr_full=NULL, weights=NULL, method="NULL", corr_lag=200){
  
  h = x
  corr = acf(h, plot = FALSE, type = "correlation")
  var_h = 1
  
  if(method=="TGS"){
    f_x =  (h * weights)
    corr = acf(f_x, plot=FALSE, type="correlation")
    var_h = run_var_est(n=10000, extra_pars, tempering, cond_distr_full)
  }
  
  var_est = var_h * (1 + 2 * sum(corr$acf))
  return(var_est)
}

plot_correlation = function(var_GS, var_TGS, ts, save=FALSE, width=8.5, height=6){
  n_plots = ncol(var_GS)
  for(i  in 1:n_plots){
    v_GS = as.data.frame(cbind(ts,var_GS[,i]))
    v_TGS = as.data.frame(cbind(ts,var_TGS[,i]))
    
    p = ggplot(v_GS, aes(x=v_GS[,1], y=v_GS[,2])) +
      geom_line(aes(colour="GS")) +
      geom_line(data=v_TGS, aes(x=v_TGS[,1], y=v_TGS[,2], colour="TGS")) + 
      labs(x= expression("1/(1 -"~rho~")"), y="Asymptotic variance", title = scenario[i], color="Method") +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +  
      scale_y_log10() + scale_x_log10()
    if(i < n_plots){
      p = p + theme(legend.position = "none")
    }
    p
    if(save){
      if(i < n_plots){
        ggsave(paste(scenario[i], ".pdf", sep=""), path="./Plots/", width=width-1, height=height, units="cm")
      } else{
        ggsave(paste(scenario[i], ".pdf", sep=""), path="./Plots/", width=width, height=height, units="cm")
      }
    }
    p
  }
}

#Scenario 1: pairwise correlation
#Scenario 2: all positive correlation
#Scenario 3: all negative correlation

d=10; burn_in=0#5000
rho = c(0, 0.9, 0.99,   0.999, 0.9999)
ts = 1/(1-rho)
reps = length(rho); n_scenarios = 3

tempering = 1-rho^2
var_TGS = matrix(NA, nrow = reps, ncol = n_scenarios)
var_GS = matrix(NA, nrow = reps, ncol = n_scenarios)
var_TGS_t =  matrix(NA, nrow = reps, ncol = n_scenarios)
lag_corr = 100; shape=0.2
iter = 5000
start_x = rep(3, d)

var_w_tgs = matrix(NA, nrow=reps, ncol=n_scenarios)
var_w_tgs_t = matrix(NA, nrow=reps, ncol=n_scenarios)

for(j in 1:n_scenarios){
  for(dim in 1:reps){
    extra_pars = create_moments(scenario=j, rho[dim], d)
    
    res_TGS = TGS(start_x = start_x, d=d, T=iter, burn_in=burn_in, tempering=tempering[dim], cond_distr = cond_distr_full_normal,
              single_cond = sample_g_cond_normal, extra_pars = extra_pars, version = "mixed")
    
    var_TGS[dim, j] = compute_var_est(res_TGS$x[,1], extra_pars=extra_pars, method="TGS", weights = res_TGS$weights,
                                      tempering = tempering[dim], cond_distr_full = cond_distr_full_normal)

    res_GS = GS(start_x = start_x, d=d, T=iter, burn_in=burn_in, cond_moments = cond_moments_normal,
                  extra_pars = extra_pars, type = "deterministic")
    var_GS[dim, j] = compute_var_est(res_GS$x[,1])
     
    # res_TGS_t = TGS(start_x = start_x, d=d, T=iter, burn_in=burn_in, tempering=tempering[dim], cond_distr = cond_distr_full_t,
    #                 single_cond = sample_g_cond_t, extra_pars = extra_pars, version = "t", shape=shape)
    # var_TGS_t[dim, j] = compute_var_est(res_TGS_t$x[,1], weights = res_TGS_t$weights, method="TGS", corr_lag = lag_corr,
    #                                     extra_pars = extra_pars, cond_distr_full = cond_distr_full_normal)
  }
}

scenario = c("Pairwise correlation", "Block positive correlation", "Block negative correlation")

plot_correlation(var_GS, var_TGS, ts, save=TRUE)
