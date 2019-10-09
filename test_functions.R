rm(list=ls())
setwd("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code")

source("tgs_normal.R")

rho_list = c(0, 0.5, .9, 0.999)
plots = rep(NA, length(rho_list))
counter = 1
for(rho in rho_list){
  extra_pars = list(mu=c(0,0), sigma=matrix(c(1,rho,rho,1), nrow = 2))
  tempering= 1-rho^2
  ps = plot_distr(extra_pars=extra_pars, cond_distr = cond_distr_full_normal,
                  single_cond = sample_g_cond_normal, tempering=tempering, version="vanilla")
  plots[counter] = ps[1]
  counter = counter + 1
  plots[counter] = ps[2]
  counter = counter + 1
}

grid.arrange(plots[[2]], plots[[4]], plots[[6]], plots[[8]], ncol=4)
p=1
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

p=2
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

p=3
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

p=4
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

p=5
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

p=6
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

p=7
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

p=8
plots[[p]]
ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")

for(p in  1:length(plots)){
  plots[[p]]
  
  ggsave(paste("plot_", p, ".pdf", sep=""), path = "./Plots/", width = 10, height = 6.5, units = "cm")
}

rho_list = c(0, 0.25, 0.5, .75, 0.999)
tempering = 1- rho_list^2
rho = 0.999
start_x = c(3,3);d=2; T=200;burn_in=0; update=T/5
plots = rep(NA, length(rho_list))
counter = 1
extra_pars = list(mu=c(0,0), sigma=matrix(c(1,rho,rho,1), nrow = 2))

for(temper in tempering){
  TGS_res_mixed = TGS(start_x, d, T, burn_in, temper, cond_distr = cond_distr_full_normal,
                      single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")
  plot_results(TGS_res_mixed, extra_pars, paste("\u03b2 = ", round(temper, 5)), width = 10, save=TRUE)
  counter = counter + 1
}

#TEST Functions
rho = 0.999
sigma=diag(10)
sigma[sigma!= 1] = rho#(c(1,rho,rho,1)
extra_pars = list(mu=rep(0,10), sigma=sigma)
d=2;start_x = rep(3, d); T=10000;burn_in=0; update=T/5; reps=1000
ESS = matrix(NA,  reps,  4)
est = matrix(NA,  reps,  4)
#Compare TGS_adaptive vs TGS
for(i in 1:reps){
  tempering = runif(1)
  
  TGS_res_mixed = TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
                      single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")
  TGS_ada_mixed = TGS_adaptive(start_x, d, T, burn_in, cond_distr = cond_distr_full_normal,update=update,
                               single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed", tempering_start = tempering)
  
  TGS_ada_annealing = TGS_adaptive(start_x, d, T, burn_in, cond_distr = cond_distr_full_normal,update=update, annealing=TRUE,tempering_adj=FALSE,
                                   single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed", tempering_start = tempering)
  
  TGS_ada_ann_temp = TGS_adaptive(start_x, d, T, burn_in, cond_distr = cond_distr_full_normal,update=update, annealing=TRUE, tempering_adj=TRUE,
                                  single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed", tempering_start = tempering)
  ESS[i, 1] = T/(1 + 2*sum(acf(TGS_res_mixed$x[,1]*TGS_res_mixed$weights, plot = FALSE)$acf))
  ESS[i, 2] = T/(1 + 2*sum(acf(TGS_ada_mixed$x[,1]*TGS_ada_mixed$weights, plot = FALSE)$acf))
  ESS[i, 3] = T/(1 + 2*sum(acf(TGS_ada_annealing$x[,1]*TGS_ada_annealing$weights, plot = FALSE)$acf))
  ESS[i, 4] = T/(1 + 2*sum(acf(TGS_ada_ann_temp$x[,1]*TGS_ada_ann_temp$weights, plot = FALSE)$acf))
  
  est[i, 1] = sum(TGS_res_mixed$x[,1]*TGS_res_mixed$weights)/sum(TGS_res_mixed$weights)
  est[i, 2] = sum(TGS_ada_mixed$x[,1]*TGS_ada_mixed$weights)/sum(TGS_ada_mixed$weights)
  est[i, 3] = sum(TGS_ada_annealing$x[,1]*TGS_ada_annealing$weights)/sum(TGS_ada_annealing$weights)
  est[i, 4] = sum(TGS_ada_ann_temp$x[,1]*TGS_ada_ann_temp$weights)/sum(TGS_ada_ann_temp$weights)
}

apply(ESS, 2, mean)
apply(ESS, 2, median)
apply(ESS, 2, max)
apply(ESS, 2, min)
apply(ESS, 2, sd)

apply(est, 2, mean)
apply(est, 2, median)
apply(est, 2, max)
apply(est, 2, min)
apply(est, 2, sd)


sum(TGS_ada_ann_temp$x[,1] * TGS_ada_ann_temp$weights)/sum(TGS_ada_ann_temp$weights)
sum(TGS_ada_annealing$x[,1] * TGS_ada_annealing$weights)/sum(TGS_ada_annealing$weights)
sum(TGS_ada_mixed$x[,1] * TGS_ada_mixed$weights)/sum(TGS_ada_mixed$weights)
sum(TGS_res_mixed$x[,1] * TGS_res_mixed$weights)/sum(TGS_res_mixed$weights)

plot_results(TGS_res_mixed, extra_pars, "TGS")
plot_results(TGS_ada_mixed, extra_pars, "TGS feedback", save = TRUE)
plot_results(TGS_ada_annealing, extra_pars, "TGS annealing", save=TRUE)
plot_results(TGS_ada_ann_temp, extra_pars, "TGS annealing + feedback", save=TRUE)

par(mfrow=c(1,4))
acf(TGS_ada_mixed$x[,1] * TGS_ada_mixed$weights, 200)
acf(TGS_res_mixed$x[,1] * TGS_res_mixed$weights, 200)
acf(TGS_ada_annealing$x[,1] * TGS_ada_mixed$weights, 200)
acf(TGS_ada_ann_temp$x[,1] * TGS_ada_ann_temp$weights, 200)


rho = 0.999
extra_pars = list(mu=rep(0,2), sigma=matrix(c(1,rho, rho, 1), nrow=2))
d=2;start_x = rep(3, d); T=200;burn_in=0; tempering=1-rho^2

TGS_res_mixed = TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
                    single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")
TGS_res_vanilla = TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
                    single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="vanilla")
GS_res = GS(start_x, d, T, burn_in, cond_moments_normal, extra_pars)

plot_results(TGS_res_vanilla, extra_pars, "TGS vanilla", save=TRUE)
plot_results(TGS_res_mixed, extra_pars, "TGS mixed", save=TRUE)
plot_results(GS_res, extra_pars, "GS", save=TRUE)

TGS_ada_mixed = TGS_adaptive(start_x, d, T, burn_in, cond_distr = cond_distr_full_normal,update=update,
                             single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed", tempering_start = tempering)

TGS_ada_annealing = TGS_adaptive(start_x, d, T, burn_in, cond_distr = cond_distr_full_normal,update=update, annealing=TRUE,tempering_adj=FALSE,
                                 single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed", tempering_start = tempering)

TGS_ada_ann_temp = TGS_adaptive(start_x, d, T, burn_in, cond_distr = cond_distr_full_normal,update=update, annealing=TRUE, tempering_adj=TRUE,
                                single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed", tempering_start = tempering)
par(mfrow=c(1,4))
plot_weights(TGS_ada_ann_temp, "Annealing + Feedback", save=TRUE)
plot_weights(TGS_ada_mixed, "Feedback", save = TRUE)
plot_weights(TGS_res_mixed, "TGS", save=TRUE)
plot_weights(TGS_ada_annealing, "Annealing", save=TRUE)