rm(list=ls())
source("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code\\tgs_normal.R")

#TEST Functions
rho = 0.999
extra_pars = list(mu=c(0,0), sigma=matrix(c(1,rho,rho,1), nrow = 2))
start_x = c(3,3);d=2; T=200;burn_in=0; 



#Compare TGS_adaptive vs TGS
tempering = runif(1)

TGS_res_mixed = TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
                    single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")
                    
TGS_ada_mixed = TGS_adaptive(start_x, d, T, burn_in, cond_distr = cond_distr_full_normal,update=75,
             single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed", tempering_start = tempering)

sum(TGS_ada_mixed$x[,1] * TGS_ada_mixed$weights)/sum(TGS_ada_mixed$weights)
sum(TGS_res_mixed$x[,1] * TGS_res_mixed$weights)/sum(TGS_res_mixed$weights)

# GS_res = GS(start_x, d, T, burn_in, cond_moments_normal, extra_pars)

# plot_results(GS_res, extra_pars, "GS")
# 
plot_results(TGS_res_mixed, extra_pars, "TGS")
plot_results(TGS_ada_mixed, extra_pars, "TGS ada")

par(mfrow=c(1,2))
acf(TGS_ada_mixed$x[,1] * TGS_ada_mixed$weights, 200)
acf(TGS_res_mixed$x[,1] * TGS_res_mixed$weights, 200)

# par(mfrow=c(1,2))
# plot(TGS_ada_mixed$weights)
# plot(TGS_res_mixed$weights)
# 
# plot_distr(extra_pars = extra_pars, cond_distr = cond_distr_full_normal, single_cond = sample_g_cond_normal,
#            tempering= runif(1), version = "mixed", algo = "TGS" )

# plot_results(TGS_ada_mixed, extra_pars, "GS random scan")
# TGS_ada_mixed$ESS
# v_i = TGS_ada_mixed$x[,1] * TGS_ada_mixed$weights
# start = seq(1, length(v), 50)
# end = seq(50, length(v), 50)
# res = rep(NA, length(start))
# tmp = 1
# frac = rep(NA, length(start))
# for(i in 1:length(start)){
#   sub = v[start[i]: end[i]]
#   a = acf(sub, plot=FALSE)
#   res[i] = sum(abs(a$acf))/var(sub)
#   frac[i] = res[i]/tmp
#   tmp = res[i]
# }
# mean(res)


# tempering = runif(1)
# TGS_res_mixed <- TGS(start_x, d, T, burn_in,tempering, cond_distr = cond_distr_full_normal,
#                      single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")
# 
# v = TGS_res_mixed$x[,1] * TGS_res_mixed$weights
# 
# start = seq(1, length(v), 50)
# end = seq(50, length(v), 50)
# res = rep(NA, length(start))
# tmp = 17
# frac = rep(NA, length(start))
# for(i in 1:length(start)){
#   sub = v[start[i]: end[i]]
#   a = acf(sub, plot=FALSE)
#   res[i] = sum(abs(a$acf))/var(sub)
#   frac[i] = res[i]/tmp
#   tmp=res[i]
# }
# frac
# # plot(res, type="l")
# res_mean[j] = mean(res)




#acf(TGS_res_mixed$x[,1] * TGS_res_mixed$weights, 200)
#TGS_res_mixed$ESS
#plot_results(TGS_res_mixed, extra_pars, "GS random scan")



# TGS_res_vanilla <- TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
#                      single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="vanilla")

# GS_res_det <- GS(start_x, d, T, burn_in, cond_moments_normal, extra_pars, type = "deterministic")
# plot_results(GS_res_det, extra_pars, "GS")
# 
# GS_res_rs <- GS(start_x, d, T, burn_in, cond_moments_normal, extra_pars, type = "random")
# plot_results(GS_res_rs, extra_pars, "GS random scan")

# TGS_res_vanilla$ESS
# plot_results(TGS_res_mixed, extra_pars, "GS random scan")
# 
# plot(TGS_res_mixed$ESS_running, type="l")


#tempering= (1-rho^2)/2

# ps = plot_distr(extra_pars=extra_pars, cond_distr = cond_distr_full_normal,
#                 single_cond = sample_g_cond_normal, tempering=tempering, version="vanilla")
# ps[1]
# ps[2]

# tempering= seq(0.0001, 1, 0.002) #runif(1)#(1-rho^2)
# res_mean = rep(NA, length(tempering))
# for(j in 1:length(tempering)){
#   TGS_res_mixed <- TGS(start_x, d, T, burn_in, tempering[j], cond_distr = cond_distr_full_normal,
#                        single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")
#   
#   v = TGS_res_mixed$x[,1] * TGS_res_mixed$weights
#   
#   start = seq(1, length(v), 50)
#   end = seq(50, length(v), 50)
#   res = rep(NA, length(start))
#   tmp = 1
#   for(i in 1:length(start)){
#     sub = v[start[i]: end[i]]
#     a = acf(sub, plot=FALSE)
#     res[i] = sum(abs(a$acf))/var(sub)
#   }
#   # plot(res, type="l")
#   res_mean[j] = mean(res)
# }

