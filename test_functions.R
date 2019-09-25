source("D:\\Studies\\PhD\\Year1\\STAT548 - Qualifying course\\Paper 1 - TGS\\Code\\tgs_normal.R")

#TEST Functions
rho = 0.999
extra_pars = list(mu=c(0,0), sigma=matrix(c(1,rho,rho,1), nrow = 2))
start_x = c(3,3);d=2; T=200;burn_in=0; tempering=1-rho^2

ps = plot_distr(extra_pars=extra_pars, cond_distr = cond_distr_full_normal,
                single_cond = sample_g_cond_normal, tempering=tempering, version="vanilla")
ps[1]
ps[2]

TGS_res_mixed <- TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
               single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")

TGS_res_vanilla <- TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
                     single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="vanilla")

GS_res_det <- GS(start_x, d, T, burn_in, cond_moments_normal, extra_pars, type = "deterministic")
plot_results(GS_res_det, extra_pars, "GS")

GS_res_rs <- GS(start_x, d, T, burn_in, cond_moments_normal, extra_pars, type = "random")
plot_results(GS_res_rs, extra_pars, "GS random scan")

TGS_res_mixed$ESS
TGS_res_vanilla$ESS