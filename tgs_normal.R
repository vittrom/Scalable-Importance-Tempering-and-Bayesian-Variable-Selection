#TGS implementation

TGS <- function(start_x, d, T, burn_in, tempering, cond_distr, single_cond, extra_pars, version="mixed"){
  mu <- extra_pars$mu
  sigma <-  extra_pars$sigma
  output_x <- matrix(NA, nrow <- T, ncol <- d)
  sample_weights <- rep(NA, T)
  x <- start_x
  p_x = cond_distr(x, mu, sigma, tempering, version)
  for(iter in 1:(burn_in + T)){
    i <- sample.int(d, 1, prob = p_x)
    #sample from g(x_i|x_{-i})
    x[i] = single_cond(1, x, i, mu, sigma, tempering, version)
    #compute weights
    p_x <- cond_distr(x, mu, sigma, tempering, version)
    if(iter > burn_in){
      output_x[iter - burn_in,] <- x
      sample_weights[iter - burn_in] <- 1/mean(p_x)
    }
  }
  ESS <- sum(sample_weights)^2/sum(sample_weights^2)
  return(list(x=output_x, weights=sample_weights, var_w=var(sample_weights), ESS=ESS))
}

GS <- function(start_x, d, T, burn_in, single_cond, extra_pars, type="deterministic"){
  mu <- extra_pars$mu
  sigma <- extra_pars$sigma
  output_x <- matrix(NA, nrow = T, ncol = d)
  x <- start_x
  for(iter in 1:(burn_in + T)){
    if(type=="deterministic"){
      for(i in 1:d){
        x[i] = single_cond(1, x, i, mu, sigma, 1)
      }
    }
    if(type=="random"){
      i <- sample.int(d, 1)
      x[i] <- single_cond(1, x, i, mu, sigma, 1)
    }
    if(iter > burn_in){
      output_x[iter - burn_in,] <- x
    }
  }
  return(list(x=output_x))
}

cond_moments_normal <- function(x, mu_full, sigma_full, dim){
  block_sigma_dim_dim <- sigma_full[dim, dim]
  block_simga_dim_rest <- sigma_full[dim, -dim]
  block_sigma_rest_dim <- sigma_full[-dim, dim]
  block_sigma_rest_rest <- sigma_full[-dim, -dim]
  inv_block_sigma_rest_rest <- solve(block_sigma_rest_rest)
  mu_cond <- mu_full[dim] + block_simga_dim_rest %*% inv_block_sigma_rest_rest %*% (x[-dim] - mu_full[-dim])
  sigma_sq_cond <- block_sigma_dim_dim - block_simga_dim_rest %*% inv_block_sigma_rest_rest %*% block_sigma_rest_dim
  
  return(list(mu_cond= c(mu_cond), sigma_sq_cond=c(sigma_sq_cond)))
}

sample_g_cond_normal <- function(n, x, dim, mu, sigma, tempering, version="mixed"){
  moments <- cond_moments_normal(x, mu, sigma, dim)
  x_samp <- rnorm(n, mean = moments$mu_cond, sd = sqrt(moments$sigma_sq_cond))
  x_samp_tempered <- rnorm(n, mean = moments$mu_cond, sd = sqrt(moments$sigma_sq_cond/tempering))
  if(version == "vanilla"){return(x_samp_tempered)}
  if(version == "mixed"){
    if(runif(1) < 0.5){return(x_samp)} else{return(x_samp_tempered)}
    # x_samp <- 0.5 * x_samp + 0.5 * x_samp_tempered
    # return(x_samp)
  } else {return(NA)}
}

cond_distr_full_normal <- function(x, mu_full, sigma_full, tempering, version="mixed"){
  dims = length(mu_full)
  p_ix = rep(NA, dims)
  for(i in 1:dims){
    moments = cond_moments_normal(x, mu_full, sigma_full, i)
    vanilla_pix = dnorm(x[i], moments$mu_cond, sqrt(moments$sigma_sq_cond/tempering)) / dnorm(x[i],moments$mu_cond, sqrt(moments$sigma_sq_cond))
    if(is.infinite(vanilla_pix)){
      vanilla_pix = 999999
    } 
    if(version=="vanilla"){
      p_ix[i] = vanilla_pix
    }
    if(version=="mixed"){
      p_ix[i] = 0.5 + 0.5*vanilla_pix
    }
  }
  return(p_ix)
}

plot_results <- function(results, extra_pars, title, limits = c(-3, 3)){
  require(ggplot2)
  require(MASS)
  x <- results$x
  which_in_range <- which(x[,1] >= limits[1] & x[,2] <= limits[2])
  x <- x[which_in_range,]
  weights = NULL
  if(!is.null(results$weights)){
    weights = results$weights[which_in_range]
  }
  
  x_distr = mvrnorm(10000, mu=extra_pars$mu, Sigma = extra_pars$sigma)
  p <- ggplot(as.data.frame(x_distr), aes(x=x_distr[,1], y=x_distr[,2])) +
    xlim(limits) + 
    ylim(limits) +
    stat_density_2d(geom="polygon", aes(fill=stat(level))) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    labs(x="", y="",  title=title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=(15)),legend.position = "none") +
    geom_point(data=as.data.frame(x), aes(x=x[,1], y=x[,2], size=weights))
  if(!is.null(results$ESS)){
    p = p + annotate("label",x=limits[1] + 1, y=limits[2] - 0.5, vjust =1, hjust=1 , label = paste("ESS:" ,round(results$ESS,2), sep = " "))
  }
  p
}

plot_distr <- function(extra_pars, cond_distr, single_cond, tempering, version, limits=c(-3,3)){
  require(MASS)
  require(ggplot2)
  x = mvrnorm(10000, mu=extra_pars$mu, Sigma=extra_pars$sigma)
  d = length(extra_pars$mu)
  
  p_0 = ggplot(as.data.frame(x), aes(x=x[,1], y=x[,2])) +
    ylim(limits) + 
    xlim(limits) + 
    stat_density_2d(geom="polygon", aes(fill=stat(level))) +
    scale_fill_distiller(palette = "YlOrRd", direction=1) + 
    theme_classic() +
    annotate("label", x=limits[1] + 1, y=limits[2] - 0.5, vjust=1, hjust=1, label="f(x)") +
    theme(legend.position = "none") +
    labs(x="", y="")
  
  TGS_vers = TGS(start_x=c(0,0), d=d, T=10000, burn_in=0, tempering=tempering, cond_distr=cond_distr, 
                 single_cond = single_cond, extra_pars=extra_pars, version=version)
  x_1 = TGS_vers$x
  p_1 = ggplot(as.data.frame(x_1), aes(x=x_1[,1], y=x_1[,2])) +
    ylim(limits) + 
    xlim(limits) + 
    stat_density_2d(geom="polygon", aes(fill=stat(level))) +
    scale_fill_distiller(palette = "YlOrRd", direction=1) + 
    theme_classic() +
    annotate("label", x=limits[1] + 1, y=limits[2] - 0.5, vjust=1, hjust=1, label="f(x)Z(x)") +
    theme(legend.position = "none") +
    labs(x="", y="")
  
  return(list(f_x = p_0, fz = p_1))
}

#TESTS

rho = 0.9
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
GS_res_det <- GS(start_x, d, T, burn_in, sample_g_cond_normal, extra_pars, type = "deterministic")
GS_res_rs <- GS(start_x, d, T, burn_in, sample_g_cond_normal, extra_pars, type = "random")

TGS_res_mixed$ESS
TGS_res_vanilla$ESS

plot_results(TGS_res_vanilla, extra_pars, "TGS_vanilla")
