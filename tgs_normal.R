#TGS implementation + additional functions

TGS <- function(start_x, d, T, burn_in, tempering, cond_distr, single_cond, extra_pars, version="mixed", shape=0.2){
  mu <- extra_pars$mu
  sigma <-  extra_pars$sigma
  output_x <- matrix(NA, nrow <- T, ncol <- d)
  sample_weights <- rep(NA, T)
  x <- start_x
  ESS_running = rep(NA, T)
  if(version=="t"){
    p_x = cond_distr(x, mu, sigma, shape)
  }else{
    p_x = cond_distr(x, mu, sigma, tempering, version)
  }
  
  for(iter in 1:(burn_in + T)){
    i <- sample.int(d, 1, prob = p_x)
    #sample from g(x_i|x_{-i})
    if(version=="t"){
      x[i] = single_cond(1, x, i, mu, sigma, shape)
    }else{
      x[i] = single_cond(1, x, i, mu, sigma, tempering, version)
    }
    #compute weights
    if(version=="t"){
      p_x = cond_distr(x, mu, sigma, shape)
    }else{
      p_x = cond_distr(x, mu, sigma, tempering, version)
    }
    if(iter > burn_in){
      output_x[iter - burn_in,] <- x
      sample_weights[iter - burn_in] <- 1/mean(p_x)
    }
    ESS_running[iter] = sum(sample_weights, na.rm = TRUE)^2/sum(sample_weights^2, na.rm = TRUE)
  }
  ESS <- sum(sample_weights)^2/sum(sample_weights^2)
  return(list(x=output_x, weights=sample_weights, var_w=var(sample_weights), ESS=ESS, ESS_running=ESS_running))
}

GS <- function(start_x, d, T, burn_in, cond_moments, extra_pars, type="deterministic"){
  mu <- extra_pars$mu
  sigma <- extra_pars$sigma
  output_x <- matrix(NA, nrow = T, ncol = d)
  x <- start_x
  for(iter in 1:(burn_in + T)){
    if(type=="deterministic"){
      for(i in 1:d){
        moments = cond_moments(x, mu, sigma, i)
        x[i] = rnorm(1, moments$mu, sqrt(moments$sigma_sq_cond))
      }
    }
    if(type=="random"){
      i <- sample.int(d, 1)
      moments = cond_moments(x, mu, sigma, i)
      x[i] = rnorm(1, moments$mu, sqrt(moments$sigma_sq_cond))
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
  } else {return(NA)}
}

sample_g_cond_t <- function(n, x, dim, mu, sigma, shape=0.2){
  moments = cond_moments_normal(x, mu, sigma, dim)
  mu_t = moments$mu_cond
  s_t = sqrt(moments$sigma_sq_cond)
  samp = rt(n, shape) * s_t + mu_t
  return(samp)
}

cond_distr_full_t <- function(x, mu_full, sigma_full, shape=0.2){
  dims = length(mu_full)
  p_ix = rep(NA, dims)
  for(i in 1:dims){
    moments = cond_moments_normal(x, mu_full, sigma_full, i)
    mu_t = moments$mu_cond
    s_t = sqrt(moments$sigma_sq_cond)
    pix = 1/s_t * dt((x[i] - mu_t)/s_t, shape) / dnorm(x[i], mu_t, s_t)
    if(is.infinite(pix)){
      pix = 999999
    } 
    p_ix[i] = pix
  }
  return(p_ix)
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

plot_results <- function(results, extra_pars, title, limits = c(-3, 3), width=8, height=6.5, save=FALSE){
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
    p = p + annotate("label",x=limits[1] + 2.5, y=limits[2] - 0.5, vjust =1, hjust=1 , label = paste("ESS:" ,round(results$ESS,2), sep = " "))
  }
  p
  if(save){
    ggsave(paste(title,".pdf", sep=""), width = width, height = width, units = "cm",path = "./Plots")
  }
  p
}

plot_distr <- function(extra_pars, cond_distr, single_cond, tempering, version, algo="TGS", limits=c(-3,3)){
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
    annotate("label", x=limits[1] + 1, y=limits[2] - 0.5, vjust=1, hjust=1, label=expression(pi~"(x)")) +
    theme(legend.position = "none") +
    labs(x="", y="")
  
  
  if(algo == "ada"){
    TGS_vers = TGS_adaptive(start_x=c(0,0), d=d, T=10000, burn_in=0, tempering_start=tempering, cond_distr=cond_distr, 
                   single_cond = single_cond, extra_pars=extra_pars, version=version)
  } else{
    TGS_vers = TGS(start_x=c(0,0), d=d, T=10000, burn_in=0, tempering=tempering, cond_distr=cond_distr, 
                   single_cond = single_cond, extra_pars=extra_pars, version=version)
  }
  
  x_1 = TGS_vers$x
  p_1 = ggplot(as.data.frame(x_1), aes(x=x_1[,1], y=x_1[,2])) +
    ylim(limits) + 
    xlim(limits) + 
    stat_density_2d(geom="polygon", aes(fill=stat(level))) +
    scale_fill_distiller(palette = "YlOrRd", direction=1) + 
    theme_classic() +
    annotate("label", x=limits[1] +2, y=limits[2] - 0.5, vjust=1, hjust=1, label=expression(pi~"(x)Z(x)")) +
    theme(legend.position = "none") +
    labs(x="", y="")
  
  return(list(f_x = p_0, fz = p_1))
}

plot_weights <- function(results,  title, width=10, height=6.5, save=FALSE){
  weights = as.data.frame(results$weights)
  n = nrow(weights)
  weights = as.data.frame(cbind( 1:n, weights))
  p = ggplot(weights, aes(x=weights[,1], y = weights[,2])) + 
    geom_point(alpha=0.1) +
    labs(x="Iteration", y="Weight", title = title) +
    theme(plot.title = element_text(hjust = 0.5))
  p
  if(save){
    ggsave(paste(title, ".pdf", sep=""), width = width, height = height, units = "cm", path = "./Plots/")
  }
  p
}

TGS_adaptive <- function(start_x, d, T, burn_in, cond_distr, single_cond, extra_pars, 
                         tempering_start, update =50, version="mixed", shape=0.2, eps=0.05, annealing=FALSE, tempering_adj = TRUE){
  mu <- extra_pars$mu
  sigma <-  extra_pars$sigma
  output_x <- matrix(NA, nrow <- T, ncol <- d)
  sample_weights <- rep(NA, T)
  x <- start_x
  lags = ceiling(10*log10(update))
  tempering = tempering_start

  tempering_seq = rep(NA, ceiling(T/update))
  tempering_seq[1] = tempering
  
  if(version=="t"){
    p_x = cond_distr(x, mu, sigma, shape)
  }else{
    p_x = cond_distr(x, mu, sigma, tempering, version)
  }
  
  tmp = lags
  for(iter in 1:(burn_in + T)){
    i <- sample.int(d, 1, prob = p_x)
    #sample from g(x_i|x_{-i})
    if(version=="t"){
      x[i] = single_cond(1, x, i, mu, sigma, shape)
    }else{
      x[i] = single_cond(1, x, i, mu, sigma, tempering, version)
    }
    
    #compute weights
    if(version=="t"){
      p_x = cond_distr(x, mu, sigma, shape)
    }else{
      p_x = cond_distr(x, mu, sigma, tempering, version)
    }
  
    if(iter > burn_in){
      output_x[iter - burn_in,] <- x
      sample_weights[iter - burn_in] <- 1/mean(p_x)
      if(annealing){
        tempering = tempering_start * (1 - (iter - burn_in)/(T+1))
      }
    }
    
    if(tempering_adj){
      if(iter - burn_in == update){
        samps = output_x[(iter- burn_in - update + 1): (iter-burn_in), i] * sample_weights[(iter - burn_in -update + 1): (iter - burn_in)]
        corr = acf(samps, plot=FALSE)

        tmp = sum(corr$acf)/var(samps)
        tempering_seq[round((iter-burn_in)/update, 0) + 1] = tempering

    }
      if((iter-burn_in) %% update == 0 & iter > (burn_in + update)){
          # max_weight = max(sample_weights, na.rm = TRUE)
          samps = output_x[(iter- burn_in - update + 1): (iter-burn_in), i] * sample_weights[(iter - burn_in - update + 1): (iter - burn_in)]
          corr = acf(samps, plot=FALSE)

          #novar implementation
          corr_sum = sum(corr$acf)/var(samps)
          ratio = corr_sum/tmp

          if(ratio > 1){
            tempering = runif(1, 1e-6, tempering)
          } else {
            tempering = min(max(runif(1,tempering*(1 - eps), tempering * (1 + eps)), 1e-6),1)
          }
          tempering_start = tempering
          tempering_seq[round((iter-burn_in)/update, 0) + 1] = tempering
          tmp = corr_sum
      }
    }
  }
  ESS <- sum(sample_weights)^2/sum(sample_weights^2)
  return(list(x=output_x, weights=sample_weights, var_w=var(sample_weights), ESS=ESS, tempering_seq=tempering_seq, 
              update=update))
}

