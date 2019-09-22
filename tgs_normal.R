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
  return(output_x)
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

rho = 0.999
extra_pars = list(mu=c(0,0), sigma=matrix(c(1,rho,rho,1), nrow = 2))
start_x = c(3,3);d=2; T=200;burn_in=0; tempering=1-rho^2
TGS_res_mixed <- TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal,
               single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="mixed")

TGS_res_vanilla <- TGS(start_x, d, T, burn_in, tempering, cond_distr = cond_distr_full_normal, 
                     single_cond = sample_g_cond_normal, extra_pars = extra_pars, version="vanilla")
GS_res_det <- GS(start_x, d, T, burn_in, sample_g_cond_normal, extra_pars, type = "deterministic")
GS_res_rs <- GS(start_x, d, T, burn_in, sample_g_cond_normal, extra_pars, type = "random")

TGS_res_mixed$ESS
TGS_res_vanilla$ESS
# beta=0.6
# x_norm_tempered_v2 = rnorm(10000, 0, 1/sqrt(beta)) * (1/sqrt(beta)) * (sqrt(2*pi))*(1/(sqrt(2*pi)))^beta
# x_norm = rnorm(1000, 0, 1)
# x_norm_tempered = sign(x_norm) * abs(x_norm)^beta
# hist(x_norm_tempered)
# hist(x_norm_tempered_v2)
