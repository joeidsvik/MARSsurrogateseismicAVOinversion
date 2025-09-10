
match_traveltimes_mean_trend <- function(tt){
  # Estimated conversion constant for timetravel to depth conversion
  convrt <- 1.055
  # Loading the depth trends
  prior_mean <- matrix(NA, nrow = prod(dim(tt)), ncol = 4)
  prior_mean[,4] <- matrix(tt*convrt, nrow = prod(dim(tt)), ncol = 1, byrow = F)
  #browser()
  for(i in 1:prod(dim(tt))){
    if(length(which(floor(depth_trend_prior[,1])==floor(prior_mean[i,4])))==0){
      d_i <- 1 
    }else{
    d_i <- which(floor(depth_trend_prior[,1])==floor(prior_mean[i,4])) }
      
    prior_mean[i,3]<- depth_trend_prior[d_i,2]
    prior_mean[i,1]<- depth_trend_prior[d_i,3]
    prior_mean[i,2]<- depth_trend_prior[d_i,4]
  }
  
  return(prior_mean)
  
  
}

prior_gen_from_tt <- function(tt, n_real){
  # Input
  # tt - traveltime matrix
  # n_real - number of realisations to be generated
  #### Prior inputs ####
  phi_cl <- 15
  phi_sat <- 15
  # estimated standard deviations for the three variables - constant across the field
  #std_cl <- 0.8
  std_cl <- 1.64
  #std_g <- 2.5
  std_g <- 2.83
  #std_o <- 2.5
  std_o <- 2.68
  # parameter for hyper-prior variance
  ksi <- 0.01
  
  
  dim_tt <- dim(tt)
  
  # match the traveltimes with the prior
  
  convrt <- 1.055
  # Loading the depth trends
  prior_mean <- matrix(NA, nrow = prod(dim(tt)), ncol = 4)
  prior_mean[,4] <- matrix(tt*convrt, nrow = prod(dim(tt)), ncol = 1, byrow = F)
  #browser()
  for(i in 1:prod(dim(tt))){
    d_i <- which(floor(depth_trend_prior[,1])==floor(prior_mean[i,4])) 
    # oil
    prior_mean[i,3]<- depth_trend_prior[d_i,2]
    # shale
    prior_mean[i,1]<- depth_trend_prior[d_i,3]
    # gas
    prior_mean[i,2]<- depth_trend_prior[d_i,4]
  }
  
  mu_cl <-  matrix(prior_mean[,1], nrow = dim_tt[2], ncol = dim_tt[1])
  mu_o <-  matrix( prior_mean[,3], nrow = dim_tt[2], ncol = dim_tt[1])
  mu_g <- matrix( prior_mean[,2], nrow = dim_tt[2], ncol = dim_tt[1])
  
  x_rf_clay <- matrix(NA, nrow = prod(dim(tt)), ncol = n_real)
  x_rf_sat_g <- matrix(NA, nrow = prod(dim(tt)), ncol = n_real)
  x_rf_sat_o <- matrix(NA, nrow = prod(dim(tt)), ncol = n_real)
  
  # Generating both clay and saturation prior unconditioned 
  #browser()
  # The values are placed in the vector column-wise
  for(i in 1:n_real){
    x_rf_clay[,i] <- fft2_rf_realis(dim_tt[1],dim_tt[2],std_cl,ksi,phi_cl) + mu_cl 
    x_rf_sat_g[,i] <- fft2_rf_realis(dim_tt[1],dim_tt[2],std_g,ksi,phi_sat) + mu_g
    x_rf_sat_o[,i] <- fft2_rf_realis(dim_tt[1],dim_tt[2],std_o,ksi,phi_sat) + mu_o
    
    
  }
  
  prior_real <- list(x_g = x_rf_sat_g, x_o = x_rf_sat_o,x_cl = x_rf_clay)
  #names(prior_real) <- 
  
  
  return(prior_real)
}









sim_forward_model <- function(x_rf_sat_g, x_rf_sat_o, x_rf_clay, traveltimes){
  # Specifying total dimension from traveltimes input
  
  if(is.null(dim(traveltimes))==TRUE){
    if(length(traveltimes)>=1){
      tot_dim <- length(traveltimes)
    }else{
    stop("Illegal size of traveltimes or something else is wrong.")
    #dim_grid <- c(1,1)
    }
  
  }else{
    tot_dim <- prod(dim(traveltimes))
    
    # Specifying [n m] dimension of the grid 
    #dim_grid <- dim(traveltimes)
}
  
  ## settting up everything needed for the data part 
  # Constant variables
  phi0_ss <- 0.4
  phi0_sh <- 0.6
  
  
  # Estimated conversion constant for timetravel to depth conversion
  convrt <- 1.055
  
  # cementation depth
  depth_cem <- 2000*convrt
  
  ref_depth <- 900*convrt
  
  
  # parameter for sand porosity
  alpha_ss <- 0.32
  
  # parameter for shale porosity
  alpha_sh <- 1.1
  kappa_ss <- 0.18
  
  # reference depth for the empirical prior for porosity
  #ref_depth <- 900*convrt
  # testing with 0 ref depth
  ref_depth <- ref_depth 
  # min and max porosity for the logistic transformation
  
  
  
  #depth included here
  n_param_mean <- 4
  
  
  poros_ss_mean <- matrix(NA, nrow=1, ncol = tot_dim)
  poros_sh_mean <- matrix(NA, nrow=1, ncol = tot_dim)
  depth_grid <- matrix(traveltimes*convrt,nrow=tot_dim, ncol =1, byrow = F)
  
  
  for(i in 1:tot_dim){
    elem <- depth_grid[i]
    
    if(elem=="NaN"){
      poros_ss_mean[1,i] <- NaN
      poros_sh_mean[1,i] <- NaN
    }else{
      poros_ss_mean[1,i] <- porosity_ss(phi0_ss, alpha_ss, elem, ref_depth , depth_cem , kappa_ss)
      poros_sh_mean[1,i] <- porosity_sh(phi0_sh, alpha_sh, elem, ref_depth)
    }
  }
  
 
  x_rf_clay_log <- exp(x_rf_clay)/(1+ exp(x_rf_clay))
  
  # x_rf_clay are flat and have only one realisations 
  mu_poros <- x_rf_clay_log*poros_sh_mean[1,] + (1-x_rf_clay_log)*poros_ss_mean[1,]
  
  
  
  sat_log_tot <- 1 + exp(x_rf_sat_g) + exp(x_rf_sat_o)
  
  x_rf_sat_g_log <- (exp(x_rf_sat_g)/(sat_log_tot))
  x_rf_sat_o_log <- (exp(x_rf_sat_o)/sat_log_tot)
  x_rf_sat_b_log <- (1/sat_log_tot)
  
  n_obs <- 2
  # simulating data part
  
  y_obs_sim_noise_free <- matrix(NA, nrow = (tot_dim*n_obs), ncol = 1)
  elast_param <- matrix(NA, nrow = (tot_dim*3), ncol = 1)
  
  
  
  x_rf_poros_ss_log <- poros_ss_mean[1,]
  x_rf_poros_sh_log <- poros_sh_mean[1,]
  x_rf_poros_log <- mu_poros
  
  
 
  
  temp_forward <- rp_likelihood(depth_cem, tot_dim, x_rf_poros_ss_log, x_rf_poros_sh_log, x_rf_sat_g_log, x_rf_sat_o_log, x_rf_sat_b_log, x_rf_clay_log,traveltimes,x_rf_poros_log) 
  elast_param <- temp_forward[[1]]
  y_obs_sim_noise_free <- temp_forward[[2]]
    

  
  return(y_obs_sim_noise_free)
}