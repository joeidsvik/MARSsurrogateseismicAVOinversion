# Likelihood model function

# Full forward model with unconsolidated sand, no slip factor, contact cement model, (not fully implemented) constant-clay model 
# and (not added yet) saturation model

source("rock_physics.R")
source("functions_synthetic_case.R")

# Heimdal data included through the well_logs file
#source("well_logs.R")

#library(fields)






rp_likelihood <- function(depth_cem, tot_dim, poros_ss, poros_sh, 
                          sat_gas, sat_oil, sat_brine, clay_content, traveltimes, 
                          poros_tot){
  #### INPUT SPECIFICATIONS #####
  ## Dimension expected of dim_grid, poros_ss, poros_sh, sat_gas, sat_oil, sat_brine, clay_content
  
  
  #################
  ####CONSTANTS####
  #################
  
  convrt <- 1.055
  
  
  phi_ssc <- 0.4
  phi_shc <- 0.6
  phi_ssc <- 0.4
  
  #phi_ssc <- 1.5
  
  
  k_quartz <- 37*10^9
  mu_quartz <- 44*10^9
  
  # Alternative could be: 28.45, average of illite and smectite 
  #k_clay <- 28.45*10^9
  # Alternative could be: 9.6, average of illite and smectite
  #mu_clay <- 9.6*10^9
  k_clay <- 25*10^9
  mu_clay <- 6*10^9
  rho_clay <- 2.65*1000
  
  
  #depth_dt <- (1700:2250)*convrt
  
  # bulk modulus: pa (kg/ms2)
  # water free
  k_brine <- 2.4*10^9
  #k_brine <- 4*10^9
  # oil
  k_fluid_o <- 10^9
  # gass
  k_fluid_g <- 0.1*10^9
  
  
  # Density unit: kg/m3
  # oil
  rho_oil <- 0.71*1000
  # water
  rho_brine <- 1.03*1000
  # gas
  rho_gas <- 0.29*1000
  
  
  rho_quartz <- 2.65*1000
  
  #cont_grain_sh <- 14
  #phi_b_cem <-0.2828692
  
  #phi_b_cem <- 0.92
  #phi_b_cem <- 0.94
  phi_b_cem <- 0.2956711
  
  # Depth around which the smoothing happens
  sm <- 0
  #K <- 2*sm
  # Linear step for the weighting of each of the models around the smoothing are, linear
  if(sm==0){step<-0}else{step <- 0.5/sm}
  
  # porosity model variables
  
  phi0_ss <- 0.4
  phi0_sh <- 0.6
  
  alpha_ss <- 0.32
  alpha_sh <- 0.6
  
  ref_depth <- 900*convrt
  
  #depth_cem <- depth_cem
  
  #kappa_ss <- 0.12
  kappa_ss <- 0.18
  
  #no_slip <- 0.1
  no_slip <- 0.4
  
  #no_slip_cem <- 0.9
  
  shear_mod <- matrix(NA, nrow = 4, ncol = tot_dim)
  alpha_cont_cem <- matrix(NA, nrow=1, ncol = tot_dim)
  # 4th row is to keep track of depth corresponding to each of the data points
  modulus_uncons_sand <- matrix(NA,  nrow = 2, ncol = tot_dim)
  mix_model <- matrix(NA,  nrow = 2, ncol = tot_dim)
  
  gassmann_dry_result_sat <- matrix(NA, nrow = 5, ncol = tot_dim)
  
  # make a vector of depths of the grid 
  
  #depth_grid <- matrix(gridded_ds$msh[[3]]*convrt,nrow=tot_dim, ncol =1, byrow = F)
  depth_grid <- matrix(traveltimes*convrt,nrow=tot_dim, ncol =1, byrow = F)
  
  # fill the 5th column of the vector with the depth of the grid (traveltimes)
  
  # Size 2, it has k_mixed and mu_mixed
  
  
  # index for elements
  i <- 1
  # linear factor for the region of mixture of the two models
  k <- 1
  
  for(i in 1:tot_dim){
    elem <- depth_grid[i]
    # laminated sand-shale model, didnt change the name
    k_mixed <- constant_clay_model(k_quartz,k_clay,mu_quartz,mu_clay,clay_content[i])
    
    if(elem=="NaN"){
      
    }else{
      #cont_grain <- 8
      cont_grain <- 20 - 34*poros_tot[i]+14*poros_tot[i]^2
      #depth_trend_el_par[,i] <- gassmann_fluid_sub(rho_1, v_p_1, v_s_1, k_mineral, k_fluid_1, k_fluid_2, rho_oil_1, rho_oil_2,poros_ss[1,i])
      k_fluid_eff <- 1/(sat_brine[i]/k_brine + sat_oil[i]/k_fluid_o + sat_gas[i]/k_fluid_g)
      # including the varying saturation
      rho_fluid_eff <- sat_brine[i]*rho_brine + sat_oil[i]*rho_oil + sat_gas[i]*rho_gas
      
      # rho_mixed in principle
      rho_mineral_eff <- clay_content[i]*rho_clay + (1 - clay_content[i])*rho_quartz
    }
    
    
    
    
    if(elem=="NaN"){
      # Set everything for this node to be NaN
      gassmann_dry_result_sat[1:4,i] <- NaN
      
      
      #print(i)
      
    }else if(elem<=((depth_cem+1)-sm)){
      #print(c(i,"elseif1"))
      mix_model[,i] <- uncons_sand(k_mixed[1], k_mixed[2], elem, rho_brine, rho_quartz,cont_grain, phi_ssc, poros_tot[i],no_slip)
      shear_mod[1,i] <- mix_model[2,i]
      
      # The Gassmann substitution dry -> mixed fluid, where the saturation AND mixed mineral is included
      gassmann_dry_result_sat[1:4,i] <- gassmann_dry_fluid_lam(k_mixed[1], poros_ss[i],poros_sh[i], k_fluid_eff, mix_model[1,i],rho_fluid_eff, rho_mineral_eff ,rho_clay, mix_model[2,i], clay_content[i], poros_tot[i])
      if(any(is.na(c(gassmann_dry_result_sat[1:4,i])))==TRUE){
        #browser()
        #debugonce()
        print(c(gassmann_dry_result_sat[1:4,i],"case1"))}
      
    }else if(elem>((depth_cem+1)-sm) & elem<=(depth_cem+sm)){
      #print(c(i,"elseif2"))
      # The smoothing part happens, transitioning between the two models, around the cementation point
      mix_model_1 <-(1-(k*step)) *uncons_sand(k_mixed[1], k_mixed[2], elem, rho_brine, rho_quartz,cont_grain,phi_ssc, poros_tot[i],no_slip)
      
      mix_model_2 <- contact_cement(poros_tot[i], phi_b_cem, k_mixed[1], k_mixed[2], k_mixed[1], k_mixed[2], cont_grain,no_slip)[1:2] + c(3.6*10^9,5.1*10^9)
      alpha_cont_cem[i]<- contact_cement(poros_tot[i], phi_b_cem, k_mixed[1], k_mixed[2], k_mixed[1], k_mixed[2], cont_grain,no_slip)[3]
      mix_model_2 <- (k*step)*mix_model_2
      if(any(is.na(c(mix_model_1,mix_model_2)))==TRUE){print(c(mix_model_1,mix_model_2))}
      mix_model[,i] <- mix_model_2 + mix_model_1
      shear_mod[1,i] <- mix_model[2,i]
      gassmann_dry_result_sat[1:4,i] <- gassmann_dry_fluid_lam(k_mixed[1], poros_ss[i],poros_sh[i],k_fluid_eff, mix_model[1,i],rho_fluid_eff, rho_quartz,rho_clay, mix_model[2,i], clay_content[i],poros_tot[i])
      if(any(is.na(c(gassmann_dry_result_sat[1:4,i])))==TRUE){print(c(gassmann_dry_result_sat[1:4,i]),"case2")}
      k <- k+1
      
    }else if(elem>(depth_cem+sm)){
      #print(c(i,"elseif3"))
      
      mix_model[,i] <- contact_cement(poros_tot[i], phi_b_cem, k_mixed[1],k_mixed[2], k_mixed[1], k_mixed[2], cont_grain,no_slip)[1:2] + c(3.6*10^9,0)# + c(3.6*10^9,5.1*10^9)
      shear_mod[1,i] <- mix_model[2,i]
      alpha_cont_cem[i]<- contact_cement(poros_tot[i], phi_b_cem, k_mixed[1], k_mixed[2], k_mixed[1], k_mixed[2], cont_grain,no_slip)[3]
      gassmann_dry_result_sat[1:4,i] <- gassmann_dry_fluid_lam(k_mixed[1], poros_ss[i], poros_sh[i], k_fluid_eff, mix_model[1,i],rho_fluid_eff, rho_quartz, rho_clay,mix_model[2,i], clay_content[i], poros_tot[i])
      if(any(is.na(c(gassmann_dry_result_sat[1:4,i])))==TRUE){print(gassmann_dry_result_sat[1:4,i])}
    }
    gassmann_dry_result_sat[5,i] <- elem
    
    
    i <- i+1
  }
  
  
  
  
  
  
 
  full_model_sat <- avo_r0_g(gassmann_dry_result_sat)
  
  return_list <- list()
  return_list[[1]] <- gassmann_dry_result_sat[1:3,]
  return_list[[2]] <- full_model_sat
  
  return(return_list)
  #return(gassmann_dry_result_sat[1,])
  
}





# Plotting 2.5D plots 

plotting_likelihood <- function(){
  
  xline<- gridded_ds$msh[[2]]
  inline <- gridded_ds$msh[[1]]
  
  gassmann_dry_result_vp_grid <- matrix(gassmann_dry_result_sat[1,], nrow = dim_grid[1], ncol = dim_grid[2], byrow = F)
  gassmann_dry_result_vs_grid <- matrix(gassmann_dry_result_sat[2,], nrow = dim_grid[1], ncol = dim_grid[2], byrow = F)
  gassmann_dry_result_rho_grid <- matrix(gassmann_dry_result_sat[3,], nrow = dim_grid[1], ncol = dim_grid[2], byrow = F)
  
  par(mfrow=c(2,2))
  image.plot(x = 1:dim_grid[2], y = 1:dim_grid[1], t(gassmann_dry_result_vp_grid), xlab = "inline", ylab="xline", main = expression("P-wave velocity"))
  image.plot(x = 1:dim_grid[2], y = 1:dim_grid[1], t(gassmann_dry_result_vs_grid), xlab = "inline", ylab="xline", main = expression("S-wave velocity"))
  image.plot(x = 1:dim_grid[2], y = 1:dim_grid[1], t(gassmann_dry_result_rho_grid), xlab = "inline", ylab="xline", main = expression("Density"))
  
  poros_grid <- matrix(poros_ss[1,], nrow = dim_grid[1], ncol = dim_grid[2], byrow = F)
  image.plot(x = 1:dim_grid[2], y = 1:dim_grid[1], t(poros_grid), xlab="inline", ylab="xline", main = expression("Porosity-mean"))
  
  poros_grid_sh <- matrix(poros_sh[1,], nrow = dim_grid[1], ncol = dim_grid[2], byrow = F)
  image.plot(x = 1:dim_grid[2], y = 1:dim_grid[1], t(poros_grid_sh), xlab="inline", ylab="xline", main = expression("Porosity shale -mean"))
  
  
  full_model_sat_grid_r0 <-  matrix(full_model_sat[1,], nrow = dim_grid[1], ncol = dim_grid[2], byrow = F)
  full_model_sat_grid_g <-  matrix(full_model_sat[2,], nrow = dim_grid[1], ncol = dim_grid[2], byrow = F)
  image.plot(x = 1:dim_grid[2], y = 1:dim_grid[1], t(full_model_sat_grid_r0), xlab = "inline", ylab="xline", main = expression("Reflection coefficient"))
  image.plot(x = 1:dim_grid[2], y = 1:dim_grid[1], t(full_model_sat_grid_g), xlab = "inline", ylab="xline", main = expression("Gradient"))
  
  
}



