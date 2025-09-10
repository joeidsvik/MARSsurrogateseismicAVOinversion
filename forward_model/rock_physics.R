set.seed(3012)

# Drafting the forward model which includes an additional level with a rock physics model


# Simple Gassmann fluid substitution model 


#############
gassmann_fluid_sub <- function(rho_1, k_sat_1, mu_1, k_mineral, k_fluid_1, k_fluid_2, rho_fluid_1, rho_fluid_2,phi){

  const <- k_sat_1/(k_mineral - k_sat_1) - k_fluid_1/(phi*(k_mineral-k_fluid_1)) + k_fluid_2/(phi*(k_mineral-k_fluid_2))

  k_sat_2 <- const*k_mineral/(1+const)

  # Step 3
  mu_2 <- mu_1

  # Step 4

  rho_2 <- rho_1 + phi*(rho_fluid_2-rho_fluid_1)
  

  v_p_2 <- sqrt((k_sat_2 + (4/3)*mu_2)/rho_2)
  v_s_2 <- sqrt(mu_2/rho_2)




  return(c(v_p_2, v_s_2, rho_2, k_sat_2))
}
#################

# Gassmann's relations dry-fluid

gassmann_dry_fluid <- function(k_mineral,phi,k_fluid,k_dry, rho_fluid, rho_mineral, mu_dry){
  
  #const <- (k_dry*0.44)/(k_mineral-k_dry*0.44) + k_fluid/(phi*(k_mineral-k_fluid))
  const <- (k_dry)/(k_mineral-k_dry) + k_fluid/(phi*(k_mineral-k_fluid))
  
  k_sat <- k_mineral*const/(1+const)
  #mu_1 <- 44*10^9
  #mu_2 <- mu_1
  mu_2 <- mu_dry
  
  rho_s <- phi*rho_fluid + (1-phi)*rho_mineral
  
  v_p_s <- sqrt((k_sat + (4/3)*mu_2)/rho_s)
  v_s_s <- sqrt(mu_2/rho_s)
  #print(c(rho_s,k_sat,const,v_p_s))
  if(any(is.na(c(rho_s,k_sat,const,v_p_s)))==TRUE){
    print("NA encountered") 
    print(c(k_dry,k_fluid/(phi*(k_mineral-k_fluid)),k_sat,mu_2,rho_s,v_p_s,const))}
  
  return(c(v_p_s,v_s_s,rho_s, k_sat))
}

# Gassmann's relations dry-fluid: shaly sand case

gassmann_dry_fluid_sh_ss <- function(k_mineral,phi_sand, phi_shale,k_fluid,k_dry, rho_fluid, rho_qz, rho_cl, mu_dry, clay_cont,phi_tot){
  
  #const <- (k_dry*0.44)/(k_mineral-k_dry*0.44) + k_fluid/(phi*(k_mineral-k_fluid))
  const <- (k_dry)/(k_mineral-k_dry) + k_fluid/(phi_tot*(k_mineral-k_fluid))
  
  k_sat <- k_mineral*const/(1+const)
  #mu_1 <- 44*10^9
  #mu_2 <- mu_1
  mu_2 <- mu_dry
  
  rho_s <- (phi_sand-clay_cont*(1-phi_shale))*rho_fluid + (1-phi_sand)*rho_qz + clay_cont*(1-phi_shale)*rho_cl
  
  v_p_s <- sqrt((k_sat + (4/3)*mu_2)/rho_s)
  v_s_s <- sqrt(mu_2/rho_s)
  #print(c(rho_s,k_sat,const,v_p_s))
  if(any(is.na(c(rho_s,k_sat,const,v_p_s)))==TRUE){
    print("NA encountered") 
    print(c(k_dry,k_fluid/(phi*(k_mineral-k_fluid)),k_sat,mu_2,rho_s,v_p_s,const))}
  
  return(c(v_p_s,v_s_s,rho_s, k_sat))
}

gassmann_dry_fluid_lam <- function(k_mineral,phi_sand, phi_shale,k_fluid,k_dry, rho_fluid, rho_qz, rho_cl, mu_dry, clay_cont,phi_tot){
  
  #const <- (k_dry*0.44)/(k_mineral-k_dry*0.44) + k_fluid/(phi*(k_mineral-k_fluid))
  const <- (k_dry)/(k_mineral-k_dry) + k_fluid/(phi_tot*(k_mineral-k_fluid))
  
  k_sat <- k_mineral*const/(1+const)
  #mu_1 <- 44*10^9
  #mu_2 <- mu_1
  mu_2 <- mu_dry
  
  rho_s <- (rep(1, length(clay_cont))-clay_cont)*((1-phi_sand)*rho_qz + phi_sand*rho_fluid) + clay_cont*((1-phi_shale)*rho_cl + phi_shale*rho_fluid)
  
  
  v_p_s <- sqrt((k_sat + (4/3)*mu_2)/rho_s)
  v_s_s <- sqrt(mu_2/rho_s)
  #print(c(rho_s,k_sat,const,v_p_s))
  if(any(is.na(c(rho_s,k_sat,const,v_p_s)))==TRUE){
    #browser()
    #debugonce(contact_cement(poros_tot[i], phi_b_cem, k_mixed[1], k_mixed[2], k_mixed[1], k_mixed[2], cont_grain,no_slip)[1:2])
    print("NA encountered") 
    print(c(k_dry,k_fluid/(phi_tot*(k_mineral-k_fluid)),k_sat,mu_2,rho_s,v_p_s,const))}
  
  
  return(c(v_p_s,v_s_s,rho_s, k_sat))
}





# Friable(unconsolidated) sand model

uncons_sand <- function(k, mu, depth, rho_brine, rho_quartz, cont_grain, phi_c, phi,f_t){
  # f_t - no-slip factor
  # effective contact radius
  r_bar <- 1
  # spherical radius
  r <- 3
  
  # rho_brine is supposed to be rho_quartz
  pressure <- depth*9.81*(rho_quartz-rho_brine)
  ni <- (3*k-2*mu)/(2*(3*k+mu))
  
  k_hm <- (pressure*(cont_grain^2*(1-phi)^2*mu^2)/(18*pi^2*(1-ni)^2))^(1/3)#*(cont_grain^2*r_bar/r)^(1/3)
  
  # Old HM bound
  #mu_hm <- (5-4*ni)/(5*(2-ni))*(pressure*(3*cont_grain^2*(1-phi_c)^2*mu^2)/(2*pi^2*(1-ni)^2))^(1/3)
  
  # No-slip factor HM bound
  mu_hm <- (1/10*((12*(1-phi)^2*mu^2)/(pi^2*(1-ni)^2))^(1/3)*(cont_grain^2*r_bar/r)^(1/3)*pressure^(1/3))+(3/10*((12*(1-phi)^2*mu^2*(1-ni))/(pi^2*(2-ni)^3))^(1/3)*(cont_grain^2*r_bar/r)^(1/3)*pressure^(1/3))*f_t
  
  
  z_var <- mu_hm/6*((9*k_hm+8*mu_hm)/(k_hm+2*mu_hm))
  
  k_dry <- ((phi/phi_c)/(k_hm+4/3*mu_hm)+(1-phi/phi_c)/(k+4/3*mu_hm))^(-1)-4/3*mu_hm
  mu_dry <- ((phi/phi_c)/(mu_hm+z_var)+(1-phi/phi_c)/(mu+z_var))^(-1)-z_var
  
  
  
  return(c(k_dry,mu_dry))
}

# Constant cement model

constant_cement <- function(phi, phi_b, k_b, mu_b, k, mu){
  
  z_var <- (mu_b/6)*((9*k_b+8*mu_b)/(k_b+2*mu_b))
  
  k_dry <- ((phi/phi_b)/(k_b+(4/3)*mu_b)+(1-phi/phi_b)/(k+(4/3)*mu_b))^(-1)-(4/3)*mu_b
  
  mu_dry <- ((phi/phi_b)/(mu_b+z_var)+(1-phi/phi_b)/(mu + z_var))^(-1)-z_var
  
  return(c(k_dry, mu_dry))
}



# Contact cement model

contact_cement <- function(phi, phi_c, k_c, mu_c, k_s, mu_s, cont_grain,f_t){
  
  m_c <- k_c + 4*mu_c/3
  
  alpha <- (((2/3)*(phi_c-phi))/(1-phi_c))^0.5
  #alpha_2 <- 2*((phi_c-phi)/(3*cont_grain*(1-phi_c)))^(1/4)
  #alpha<- alpha_2
  
  ni_c <- 0.5*((k_c/mu_c)-2/3)/((k_c/mu_c)+1/3)
  
  ni_s <- 0.5*((k_s/mu_s)-2/3)/((k_s/mu_s)+1/3)
  
  lambda_n <- (2*mu_c*(1-ni_s)*(1-ni_c))/(pi*mu_s*(1-2*ni_c))
  
  c_n_lambda_n <- 0.00024649*lambda_n^(-1.9864)
  
  b_n_lambda_n <- 0.20405*lambda_n^(-0.89008)
  
  a_n_lambda_n <- -0.024153*lambda_n^(-1.3646)
  
  s_n <- a_n_lambda_n*alpha^2 + b_n_lambda_n*alpha + c_n_lambda_n
  
  lambda_tau <- mu_c/(pi*mu_s)
  
  c_tau_lambda_tau <- 10^(-4)*(9.654*ni_s^2+4.945*ni_s+3.1)*lambda_tau^(0.01867*ni_s^2+0.4011*ni_s-1.8186)
  
  b_tau_lambda_tau <- (0.0573*ni_s^2+0.0937*ni_s+0.202)*lambda_tau^(0.0274*ni_s^2+0.0529*ni_s-0.8765)
  
  a_tau_lambda_tau <- -10^(-2)*(2.26*ni_s^2+2.07*ni_s+2.3)*lambda_tau^(0.079*ni_s^2+0.1754*ni_s-1.342)
  
  s_tau <- a_tau_lambda_tau*alpha^2 + b_tau_lambda_tau*alpha^2 + c_tau_lambda_tau
  
  #s_tau <- s_tau*f_t
  
  k_dry <- (cont_grain*(1-phi_c)*m_c*s_n)/6
  
  mu_dry <- (3*k_dry)/5 + (3*cont_grain*(1-phi_c)*mu_c*s_tau)/20
  #print(c(k_dry, mu_dry,alpha))
  return(c(k_dry, mu_dry))
  
  
}

hertz_mindlin <- function(rho_mineral, rho_fluid, k, mu, cont_grain, phi_c, depth, phi){
  #rho_fluid <- 0
  pressure <- depth*9.81*(rho_mineral-rho_fluid)
  ni <- (3*k-2*mu)/(2*(3*k+mu))
  
  
  k_hm <- (pressure*(cont_grain^2*(1-phi_c)^2*mu^2)/(18*pi^2*(1-ni)^2))^(1/3)
  mu_hm <- (5-4*ni)/(5*(2-ni))*(pressure*(3*cont_grain^2*(1-phi_c)^2*mu^2)/(2*pi^2*(1-ni)^2))^(1/3)
  
  #mu_2 <- mu_hm
  #k_sat <- k_hm
  
  z_var <- mu_hm/6*((9*k_hm+8*mu_hm)/(k_hm+2*mu_hm))
  
  k_dry <- ((phi/phi_c)/(k_hm+4/3*mu_hm)+(1-phi/phi_c)/(k+4/3*mu_hm))^(-1)-4/3*mu_hm
  mu_dry <- ((phi/phi_c)/(mu_hm+z_var)+(1-phi/phi_c)/(mu+z_var))^(-1)-z_var
  k_mineral <- 25*10^9
  k_fluid <-  2.4*10^9
  
  const <- (k_dry)/(k_mineral-k_dry) + k_fluid/(phi*(k_mineral-k_fluid))
  
  k_sat <- k_mineral*const/(1+const)
  #mu_1 <- 44*10^9
  #mu_2 <- mu_1
  
  
  mu_2 <- mu_dry
  #k_sat <- k_dry
  
  rho_s <- phi*rho_fluid + (1-phi)*rho_mineral
  
  v_p_s <- sqrt((k_sat + (4/3)*mu_2)/rho_s)
  v_s_s <- sqrt(mu_2/rho_s)
  
  return(c(v_p_s, v_s_s, rho_s,k_sat))
}

constant_clay_model <- function(k_qz,k_clay,mu_qz,mu_clay,cem_frac){
  k_q <- k_qz + 4/3*mu_qz
  k_c <- k_clay + 4/3*mu_clay

  k_mix_inv <- (1-cem_frac)/k_q + cem_frac/k_c
  mu_mix_inv <- (1-cem_frac)/mu_qz + cem_frac/mu_clay
  
  return(c(1/k_mix_inv,1/mu_mix_inv))
}

porosity_sh <- function(phi_sh0, alpha_sh, depth_t, ref_depth){
  return(phi_sh0*exp(-alpha_sh*(depth_t-ref_depth)/1000))
}

porosity_ss <- function(phi_ss0,alpha_ss, depth_t, ref_depth, depth_tc, k_ss){
  if(depth_t <= depth_tc){
    #print(1)
    return(phi_ss0*exp(-alpha_ss*(depth_t-ref_depth)/1000))
  }else{
    
    return(phi_ss0*exp(-alpha_ss*(depth_tc-ref_depth)/1000)-k_ss*(depth_t-depth_tc)/1000)
  }
}


avo_r0_g <- function(elastic_p){
  refl_coeff <- matrix(NA, nrow=3, ncol=length(elastic_p[1,]))
  #v_p_caprock <- matrix(2592, nrow = length(elastic_p[1,]), ncol = 1)
  #v_s_caprock <- matrix(1150, nrow = length(elastic_p[1,]), ncol = 1)
  #rho_caprock <- matrix(2350, nrow = length(elastic_p[1,]), ncol=1)
  
  # 2650
  
  v_p_caprock <- matrix((2650), nrow = length(elastic_p[1,]), ncol = 1)
  v_s_caprock <- matrix((1150), nrow = length(elastic_p[1,]), ncol = 1)
  rho_caprock <- matrix((2350), nrow = length(elastic_p[1,]), ncol=1)
  
  
  delta_v_p <- elastic_p[1,] - v_p_caprock[,1]
  delta_v_s <- elastic_p[2,]- v_s_caprock[,1]
  delta_rho <- elastic_p[3,] - rho_caprock[,1]
  
  bar_v_p <- 0.5*(elastic_p[1,] + v_p_caprock[,1])
  bar_v_s <- 0.5*(elastic_p[2,]+ v_s_caprock[,1])
  bar_rho <- 0.5*(elastic_p[3,] + rho_caprock[,1])
  
  
  refl_coeff[1,] <- 0.5*(delta_v_p/bar_v_p + delta_rho/bar_rho)
  
  
  refl_coeff[2,] <- 0.5*(delta_v_p/bar_v_p) - 2*(bar_v_s^2/bar_v_p^2)*(delta_rho/bar_rho + 2*delta_v_s/bar_v_s)
  
  refl_coeff[3,] <- elastic_p[5,]
  nrow_dim <- length(refl_coeff[1,])+ length(refl_coeff[2,])
  ret_var <- matrix(NA, nrow=nrow_dim, ncol = 1)
  ret_var[seq(1,nrow_dim,2),1] <-refl_coeff[1,]
  ret_var[seq(2,nrow_dim,2),1] <- refl_coeff[2,]
 
  return(ret_var) 
}




