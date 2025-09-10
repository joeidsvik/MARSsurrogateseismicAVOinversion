#### Ternary plots ####
# install.packages("Ternary")
library(Ternary)
library(R.matlab)

load("benchmark_case_no_obs_noise.RData")


# Additional variables needed to run the file independently of other files
gridded_ds <- readMat("Alvheim_horizontal_classify_with_griddata.mat")

# data: removing some extra rows and columns on top so that the number of rows and columns is a even number
# in addition to removing the rows including NA so that the inversion can be performed
r0_data <- gridded_ds$msh[[4]][-(1:2),-(1:2)]
r0_data <- r0_data[-nrow(r0_data), -ncol(r0_data)]
g_data <- gridded_ds$msh[[5]][-(1:2),-(1:2)]
g_data <- g_data[-nrow(g_data), -ncol(g_data)]

traveltimes <-gridded_ds$msh[[3]][-(1:2),-(1:2)]
traveltimes <- traveltimes[-nrow(traveltimes), -ncol(traveltimes)]

inline <- gridded_ds$msh[[1]][-(1:2),-(1:2)]
inline <- inline[-nrow(inline), -ncol(inline)]

xline <- gridded_ds$msh[[2]][-(1:2),-(1:2)]
xline <-  xline[-nrow(xline), -ncol(xline)]

# full dimension of the grid
tot_dim <- prod(dim(traveltimes))
# dimension of the full grid of format [n m]
dim_grid <- dim(traveltimes)
# inline
x_grid_full <- inline
# xline
y_grid_full <- xline
# xline - y axis - row
nt <- dim_grid[1]

# inline - x axis - column
nx <- dim_grid[2]
# conversion constant for timetravel to depth conversion
convrt <- 1.055
# base case: 2000
# cementation depth: 1950,2020,2050 - maybe change to 1980??
depth_cem <- 2000*convrt
#depth_cem <- 2050*convrt
#ref_depth <- 0
ref_depth <- 900*convrt
#### Prior inputs ####
#phi_cl <- 15
#phi_sat <- 15

phi_cl <- 15
phi_sat <- 15

depth_full <- traveltimes*convrt
#default
uncond <- FALSE
#uncond <- TRUE
depth_trend_prior <- read.table("AlvheimPrior.csv")
write_to_file <- FALSE
# Whether prior realisations should be generated with spatial correlation
corr <- TRUE
# observation error
cov_obs <- matrix(c(0.003,-0.006,-0.006,0.03), nrow=2, ncol=2, byrow = F)
t_chol <- t(chol(cov_obs))
Rloc <- cov_obs




#### LETKF constants #####

# number of ensembles
n_e <- 100

maxIter <- 2
# number of parameters starting out with, which are inverted for
n_param <- 3
# number of observations at each location
n_obs <- 2

# Layers in the model
nLayer <-  1

# Defining the size of patch for the variables and observations
#patch_size <- 6
#obs_patch <- 16
patch_size <- 6
obs_patch <- 16

# enkf <- 1 ; ienks
# enkf <- 2 ; etkf

enkf <- 2

#### Grid reshaping and localisation ####
# Number of nodes in one direction - for one parameter
n_nodes_col <- dim_grid[2]-(obs_patch-patch_size)
n_nodes_row <- dim_grid[1] - (obs_patch - patch_size)

# If we want to look at a square grid
#n_nodes <- min(n_nodes_row, n_nodes_col)

# Nodes in which the observations are made
# n_nodes <- dim_grid[1]

n_nodes_obs_row <- dim_grid[1]
n_nodes_obs_col <- dim_grid[2]

# need to be generalized
diff_nodes <- (obs_patch-patch_size)/2

# calculating different sizes for the sake of simplicity
tot_dim_loc <- n_nodes_col*n_nodes_row

#dim_grid_loc <- n_nodes_obs
tot_dim_obs <- n_nodes_obs_row*n_nodes_obs_col

# Dimension variables for flexibility
ncol_grid <- n_e
nrow_one_param <- n_nodes_row*n_nodes_col
nrow_one_obs <- n_nodes_obs_row*n_nodes_obs_col
nrow_grid <- nrow_one_param*n_param



x_rf_sat_g_log_loc <- prior_real_list[[3]]
x_rf_sat_o_log_loc <- prior_real_list[[2]]
x_rf_sat_b_log_loc <- prior_real_list[[19]]
x_rf_clay_log_loc <- prior_real_list[[4]]


ens_grid_sat_g <- post_list[[1]]
ens_grid_sat_o <- post_list[[2]]
ens_grid_clay <- post_list[[3]]



# Flattening out the ensembles before 
ens_sat_g_log_flat <- matrix(NA, nrow = (n_nodes_row*n_nodes_col), ncol = n_e)
ens_sat_o_log_flat <- matrix(NA, nrow = (n_nodes_row*n_nodes_col), ncol = n_e)
ens_sat_b_log_flat <- matrix(NA, nrow = (n_nodes_row*n_nodes_col), ncol = n_e)
ens_clay_log_flat <- matrix(NA, nrow = (n_nodes_row*n_nodes_col), ncol = n_e)

ens_sat_g_flat <- matrix(NA, nrow = (n_nodes_row*n_nodes_col), ncol = n_e)
ens_sat_o_flat <- matrix(NA, nrow = (n_nodes_row*n_nodes_col), ncol = n_e)


ens_clay_flat <- matrix(NA, nrow = (n_nodes_row*n_nodes_col), ncol = n_e)




for(e in 1:n_e){
  
  ens_sat_g_flat[,e] <- matrix(ens_grid_sat_g[[e]], nrow=(n_nodes_col*n_nodes_row),ncol=1, byrow = F)
  ens_sat_o_flat[,e] <- matrix(ens_grid_sat_o[[e]], nrow=(n_nodes_col*n_nodes_row),ncol=1, byrow = F)
  ens_clay_flat[,e] <- matrix(ens_grid_clay[[e]], nrow=(n_nodes_col*n_nodes_row),ncol=1, byrow = F)
  
  
  
  sum_tot <- 1+ exp(ens_grid_sat_g[[e]]) + exp(ens_grid_sat_o[[e]])
  sat_g_grid_log_ens <- exp(ens_grid_sat_g[[e]])/sum_tot
  sat_o_grid_log_ens <- exp(ens_grid_sat_o[[e]])/sum_tot
  sat_b_grid_log_ens <- 1/sum_tot
  
  #poros_grid_log_ens <- (exp(ens_grid_poros[[1]])/(1+exp(ens_grid_poros[[1]]))*por_max + 1/(1+exp(ens_grid_poros[[1]]))*por_min)
  clay_grid_log_ens <- exp(ens_grid_clay[[e]])/(1 + exp(ens_grid_clay[[e]]))
  
  ens_sat_g_log_flat[,e] <- matrix(sat_g_grid_log_ens, nrow=(n_nodes_row*n_nodes_col),ncol=1, byrow = F)
  ens_sat_o_log_flat[,e] <- matrix(sat_o_grid_log_ens, nrow = (n_nodes_row*n_nodes_col), ncol = 1, byrow = F)
  ens_sat_b_log_flat[,e] <- matrix(sat_b_grid_log_ens, nrow = (n_nodes_row*n_nodes_col), ncol = 1, byrow = F)
  
  ens_clay_log_flat[,e] <- matrix(clay_grid_log_ens, nrow = (n_nodes_row*n_nodes_col), ncol = 1, byrow = F)
  
}


# Saving data for the ternary plots without correlation
#save(x_rf_sat_g_log_loc, x_rf_sat_o_log_loc, x_rf_sat_b_log_loc, x_rf_clay_log_loc, ens_clay_log_flat, ens_sat_g_log_flat,ens_sat_o_log_flat, ens_sat_b_log_flat , file = "ternary_no_corr.RData")
# Saving data for the ternary plots for benchmark case
##save(x_rf_sat_g_log_loc, x_rf_sat_o_log_loc, x_rf_sat_b_log_loc, x_rf_clay_log_loc, ens_clay_log_flat, ens_sat_g_log_flat,ens_sat_o_log_flat, ens_sat_b_log_flat , file = "ternary_benchmark.RData")



#s_indx <- 16875
#s_indx <- 12533

# inline436 xline4532 , [2,28]
# deep location index
d_indx <- 4538
#(69,119), inline 800, xline 4800
#(118*168+69), 
# shallow location index
s_indx <- 19893
# inline 412 xline 4800, , [69,22]
# s_indx <- 3597



# Prior and psoterior, shallow and deep locations #
#### Prior and posterior ternary plots for benchmark case ####

# Prior ternary plots
setEPS()
postscript(paste0("minor_revision_figures/ternary_prior_new_s",s_indx,"_d",d_indx,".eps", collapse=","), width = 12, height = 6.9, pointsize = 13, colormodel = "cmyk")
#png(paste0("minor_revision_figures/ternary_prior_new_s",s_indx,"_d",d_indx,".png", collapse=","),units="px", width = 2083, height=1200, pointsize = 35)
par(mfrow=c(1,2),mar =c(1, 1, 0, 1) + 0.1)
# indexes: depth_loc[75,101], xline 4824, inline 728
# Calculated 16785-(floor(16785/168)*168), ceil(16785/168)
# 2124.398
tern_sat_g <- x_rf_sat_g_log_loc[s_indx,]
tern_sat_o <- x_rf_sat_o_log_loc[s_indx,]
tern_sat_b <-  x_rf_sat_b_log_loc[s_indx,]


tern_sat <- cbind(tern_sat_g, tern_sat_o, tern_sat_b)

TernaryPlot(axis.labels = seq(0, 1, by = 0.1), alab = "gas", blab="oil", clab="brine", axis.cex = c(1.3,1.3,1.3), lab.cex=c(1.5,1.5,1.5))
TernaryPoints(tern_sat, cex = 1.1, col = "black", lwd = 3)
mtext(~ bold("a)"), side = 3, line = 0, padj = 1, adj = 0, cex = 1.4)

# depth: 2200.892 
# indexes: depth_loc[2,28]

tern_sat_g <- x_rf_sat_g_log_loc[d_indx,]
tern_sat_o <- x_rf_sat_o_log_loc[d_indx,]
tern_sat_b <-  x_rf_sat_b_log_loc[d_indx,]


tern_sat <- cbind(tern_sat_g, tern_sat_o, tern_sat_b)

TernaryPlot(axis.labels = seq(0, 1, by = 0.1), alab = "gas", blab="oil", clab="brine", axis.cex = c(1.3,1.3,1.3), lab.cex=c(1.5,1.5,1.5))
TernaryPoints(tern_sat, cex = 1.1, col = "black", lwd = 3)
mtext(~ bold("b)"), side = 3, line = 0, padj = 1, adj = 0, cex = 1.4)

dev.off()






### Posterior ternary plots ######
setEPS()
postscript(paste0("minor_revision_figures/ternary_post_new_s",s_indx,"_d",d_indx,".eps", collapse=","),width = 12, height = 6.9, pointsize = 13, colormodel = "cmyk")
#png(paste0("minor_revision_figures/ternary_post_new_s",s_indx,"_d",d_indx,".png", collapse=","),units="px", width = 2083, height=1200, pointsize = 35)
par(mfrow=c(1,2),mar =c(1, 1, 0, 1) + 0.1)

# 2124.398
tern_sat_g <- ens_sat_g_log_flat[s_indx,]
tern_sat_o <- ens_sat_o_log_flat[s_indx,]
tern_sat_b <-  ens_sat_b_log_flat[s_indx,]


tern_sat <- cbind(tern_sat_g, tern_sat_o, tern_sat_b)

TernaryPlot(axis.labels = seq(0, 1, by = 0.1), alab = "gas", blab="oil", clab="brine", axis.cex = c(1.3,1.3,1.3), lab.cex=c(1.5,1.5,1.5))
TernaryPoints(tern_sat, cex = 1.1, col = "black", lwd = 3)
mtext(~ bold("a)"), side = 3, line = 0, padj = 1, adj = 0, cex = 1.4)

# depth: 2200.892 

tern_sat_g <- ens_sat_g_log_flat[d_indx,]
tern_sat_o <- ens_sat_o_log_flat[d_indx,]
tern_sat_b <- ens_sat_b_log_flat[d_indx,]


tern_sat <- cbind(tern_sat_g, tern_sat_o, tern_sat_b)

TernaryPlot(axis.labels = seq(0, 1, by = 0.1), alab = "gas", blab="oil", clab="brine", axis.cex = c(1.3,1.3,1.3), lab.cex=c(1.5,1.5,1.5))
TernaryPoints(tern_sat, cex = 1.1, col = "black", lwd = 3)
mtext(~ bold("b)"), side = 3, line = 0, padj = 1, adj = 0, cex = 1.4)

dev.off()


#load("benchmark_case.RData")


# New plots with new figures for the paper

#[99.174]
new_index <- 98*178+174
new_index_letkf <- 98*168+174
mcmc_sat_g <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_g_17618.rds")
mcmc_sat_o <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_o_17618.rds")


#[95,174]
new_index <- 94*178+174
new_index_letkf <- 94*168+174
mcmc_sat_g <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_g_16906.rds")
mcmc_sat_o <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_o_16906.rds")


#[[95,199]]

new_index <- 94*178+199
new_index_letkf <- 94*168+199

mcmc_sat_g <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_g_16931.rds")
mcmc_sat_o <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_o_16931.rds")

#[103,232] *
new_index <- 102*178+232
new_index_letkf <- 102*168+232

mcmc_sat_g_18388 <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_g_18388.rds")
mcmc_sat_o_18388 <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_o_18388.rds")

#[148,174]
new_index <- 148*178+174
new_index_letkf <- 148*168+174

mcmc_sat_g_26518 <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_g_26518.rds")
mcmc_sat_o_26518 <- readRDS("/Volumes/work/minasp/work_w_Karen_TT/mcmc_sat_o_26518.rds")



transf_func_o <- function(sat_g,sat_o){
  sum_tot <- 1 + exp(sat_g) + exp(sat_o)
  log_o <- exp(sat_o)/sum_tot
  return(log_o)
}

transf_func_g <- function(sat_g,sat_o){
  sum_tot <- 1 + exp(sat_g) + exp(sat_o)
  log_g <- exp(sat_g)/sum_tot
  return(log_g)
}
new_index <- 148*178+174
new_index_letkf <- 148*168+174

mcmc_sat_g <- mcmc_sat_g_26518
mcmc_sat_o <- mcmc_sat_o_26518


tern_log_sat_g_mcmc <- mapply(transf_func_g, mcmc_sat_g, mcmc_sat_o)
tern_log_sat_o_mcmc <- mapply(transf_func_o, mcmc_sat_g, mcmc_sat_o)
tern_log_sat_b_mcmc <- (1 - tern_log_sat_o_mcmc - tern_log_sat_g_mcmc)
# transform the variables and make a brine vector

tern_sat_mcmc <- cbind(tern_log_sat_g_mcmc, tern_log_sat_o_mcmc, tern_log_sat_b_mcmc)

new_tern_sat_g_letkf <- ens_sat_g_log_flat[new_index_letkf,]
new_tern_sat_o_letkf <- ens_sat_o_log_flat[new_index_letkf,]
new_tern_sat_b_letkf <- ens_sat_b_log_flat[new_index_letkf,]


tern_sat_letkf <- cbind(new_tern_sat_g_letkf, new_tern_sat_o_letkf, new_tern_sat_b_letkf)


setEPS()
postscript(paste0("/Users/minasp/Library/CloudStorage/OneDrive-NTNU/work with Mina/latex/Images/ternary_post_new_mcmc",new_index,".eps", collapse=","),width = 6.9, height = 6.9, pointsize = 13, colormodel = "cmyk")
#png(paste0("minor_revision_figures/ternary_post_new_s",s_indx,"_d",d_indx,".png", collapse=","),units="px", width = 2083, height=1200, pointsize = 35)
#par(mfrow=c(1,2),mar =c(1, 1, 0, 1) + 0.1)


TernaryPlot(axis.labels = seq(0, 1, by = 0.1), alab = "gas", blab="oil", clab="brine", axis.cex = c(1.3,1.3,1.3), lab.cex=c(1.5,1.5,1.5))
TernaryPoints(tern_sat_mcmc, cex = 1.1, col = "black", lwd = 3)
#mtext(~ bold("a)"), side = 3, line = 0, padj = 1, adj = 0, cex = 1.4)
dev.off()

setEPS()
postscript(paste0("/Users/minasp/Library/CloudStorage/OneDrive-NTNU/work with Mina/latex/Images/ternary_post_new_enkf",new_index,".eps", collapse=","),width = 6.9, height = 6.9, pointsize = 13, colormodel = "cmyk")
TernaryPlot(axis.labels = seq(0, 1, by = 0.1), alab = "gas", blab="oil", clab="brine", axis.cex = c(1.3,1.3,1.3), lab.cex=c(1.5,1.5,1.5))
TernaryPoints(tern_sat_letkf, cex = 1.1, col = "black", lwd = 3)
#mtext(~ bold("b)"), side = 3, line = 0, padj = 1, adj = 0, cex = 1.4)

dev.off()


