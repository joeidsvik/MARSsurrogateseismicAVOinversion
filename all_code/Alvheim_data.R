

# Copy of Minas code in main_fm.R 


library(R.matlab)

library(geoR)

library(MASS)
#library(earth)
library(ggplot2)
library(fields)
library(akima)
#library(mvnfast)

source("rock_physics.R")
source("function_full_likelihood_laminated_shale_sand.R")
source("sim_fm.R")
source("function_fft2_rf_realis.R")


##### IF THE DATA IS TO BE LOADED #######
gridded_ds <- readMat("Alvheim_horizontal_classify_with_griddata.mat") 
# Loading the prior mean depth trend
depth_trend_prior <- read.table("AlvheimPrior.csv") 

# DONT CHANGE LOADING OF DATA BELOW UNLESS NAs purposefully being included in the code
# data: removing some extra rows and columns on top so that the number of rows and columns is a even number
# in addition to removing the rows including NA so that the inversion can be performed
r0_data <- gridded_ds$msh[[4]][-(1:2),-(1:2)] 
r0_data <- r0_data[-nrow(r0_data), -ncol(r0_data)] 
g_data <- gridded_ds$msh[[5]][-(1:2),-(1:2)] 
g_data <- g_data[-nrow(g_data), -ncol(g_data)] 
traveltimes <- gridded_ds$msh[[3]][-(1:2),-(1:2)]
traveltimes <- traveltimes[-nrow(traveltimes), -ncol(traveltimes)]

# x-axis
inline <- gridded_ds$msh[[1]][-(1:2),-(1:2)]
inline <- inline[-nrow(inline), -ncol(inline)]

# y-axis
xline <- gridded_ds$msh[[2]][-(1:2),-(1:2)]
xline <-  xline[-nrow(xline), -ncol(xline)]

