
### Packages necessary for the code and which parts for what, specify ####
library(R.matlab)

#library(geoR) erstattes av
library(sf)
library(sp)
library(mapview)
library(spdep)
library(rgdal)
library(spatialreg)
library(GWmodel)
library(spgwr)
library(gwrr)
library(tmap)

library(MASS)
library(ggplot2)
library(fields)
library(akima)

set.seed(357)

source("rock_physics.R")
source("function_full_likelihood_laminated_shale_sand.R")
source("sim_fm.R")


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

dim(g_data)

traveltimes <-gridded_ds$msh[[3]][-(1:2),-(1:2)]
traveltimes <- traveltimes[-nrow(traveltimes), -ncol(traveltimes)]

dim(traveltimes)

# x-axis
inline <- gridded_ds$msh[[1]][-(1:2),-(1:2)]
inline <- inline[-nrow(inline), -ncol(inline)]

# y-axis
xline <- gridded_ds$msh[[2]][-(1:2),-(1:2)]
xline <-  xline[-nrow(xline), -ncol(xline)]


# after this step you can additionally extract an area from the traveltime matrix
tt_area <- traveltimes[100:150, 1:50]


# if depth and not traveltime is to be plotted use conversion factor below
# Estimated conversion constant for timetravel to depth conversion
convrt <- 1.055

# Plotting of the traveltimes
plot.travel.times <- function(){
  par(mfrow=c(1,1), mar=c(4,4,2,2), oma =c(0,0,0,0))
  image.plot(x=inline, y = -xline, traveltimes*convrt, xlab="inline", ylab="xline",
             cex.lab=1.3, cex.axis=1.2, legend.width = 1.5, legend.shrink=1, axis.args=list(cex.axis=1.2),
             legend.args=list( text="Depth[m]"), col = rev(viridis(256)), asp = 1, axes = F)
  axis(1, cex.axis=1.2)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.2)
  #dev.off()
}

# Plotting of the data

plot.R0 <- function(){
  image.plot(x = inline,y = -xline, r0_data,  main = "", xlab = "inline", 
             ylab = "xline", col = viridis(256), cex.lab=1.4, cex.axis=1.3, 
             legend.width = 1.5, legend.shrink=1, axis.args=list(cex.axis=1.3), asp = 1, axes = F)
  box(lwd=0.1)
  axis(1, cex.axis=1.3)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.3)
  #dev.off()  
}


plot.G <- function(){
  image.plot( x = inline,y = -xline, g_data, main = "",  xlab = "inline", 
              ylab = "xline", col = viridis(256), cex.lab=1.4, cex.axis=1.3, 
              legend.width = 1.5, legend.shrink=1, axis.args=list(cex.axis=1.3), asp = 1, axes = F)
  axis(1, cex.axis=1.3)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.3)
  #dev.off()
}



####################
# Function call for finding the mean values corresponding to the traveltimes.
# Needs to include the prior depth trend file as well
# The returned variable is of dimension (n x m) x 4
# Where columns are 1 = x_g, 2 = x_o, 3 = v_clay, 4 = depth (here depth is convr*traveltimes)
#prior_mean <- match_traveltimes_mean_trend(tt_area)

########### ########
# Function call with all the inputs and parameters already given beforehand.
# Input details
# traveltimes should be a matrix in the input [n x m] and should either be extracted from the data file, or have a sensible value if testing
# x_g, x_o and x_clay flat vectors of dimension [(n x m) x 1]
# Output
# (R0,G)
# numbers taken out as an example from a prior mean trend estimate for some depth
#sim_data <- sim_forward_model(prior_mean[,1], prior_mean[,2], prior_mean[,3], tt_area)


