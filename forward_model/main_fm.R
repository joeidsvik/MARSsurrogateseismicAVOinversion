
### Packages necessary for the code and which parts for what, specify ####
library(R.matlab)
library(geoR)
library(MASS)
library(ggplot2)
library(fields)
library(akima)

options(error = recover)
set.seed(357)

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

traveltimes <-gridded_ds$msh[[3]][-(1:2),-(1:2)]
traveltimes <- traveltimes[-nrow(traveltimes), -ncol(traveltimes)]

# x-axis
inline <- gridded_ds$msh[[1]][-(1:2),-(1:2)]
inline <- inline[-nrow(inline), -ncol(inline)]

# y-axis
xline <- gridded_ds$msh[[2]][-(1:2),-(1:2)]
xline <-  xline[-nrow(xline), -ncol(xline)]





# if depth and not traveltime is to be plotted use conversion factor below
# Estimated conversion constant for timetravel to depth conversion
convrt <- 1.055

# Plotting of the traveltimes
par(mfrow=c(1,1), mar=c(4,4,2,2), oma =c(0,0,0,0))
image.plot(x=inline, y = -xline, traveltimes*convrt, xlab="inline", ylab="xline",legend.args=list( text="Depth[m]"), col = rev(viridis(256)), asp = 1, axes = F)
axis(1)
box()
ylabs <- axis(2, labels = FALSE, tick = FALSE)
axis(2, at = ylabs, labels=-(ylabs))

# Plotting of the data
par(mfrow=c(1,2), mar=c(4,4.5,2,1.5), oma =c(0,0,0,1))
image.plot(x = inline,y = -xline, r0_data,  main = "Zero-offset", xlab = "inline", ylab = "xline", col = viridis(256), asp = 1, axes = F)
box(lwd=0.1)
axis(1)
box()
ylabs <- axis(2, labels = FALSE, tick = FALSE)
axis(2, at = ylabs, labels=-(ylabs))
image.plot( x = inline,y = -xline, g_data, main = "Gradient",  xlab = "inline", ylab = "xline", col = viridis(256), asp = 1, axes = F)
axis(1)
box()
ylabs <- axis(2, labels = FALSE, tick = FALSE)
axis(2, at = ylabs, labels=-(ylabs))


# Example of a small region of traveltimes extracted (corresponds to depths of a region)
# after this step you can additionally extract an area from the traveltime matrix
r_area <- c(6:35)
c_area <- c(171:200)
tt_area <- traveltimes[r_area,c_area]


depth_test <- traveltimes*convrt

par(mfrow=c(2,1), mar=c(4,4.5,2,1.5), oma =c(0,0,0,1))

image.plot(x = inline[r_area,c_area], y = -xline[r_area,c_area], depth_test[r_area,c_area], col = rev(viridis(256)), zlim = c(min(depth_test),max(depth_test)))
# Negative xline value, only due to orientation of the y-axis in the figure
points(x = 1026, y = -4924, col = "#21918c", pch = 17, cex = 1.5)


####################
# Function call for finding the mean values corresponding to the traveltimes.
# Needs to include the prior depth trend file as well
# The returned variable is of dimension (n x m) x 3
# Where columns are 1 = x_g, 2 = x_o, 3 = v_clay
# remember the output is not on a 0-1 range. Apply the 
prior_mean <- match_traveltimes_mean_trend(tt_area)


# Input: 
# tt: the region of travel times wished to focus on 
# n_real: number of realisations of the prior desired
#prior_mean <- prior_gen_from_tt(tt = tt_area, n_real  = 3)

########### ########
# Function call with all the inputs and parameters already given beforehand.
# Input details
# traveltimes should be a matrix in the input [n x m] and should either be extracted from the data file, or have a sensible value if testing
# x_g, x_o and x_clay flat vectors of dimension [(n x m) x 1]
# Output
# (R0,G)
# numbers taken out as an example from a prior mean trend estimate for some depth

# If it is preferred to run the simulation of the forward model for multiple realisations of the prior
# then the for-loop can be implemented for the sim_forward_model
sim_data <- sim_forward_model(prior_mean$x_g[,1], prior_mean$x_o[,1], prior_mean$x_cl[,1], tt_area)



