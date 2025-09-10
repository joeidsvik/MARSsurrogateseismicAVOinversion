

# File to run MCMC with random walk proposal with either no covariance (q1) 
# or covariance Sigma from prior (q2). Both symmetric. 
# forward model: approximation using MARS


source("Alvheim_data.R")
source("Alvheim_prior.R")
source("Alvheim_likelihood.R")
source("Alvheim_plot_aux.R")

library(tictoc)
library(parallel)
library(Matrix)
library(nlme)
library(mcmcse)
library(mvnfast)

load("h_hat_models/mars_R0_bf.Rda")
load("h_hat_models/mars_G_bf.Rda")

# Proposal distributions ----

# Normal random walk with no covariance, q1 (symmetric) 
q1 <- function(x.t, s, N, tt_area = 0, cm.g = 0, cm.o = 0, cm.c = 0){
  x.prop <- cbind(x.t[,1] + s * rnorm(N), x.t[,2] + s * rnorm(N), x.t[,3] + s * rnorm(N))
  return(x.prop) 
}


# Normal random walk wits covariance, q2 (symmetric) 
q2 <- function(x.t, s, N, tt_area, cm.g, cm.o, cm.c){
  z.prop.list <- prior_gen_mean_tt(tt = tt_area, n_real=1, prior_mean = x.t-x.t, h = 1)
  x.prop.g <- x.t[,1] + s*z.prop.list$x_g
  x.prop.o <- x.t[,2] + s*z.prop.list$x_o
  x.prop.c <- x.t[,3] + s*z.prop.list$x_cl
  
  x.prop <- cbind(x.prop.g, x.prop.o, x.prop.c)
  return(x.prop)
}



# Markov chain Monte Carlo ----

MCMC.MARS.RW <- function(x0, xline.idx, inline.idx, IOmega, s, m, save.every, q, h.hat.R0, h.hat.G){
  set.seed(2)

  # Travel times at selected area
  tt_area <- traveltimes[xline.idx, inline.idx] 
  

  # Expected value of prior
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  
  # standard deviations for prior
  sd.g <- 2.83
  sd.o <- 2.68
  sd.c <- 1.64 

  # Torus for FFT
  
  cm.g <- build_base_mx(n_xline = 178, n_inline = 248, range_phi = 15, sig_std = sd.g, delta_x = 1, delta_y = 1, corr_func_type = "gauss")
  cm.o <- build_base_mx(n_xline = 178, n_inline = 248, range_phi = 15, sig_std = sd.o, delta_x = 1, delta_y = 1, corr_func_type = "gauss")
  cm.c <- build_base_mx(n_xline = 178, n_inline = 248, range_phi = 15, sig_std = sd.c, delta_x = 1, delta_y = 1, corr_func_type = "gauss")


  # Depths at selected area
  depth_area <- mu.prior[,4]
  
  
  # Wells 
  r <- length(xline.idx)
  k <- length(inline.idx)
  N=r*k

  # well 1 (gas) (484,4798) = (inline, xline)
  j.well.1 <- which(inline[1,] == 484)
  i.well.1 <- c(which(xline[,1] == 4796), which(xline[,1] == 4800))
  j.1 <- which(inline.idx == j.well.1)
  i.1 <- c(which(xline.idx == i.well.1[1]), which(xline.idx == i.well.1[2]))
  idx.well.1 <- c(((j.1-1)*r+i.1),((j.1)*r+i.1))

  # well 2 (gas) (716,4914)
  j.well.2 <- which(inline[1,] == 716)
  i.well.2 <- c(which(xline[,1] == 4912), which(xline[,1] == 4916))
  j.2 <- which(inline.idx == j.well.2)
  i.2 <- c(which(xline.idx == i.well.2[1]), which(xline.idx == i.well.2[2]))
  idx.well.2 <- c(((j.2-1)*r+i.2),((j.2)*r+i.2))

  # well 3 (gas) (1300,5150)
  j.well.3 <- which(inline[1,] == 1296) # highest inline
  i.well.3 <- c(which(xline[,1] == 5148), which(xline[,1] == 5152))
  j.3 <- which(inline.idx == j.well.3)
  i.3 <- c(which(xline.idx == i.well.3[1]), which(xline.idx == i.well.3[2]))
  idx.well.3 <- c(((j.3-1)*r+i.3),((j.3)*r+i.3))[1:2] 

  # well 4 (oil) (1026,4924)
  j.well.4 <- c(which(inline[1,] == 1024), which(inline[1,] == 1028))
  i.well.4 <- which(xline[,1] == 4924)
  j.4 <- c(which(inline.idx == j.well.4[1]), which(inline.idx == j.well.4[2]))
  i.4 <- which(xline.idx == i.well.4)
  idx.well.4 <- c(((j.4-1)*r+i.4),((j.4)*r+i.4))
  

  # Data R0 G data 
  R0.G <- getR0G(xline.idx, inline.idx)
  
  
  # log likelihood and  log prior 
  log.like.y <- log.likelihood.avo.mars(R0.G, x0, depth_area, IOmega, h.hat.R0, h.hat.G)
  log.like.w <- log.likelihood.wells(x0, xline.idx, inline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4)
  log.pri <- log.prior.total.area(x0, mu.prior, N, cm.g, cm.o, cm.c)
  
  l0 <- as.numeric(log.like.y + log.like.w + log.pri)
  
  j <- m/save.every 
  
  x.data.gas <- matrix(nrow = j, ncol = N)
  x.data.gas[1,] <- x0[,1]
  
  x.data.oil <- matrix(nrow = j, ncol = N)
  x.data.oil[1,] <- x0[,2] 
  
  x.data.clay <- matrix(nrow = j, ncol = N) 
  x.data.clay[1,] <- x0[,3]
  
  accepted <- 0
  
  x.t <- x0
  l.t <- l0
  
  idx <- 1
  tic("MCMC loop")
  start.time <- Sys.time()
  for(t in 1:m){ # 2:m){ # (m-1)){
    # get x.prop from proposal density q
    x.prop <- q(x.t, s, N, tt_area, cm.g, cm.o, cm.c)

    # log likelihood and log prior 
    log.like.y <- log.likelihood.avo.mars(R0.G, x.prop, depth_area, IOmega, h.hat.R0, h.hat.G)
    log.like.w <- log.likelihood.wells(x.prop, xline.idx, inline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4)
    log.pri <- log.prior.total.area(x.prop, mu.prior, N, cm.g, cm.o, cm.c)
    
    
    l.prop <- as.numeric(log.like.y + log.like.w + log.pri)
    
    # get acceptance probaility 
    alpha <- min(1, exp((l.prop - l.t))) # Symmetric 
    u <- runif(1)
    
    if (alpha > u){
      accepted = accepted + 1
      x.t <- x.prop
      l.t <- l.prop
      
      if (!(t %% save.every)){
        x.data.gas[idx,] <- x.prop[,1]
        x.data.oil[idx,] <- x.prop[,2]
        x.data.clay[idx,] <- x.prop[,3]
        idx <- idx + 1 
      }

    }
    else if (!(t %% save.every)){
      x.data.gas[idx,] <- x.t[,1]
      x.data.oil[idx,] <- x.t[,2]
      x.data.clay[idx,] <- x.t[,3]
      idx <- idx + 1
    }
    if(!(t %% 5000)) print(t)
  }
  diff.time <- difftime(Sys.time(), start.time, units = "secs")[[1]]
  data <- list("gas" = x.data.gas, "oil" = x.data.oil, "clay" = x.data.clay, "accepted" = accepted, "time" = diff.time)
  toc()
  return(data)
}
MCMC.MARS.RW <- compiler::cmpfun(MCMC.MARS.RW)


# Prepare ----
set.seed(2)

# all area (skipping some data points) ----
inline.idx = seq(1, 248, 1)
xline.idx = seq(1, 178, 1)


# Covariance matrix for likelihood
#tic("Omega")
#Omega0 <- create.Omega0()

#Omega <- create.Omega(Omega0, N = N) 
#IOmega <- solve(Omega)
#save(IOmega, file="IOmega.Rda")
load("IOmega.Rda")
#tic()


tt_area <- traveltimes[xline.idx, inline.idx] # Travel time of selected area
mu.prior <- match_traveltimes_mean_trend(tt_area) # Expected value of prior

# start in the mean value of prior
x0 <- cbind(mu.prior[,1] , mu.prior[,2], mu.prior[,3])


# Run MCMC ----

run.MCMC.1 <- function(){
  m <- 500000
  save.every <- 10 
  b <- m/save.every - (m/save.every)/2


  ### main MCMC run ###
  # s = 0.006 -> 26.8 %
  s = 0.00045

  system.time(
    d1 <- MCMC.MARS.RW(x0, xline.idx, inline.idx, IOmega, s, m, save.every, q1, h.hat.R0 = mars.R0.bf, h.hat.G = mars.G.bf)
  )
  print(paste("s = ", s))
  print(paste("acceptance rate = ",  d1$accepted/m * 100))

  pdf("MCMC_runs/plots/q1_trace_gas50.pdf", width=6, height=6)
  par(mar=c(5,5,4,2))
  plot(d1$gas[,50], type="l", xlab="", ylab="", cex.axis =1.3, cex.lab=1.4, ylim=c(-8,-4))
  dev.off()

  plot.Sg(d1, burn.in = b, filename="MCMC_runs/plots/q1_00045_gas.pdf") 
  plot.So(d1, burn.in = b, filename="MCMC_runs/plots/q1_00045_oil.pdf") 
  plot.Sb(d1, burn.in = b, filename="MCMC_runs/plots/q1_00045_brine.pdf") 
  plot.Vcl(d1, burn.in = b, filename="MCMC_runs/plots/q1_00045_clay.pdf") 

  filename <- "MCMC_runs/plots/q1_uncertainty_00045"
  percentile(d=d1, inline.idx, xline.idx, bi=b, p=80, filename=filename)

  #save(d1, file="MCMC_runs/total_area/RW_MARS.Rda")
}

run.MCMC.2 <- function(){
  m <- 200000 # 500000
  save.every <- 10 
  b <- m/save.every - (m/save.every)/2


  s = 0.0088
  system.time(
    d2 <- MCMC.MARS.RW(x0, xline.idx, inline.idx, IOmega, s, m, save.every, q2, h.hat.R0 = mars.R0.bf, h.hat.G = mars.G.bf)
  )
  print(paste("s = ", s))
  print(paste("acceptance rate = ",  d2$accepted/m * 100))

  pdf("MCMC_runs/plots/q2_trace_gas50.pdf", width=6, height=6)
  par(mar=c(5,5,4,2))
  plot(d2$gas[,50], type="l", xlab="", ylab="", cex.axis =1.3, cex.lab=1.4, ylim=c(-8,-4))
  dev.off()

  plot.Sg(d2, burn.in = b, filename="MCMC_runs/plots/q2_0088_gas.pdf") 
  plot.So(d2, burn.in = b, filename="MCMC_runs/plots/q2_0088_oil.pdf") 
  plot.Sb(d2, burn.in = b, filename="MCMC_runs/plots/q2_0088_brine.pdf") 
  plot.Vcl(d2, burn.in = b, filename="MCMC_runs/plots/q2_0088_clay.pdf") 

  filename <- "MCMC_runs/plots/q2_uncertainty_0088"
  percentile(d=d2, inline.idx, xline.idx, bi=b, p=80, filename=filename)

  #save(d2, file="MCMC_runs/total_area/RW2.Rda")
}

#run.MCMC.2()
#run.MCMC.1()



