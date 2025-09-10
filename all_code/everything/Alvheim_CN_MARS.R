

# File to run MCMC with CN proposal q3 
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

load("h_hat_models/mars_R0_bf.Rda")
load("h_hat_models/mars_G_bf.Rda")


q3 <- function(x.t, mu.prior, s, N, chol.Sig.g, chol.Sig.o, chol.Sig.c){
  # mean of proposed x
  mu.g <- mu.prior[,1] + sqrt(1-s^2)*(x.t[,1] - mu.prior[,1])
  mu.o <- mu.prior[,2] + sqrt(1-s^2)*(x.t[,2] - mu.prior[,2])
  mu.c <- mu.prior[,3] + sqrt(1-s^2)*(x.t[,3] - mu.prior[,3])
  
  # proposed x
  x.prop.g <- mu.g + s * chol.Sig.g %*% rmvn(N, 0, 1) 
  x.prop.o <- mu.o + s * chol.Sig.o %*% rmvn(N, 0, 1) 
  x.prop.c <- mu.c + s * chol.Sig.c %*% rmvn(N, 0, 1) 
  
  x.prop <- cbind(x.prop.g, x.prop.o, x.prop.c)
  
  mu.g <- mu.prior[,1] + sqrt(1-s^2)*(x.prop.g - mu.prior[,1])
  mu.o <- mu.prior[,2] + sqrt(1-s^2)*(x.prop.o - mu.prior[,2])
  mu.c <- mu.prior[,3] + sqrt(1-s^2)*(x.prop.c - mu.prior[,3])

  return("x.prop" = x.prop)
}

# Markov chain Monte Carlo ----

MCMC.MARS.CN <- function(x0, xline.idx, inline.idx, s, m, save.every, q, h.hat.R0, h.hat.G, LOAD = FALSE){
  set.seed(2)
  
  # Number of grid points
  ny <- length(xline.idx)
  nx <- length(inline.idx)
  N <- nx*ny
  
  
  # Travel times at selected area
  tt_area <- traveltimes[xline.idx, inline.idx] 
  
  # Expected value of prior
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  
  
  # Depths at selected area
  depth_area <- mu.prior[,4]
  
  
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  # well 4 (oil) (1026,4924)
  j.well.4 <- c(which(inline[1,] == 1024), which(inline[1,] == 1028))
  i.well.4 <- which(xline[,1] == 4924)
  j.4 <- c(which(inline.idx == j.well.4[1]), which(inline.idx == j.well.4[2]))
  i.4 <- which(xline.idx == i.well.4)
  idx.well.4 <- c(((j.4-1)*r+i.4),((j.4)*r+i.4))
  
  
  # Covariance matrix for prior and proposal
  Sigma.list = getSigmas(xline.idx, inline.idx, LOAD = LOAD)
  ISig.g = Sigma.list$Ig
  ISig.o = Sigma.list$Io
  ISig.c = Sigma.list$Ic
  chol.Sig.g = Sigma.list$chol.g
  chol.Sig.o = Sigma.list$chol.o
  chol.Sig.c = Sigma.list$chol.c
  
  
  # Data R0 G data 
  R0.G <- getR0G(xline.idx, inline.idx)
  
  
  # Covariance matrix for likelihood
  Omega0 <- create.Omega0()
  
  Omega <- create.Omega(Omega0, N = N) 
  IOmega <- solve(Omega)
  
  # log likelihood
  log.like.y <- log.likelihood.avo.mars(R0.G, x0, depth_area, IOmega, h.hat.R0, h.hat.G)
  log.like.w <- log.likelihood.well.4(x0, xline.idx, inline.idx, idx.well.4)
  
  l0 <- as.numeric(log.like.y + log.like.w)
  
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
  
  idx <- 2
  tic("MCMC loop")
  for(t in 2:(m-1)){
    # get x.prop from proposal density q
    x.prop <- q(x.t, mu.prior, s, N, chol.Sig.g, chol.Sig.o, chol.Sig.c)
    
    # log likelihood 
    log.like.y <- log.likelihood.avo.mars(R0.G, x.prop, depth_area, IOmega, h.hat.R0, h.hat.G)
    log.like.w <- log.likelihood.well.4(x.prop, xline.idx, inline.idx, idx.well.4) 
    
    l.prop <- as.numeric(log.like.y + log.like.w)
    
    # get acceptance probaility  
    alpha <- min(1, exp((l.prop - l.t))) 
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
  data <- list("gas" = x.data.gas, "oil" = x.data.oil, "clay" = x.data.clay, "accepted" = accepted)
  toc()
  return(data)
}
MCMC.MARS.CN <- compiler::cmpfun(MCMC.MARS.CN)


# Prepare ----
set.seed(2)

# 30x30 "area of interest" ----
inline.idx = seq(156, 215, 2)
xline.idx = seq(71, 130, 2)

N <- length(inline.idx) * length(xline.idx)

tt_area <- traveltimes[xline.idx, inline.idx] # Travel time of selected area
mu.prior <- match_traveltimes_mean_trend(tt_area) # Expected value of prior

# start in the mean value of prior
x0 <- cbind(mu.prior[,1] , mu.prior[,2], mu.prior[,3])


# Run MCMC ----

m <- 500000
save.every <- 10 # There is a bug when setting save.every = 1. 
m/save.every



### main MCMC run ###

system.time(
  d3 <- MCMC.MARS.CN(x0, xline.idx, inline.idx, s = 0.045, m, save.every, q3, h.hat.R0 = mars.R0.bf, h.hat.G = mars.G.bf, LOAD=TRUE)
)
d3$accepted/m
save(d3, file="MCMC_runs/small_area/CN_MARS_045.Rda")

#plot.So(d3, burn.in = b, filename="MCMC_runs/plots/q3_045.pdf") 

#pdf("MCMC_runs/plots/q3_trace.pdf")
#plot(d3$oil[,400], type = "l")
#dev.off()


# Visualization of MCMC runs ----
#b <- m/save.every - (m/save.every)/2

#plot.Sg(d3, burn.in = b) 
#plot.So(d3, burn.in = b) 
#plot.Sb(d3, burn.in = b)
#plot.Vcl(d3, burn.in = b)


#par(mfrow = c(2, 2))
#plot(d3$gas[,200], type = "l")
#plot(d3$gas[,400], type = "l")
#plot(d3$gas[,600], type = "l")
#plot(d3$gas[,800], type = "l")



