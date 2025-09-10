

# File to run MCMC with MALA 
# forward model: approximation using MARS
# derivatives approximated (directly) using MARS (not derivative of h_hat)


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
load("h_hat_models/D_mars.Rda")


# Proposal distribution (MALA med MARS) ----




q4 <- function(x.t, mu.prior, s, N, depth_area, IOmega, ISig.g, ISig.o, ISig.c, log.derivative.t, h.hat.R0, h.hat.G, D.hat, R0.G, idx.well.4){
  mu.g.t <- x.t[,1] + (0.5 * s^2) * log.derivative.t[1,] 
  mu.o.t <- x.t[,2] + (0.5 * s^2) * log.derivative.t[2,] 
  mu.c.t <- x.t[,3] + (0.5 * s^2) * log.derivative.t[3,]
  
  x.prop.g <- mu.g.t + s* rnorm(N)
  x.prop.o <- mu.o.t + s* rnorm(N)
  x.prop.c <- mu.c.t + s* rnorm(N)
  
  x.prop <- cbind(x.prop.g, x.prop.o, x.prop.c)
  
  # log(q(x.prop|x.t))
  q.prop <- -(0.5 / s^2) * (sum((x.prop.g-mu.g.t)^2) 
                            + sum((x.prop.o-mu.o.t)^2) 
                            + sum((x.prop.c-mu.c.t)^2)) 
  
  # nabla log(likelihood.avo(x.prop))
  log.derivative.likelihood.prop <- derivative.log.likelihood.avo(x.prop, depth_area, h.hat.R0, h.hat.G, D.hat, IOmega, R0.G)

  # nabla log(likelihood.well(x.prop))
  log.derivative.well.prop <- derivative.log.likelihood.well.4(x.prop, inline.idx, xline.idx, idx.well.4)

  # nabla log(prior(x.prop))
  log.derivative.prior.prop <- derivative.log.prior(x.prop, mu.prior, ISig.g, ISig.o, ISig.c)
  
  # nabla log(post(x.prop)) = log(likelihood(x.prop)) + log(prior(x.prop))
  log.derivative.prop <- log.derivative.likelihood.prop + log.derivative.prior.prop + log.derivative.well.prop

  mu.g.prop <- x.prop.g + (0.5 * s^2) * log.derivative.prop[1,]
  mu.o.prop <- x.prop.o + (0.5 * s^2) * log.derivative.prop[2,]
  mu.c.prop <- x.prop.c + (0.5 * s^2) * log.derivative.prop[3,]
  
  # log(q(x.t|x.prop))
  q.t <- - (0.5 / s^2) * (sum((x.t[,1]-mu.g.prop)^2) 
                         + sum((x.t[,2]-mu.o.prop)^2) 
                         + sum((x.t[,3]-mu.c.prop)^2) )
  
  return(list("x.prop" = x.prop, "q.t" = q.t, "q.prop" = q.prop))
}




# Markov chain Monte Carlo ----

MCMC.MARS.MALA <- function(x0, xline.idx, inline.idx, s, m, save.every, q, h.hat.R0, h.hat.G, D.hat, LOAD = FALSE){
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
  
  # log likelihood and log prior 
  log.like.avo <- log.likelihood.avo.and.derivative.mars(R0.G, x0, depth_area, IOmega, h.hat.R0, h.hat.G, D.hat)
  log.like.w <- log.likelihood.well.4.and.derivative(x0, xline.idx, inline.idx, idx.well.4)
  log.prior <- log.prior.and.derivative(x0, mu.prior, ISig.g, ISig.o, ISig.c)
    
  l0 <- as.numeric(log.like.avo$p + log.like.w$p + log.prior$p)
  
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

  log.derivative.t <- log.like.avo$d + log.like.w$d + log.prior$d 

  idx <- 2
  tic("MCMC loop")
  for(t in 2:(m-1)){
    # get x.prop from proposal density q
    langevin <- q(x.t, mu.prior, s, N, depth_area, IOmega, ISig.g, ISig.o, ISig.c, log.derivative.t, h.hat.R0, h.hat.G, D.hat, R0.G, idx.well.4)
    x.prop <- langevin$x.prop
    q.prop <- langevin$q.prop
    q.t <- langevin$q.t

    # log likelihood and log prior 
    log.like.avo <- log.likelihood.avo.and.derivative.mars(R0.G, x.prop, depth_area, IOmega, h.hat.R0, h.hat.G, D.hat)
    log.like.w <- log.likelihood.well.4.and.derivative(x.prop, xline.idx, inline.idx, idx.well.4)
    log.prior <- log.prior.and.derivative(x.prop, mu.prior, ISig.g, ISig.o, ISig.c)
    
    l.prop <- as.numeric(log.like.avo$p + log.like.w$p + log.prior$p)

    
    # get acceptance probaility 
    alpha <- min(1, exp((l.prop  + q.t - l.t - q.prop))) 
    u <- runif(1)
    
    if (alpha > u){
      accepted = accepted + 1
      x.t <- x.prop
      l.t <- l.prop
      log.derivative.t <- log.like.avo$d + log.like.w$d + log.prior$d 
      
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
    if(!(t %% 5)) print(t)
  }
  data <- list("gas" = x.data.gas, "oil" = x.data.oil, "clay" = x.data.clay, "accepted" = accepted)
  toc()
  return(data)
}
MCMC.MARS.MALA <- compiler::cmpfun(MCMC.MARS.MALA)


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

m <- 50
save.every <- 10 # There is a bug when setting save.every = 1. 
m/save.every




### main MCMC run ###
# LOAD = TRUE to use same sigmas again when running many times on same area

system.time(
  d4 <- MCMC.MARS.MALA(x0, xline.idx, inline.idx, s = 0.05, m, save.every, q = q4, h.hat.R0 = mars.R0.bf, h.hat.G = mars.G.bf, D.hat = D.mars, LOAD=TRUE)
)
d4$accepted/m

#save(d4, file = "MCMC_runs/small_area/MALA_05.Rda")
#plot.So(d4, burn.in = b, filename="MCMC_runs/plots/q4_05.pdf") 

#pdf("MCMC_runs/plots/q4_trace.pdf")
#plot(d4$oil[,400], type = "l")
#dev.off()


# Visualization of MCMC runs ----

#b <- m/save.every - (m/save.every)/2
#plot.Sg(d4, burn.in = b) 
#plot.So(d4, burn.in = b)
#plot.Sb(d4, burn.in = b)
#plot.Vcl(d4, burn.in = b)

#par(mfrow = c(2, 2))
#plot(d4$gas[,200], type = "l")
#plot(d4$gas[,400], type = "l")
#plot(d4$gas[,600], type = "l")
#plot(d4$gas[,800], type = "l")


