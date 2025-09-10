

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
library(mvnfast)

load("h_hat_models/mars_R0_bf.Rda")
load("h_hat_models/mars_G_bf.Rda")
load("h_hat_models/D_mars.Rda")



# Proposal distribution (MALA med MARS) ----




q4 <- function(x.t, mu.prior, s, N, depth_area, IOmega, log.derivative.t, h.hat.R0, h.hat.G, D.hat, R0.G, cm.g, cm.o, cm.c, idx.well.1, idx.well.2, idx.well.3, idx.well.4){
  # mean of proposed x
  mu.g.t <- x.t[,1] + (0.5 * s^2) * log.derivative.t[1,] 
  mu.o.t <- x.t[,2] + (0.5 * s^2) * log.derivative.t[2,] 
  mu.c.t <- x.t[,3] + (0.5 * s^2) * log.derivative.t[3,]
  
  # proposed x
  z.prop.list <- prior_gen_mean_tt(tt_area, n_real=1, mu.prior - mu.prior, h = 1)
  
  x.prop.g <- mu.g.t + s*z.prop.list$x_g
  x.prop.o <- mu.o.t + s*z.prop.list$x_o
  x.prop.c <- mu.c.t + s*z.prop.list$x_cl
  
  x.prop <- cbind(x.prop.g, x.prop.o, x.prop.c)
  
  # log(q(x.prop|x.t))
  q.prop <- -(0.5 / s^2) * (sum((x.prop.g-mu.g.t)^2) 
                            + sum((x.prop.o-mu.o.t)^2) 
                            + sum((x.prop.c-mu.c.t)^2)) 
  
  # nabla log(likelihood.avo(x.prop))
  log.derivative.likelihood.prop <- derivative.log.likelihood.avo(x.prop, N, depth_area, h.hat.R0, h.hat.G, D.hat, IOmega, R0.G)
  
  # nabla log(likelihood.well(x.prop))
  log.derivative.well.prop <- derivative.log.likelihood.well(x.prop, inline.idx, xline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4)

  # nabla log(prior(x.prop))
  log.derivative.prior.prop <- derivative.log.prior.total.area(x.prop, mu.prior, N, cm.g, cm.o, cm.c)
  
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

MCMC.MARS.MALA <- function(x0, xline.idx, inline.idx, s, m, save.every, q, IOmega, h.hat.R0, h.hat.G, D.hat){
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
  
  # log likelihood and log prior 
  log.like.avo <- log.likelihood.avo.and.derivative.mars(R0.G, x0, N, depth_area, IOmega, h.hat.R0, h.hat.G, D.hat)
  log.like.w <- log.likelihood.and.deriavtive.wells(x0, xline.idx, inline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4)
  log.prior <- log.prior.and.derivative.total.area(x0, mu.prior, N, cm.g, cm.o, cm.c)
  

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
  

  idx <- 1
  tic("MCMC loop")
  start.time <- Sys.time()
  for(t in 1:m){ # 2:m){ # (m-1)){
    # get x.prop from proposal density q
    langevin <- q(x.t, mu.prior, s, N, depth_area, IOmega, log.derivative.t, h.hat.R0, h.hat.G, D.hat, R0.G, cm.g, cm.o, cm.c, idx.well.1, idx.well.2, idx.well.3, idx.well.4)
    x.prop <- langevin$x.prop
    q.prop <- langevin$q.prop
    q.t <- langevin$q.t
    
    # log likelihood and log prior 
    log.like.avo <- log.likelihood.avo.and.derivative.mars(R0.G, x.prop, N, depth_area, IOmega, h.hat.R0, h.hat.G, D.hat)
    log.like.w <- log.likelihood.and.deriavtive.wells(x.prop, xline.idx, inline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4)
    log.prior <- log.prior.and.derivative.total.area(x.prop, mu.prior, N, cm.g, cm.o, cm.c)
    
    print(paste("# nan in prior = ", sum(is.na(log.prior$p))))
    print(paste("# nan in derivative prior = ", sum(is.na(log.prior$d))))

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
    if(!(t %% 5000)) print(t)
  }
  diff.time <- difftime(Sys.time(), start.time, units = "secs")[[1]]
  data <- list("gas" = x.data.gas, "oil" = x.data.oil, "clay" = x.data.clay, "accepted" = accepted, "time" = diff.time)
  toc()
  return(data)
}
MCMC.MARS.MALA <- compiler::cmpfun(MCMC.MARS.MALA)


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

run.MCMC <- function(){
  m <- 260000 # 500000
  save.every <-  10 
  b <- m/save.every - (m/save.every)/2

  # s = 0.001 -> 70%
  # s = 0.0015 -> ?%
  # s = 0.002 -> 5.6%
  s = 0.002

  system.time(
    d4 <- MCMC.MARS.MALA(x0, xline.idx, inline.idx, s, m, save.every, q = q4, IOmega, h.hat.R0 = mars.R0.bf, h.hat.G = mars.G.bf, D.hat = D.mars)
  )

  print(paste("s = ", s))

  print(d4$accepted/m * 100)

  plot.Sg(d4, burn.in = b, filename="MCMC_runs/plots/q4_002_gas.pdf") 
  plot.So(d4, burn.in = b, filename="MCMC_runs/plots/q4_002_oil.pdf") 
  plot.Sb(d4, burn.in = b, filename="MCMC_runs/plots/q4_002_brine.pdf") 
  plot.Vcl(d4, burn.in = b, filename="MCMC_runs/plots/q4_002_clay.pdf") 

  filename <- "MCMC_runs/plots/q4_uncertainty_002"
  percentile(d=d4, inline.idx, xline.idx, bi=b, p=80, filename=filename)

  pdf("MCMC_runs/plots/q4_trace_gas.pdf", width=6, height=6)
  par(mar=c(5,5,4,2))
  plot(d4$gas[,50], type="l", xlab="", ylab="", cex.axis =1.3, cex.lab=1.4, ylim=c(-8,-4))
  dev.off()

  #save(d4, file="MCMC_runs/total_area/CN_MARS_05.Rda")

}


#run.MCMC()






