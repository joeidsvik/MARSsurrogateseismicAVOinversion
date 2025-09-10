


# File to run MCMC with CN proposal with covariance Sigma from prior (q2). 
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

q3 <- function(x.t, mu.prior, s, N, tt_area){
  # mean of proposed x
  mu.g <- mu.prior[,1] + sqrt(1-s^2)*(x.t[,1] - mu.prior[,1])
  mu.o <- mu.prior[,2] + sqrt(1-s^2)*(x.t[,2] - mu.prior[,2])
  mu.c <- mu.prior[,3] + sqrt(1-s^2)*(x.t[,3] - mu.prior[,3])
  
  # proposed x
  z.prop.list <- prior_gen_mean_tt(tt_area, n_real=1, mu.prior - mu.prior, h = 1)
  x.prop.g <- mu.g + s*z.prop.list$x_g
  x.prop.o <- mu.o + s*z.prop.list$x_o
  x.prop.c <- mu.c + s*z.prop.list$x_cl
  
  x.prop <- cbind(x.prop.g, x.prop.o, x.prop.c)
  
  return(x.prop)
}

  


# Markov chain Monte Carlo ----

MCMC.MARS.CN <- function(x0, xline.idx, inline.idx, s, m, save.every, IOmega, q, h.hat.R0, h.hat.G){
  set.seed(2)

  # Travel times at selected area
  tt_area <- traveltimes[xline.idx, inline.idx] 
  

  # Expected value of prior
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  
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
  
  idx <- 1
  tic("MCMC loop")
  start.time <- Sys.time()
  for(t in 1:m){ # 2:m){ # (m-1)){
    # get x.prop from proposal density q
    x.prop <- q(x.t, mu.prior, s, N, tt_area)
    
    # log likelihood and log prior 
    log.like.y <- log.likelihood.avo.mars(R0.G, x.prop, depth_area, IOmega, h.hat.R0, h.hat.G)
    log.like.w <- log.likelihood.wells(x.prop, xline.idx, inline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4) # total area
    
    #l.prop <- as.numeric(log.like.y + log.like.w + log.pri)
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
  diff.time <- difftime(Sys.time(), start.time, units = "secs")[[1]]
  data <- list("gas" = x.data.gas, "oil" = x.data.oil, "clay" = x.data.clay, "accepted" = accepted, "time" = diff.time)
  toc()
  return(data)
}
MCMC.MARS.CN <- compiler::cmpfun(MCMC.MARS.CN)


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

m <- 500000
save.every <- 10 
b <- m/save.every - (m/save.every)/2


### main MCMC run ###

run.MCMC <- function(){
  # s = 0.01 -> 19.3 %
  # s = 0.009 -> 23.4 %
  s = 0.0088

  system.time(
    d3 <- MCMC.MARS.CN(x0, xline.idx, inline.idx, s, m, save.every, IOmega, q3, h.hat.R0 = mars.R0.bf, h.hat.G = mars.G.bf)
  )

  # save x.0.posterior
  #l <- m / save.every
  #d.mean.gas <- apply(d3$gas[b:l,], 2, mean)
  #d.mean.oil <- apply(d3$oil[b:l,], 2, mean)
  #d.mean.clay <- apply(d3$clay[b:l,], 2, mean)
  #x0.posterior <- cbind(d.mean.gas, d.mean.oil, d.mean.clay)
  #save(x0.posterior, file="x0_posterior.Rda")



  print(paste("s = ", s))
  print(paste("acceptance rate = ",  d3$accepted/m * 100))

  plot.Sg(d3, burn.in = b, filename="MCMC_runs/plots/q3_0088_gas.pdf") 
  plot.So(d3, burn.in = b, filename="MCMC_runs/plots/q3_0088_oil.pdf")
  plot.Sb(d3, burn.in = b, filename="MCMC_runs/plots/q3_0088_brine.pdf") 
  plot.Vcl(d3, burn.in = b, filename="MCMC_runs/plots/q3_0088_clay.pdf") 

  filename <- "MCMC_runs/plots/q3_uncertainty_0088"
  percentile(d=d3, inline.idx, xline.idx, bi=b, p=80, filename=filename)

  pdf("MCMC_runs/plots/q3_trace_gas50.pdf", width=6, height=6)
  par(mar=c(5,5,4,2))
  plot(d3$gas[,50], type="l", xlab="", ylab="", cex.axis =1.3, cex.lab=1.4)
  dev.off()

  #save(d3, file="MCMC_runs/total_area/CN_MARS_015.Rda")


}

run.MCMC()
