### Alvheim
# Karen

rm(list = ls())

source("Alvheim_prior.R")
source("Alvheim_likelihood.R")
source("Alvheim_proposal.R")
source("Alvheim_transformations.R")
source("Alvheim_visualize.R")

library(Matrix)
library(nlme)
library(mcmcse)




# Markov chain Monte Carlo ----
MCMC <- function(x0, inline.idx, xline.idx, h, m, save.every, propose.q, log.dens.q){
  set.seed(2)
  # Number of grid points
  N <- length(inline.idx) * length(xline.idx)
  
  
  # Depth of selected area
  tt_area <- traveltimes[xline.idx, inline.idx]
  
  
  # If the well is not part of the selected area
  well.inline.idx <- c(180:181)
  well.xline.idx <- c(105)
  if (!(all(well.inline.idx %in% inline.idx, TRUE)) && !(well.xline.idx %in% xline.idx)){
    log.likelihood.well <- function(x.prop, inline.idx, xline.idx){
      return(0)
    }
  }
  
  # Data R0 G data with length 2N
  R0.G <- getR0G(inline.idx, xline.idx)
  
  
  # Expected value of prior
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  
  
  # Covariance matrix for prior and proposal
  Sig.g <- create.Sigma.g(inline.idx, xline.idx)
  ISig.g <- solve(Sig.g)
  
  Sig.o <- create.Sigma.o(inline.idx , xline.idx)
  ISig.o <- solve(Sig.o)
  
  Sig.cl <- create.Sigma.cl(inline.idx, xline.idx)
  ISig.cl <- solve(Sig.cl)
  
  
  # Covariance matrix for likelihood
  Omega0 <- create.Omega0()
  Omega <- create.Omega(Omega0, N = length(inline.idx)*length(xline.idx))
  IOmega <- solve(Omega)
  
  l0 <- as.numeric(log.likelihood(R0.G, x0, tt_area, IOmega) 
                   + log.likelihood.well(x0, inline.idx, xline.idx)
                   + log.prior(x0, mu.prior, ISig.g, ISig.o, ISig.cl))
      
  j <- m/save.every + 1
  x.data.gas <- data.frame(matrix(nrow = j, ncol = N)) 
  x.data.gas[1,] <- x0[,1]
  
  x.data.oil <- data.frame(matrix(nrow = j, ncol = N)) 
  x.data.oil[1,] <- x0[,2] 
  
  x.data.clay <- data.frame(matrix(nrow = j, ncol = N)) 
  x.data.clay[1,] <- x0[,3]
  
  l.data <- data.frame(posterior = rep(NaN, j))
  l.data[1,] <- l0
  
  accepted <- data.frame(accepted = rep(NaN, j))
  accepted[1,] <- TRUE
  
  x.i <- x0
  l.i <- l0
  
  idx <- 2
  #accepted <- 0
  for(i in 2:m){
    x.prop <- propose.q(x.i, mu.prior, h, N, Sig.g, Sig.o, Sig.cl)
    
    l.prop <- as.numeric(log.likelihood(R0.G, x.prop, tt_area, IOmega)
                         + log.likelihood.well(x.prop, inline.idx, xline.idx)
                         + log.prior(x.prop, mu.prior, ISig.g, ISig.o, ISig.cl))
    
    q.i <- log.dens.q(x.i, x.prop, h, ISig.g, ISig.o, ISig.cl, mu.prior)
    q.prop <- log.dens.q(x.prop, x.i, h, ISig.g, ISig.o, ISig.cl, mu.prior)

    alpha <- min(1, exp((l.prop + q.i - l.i - q.prop)) )
    
    u <- runif(1)
    if (alpha > u){
      #accepted <- accepted + 1
      x.i <- x.prop
      l.i <- l.prop
      
      if (!(i %% save.every)){
        x.data.gas[idx,] <- x.prop[,1]
        x.data.oil[idx,] <- x.prop[,2]
        x.data.clay[idx,] <- x.prop[,3]
        l.data[idx,] <- l.prop
        
        accepted[idx,] <- TRUE 
        idx <- idx + 1 
      }

    }
    else if (!(i %% save.every)){
        x.data.gas[idx,] <- x.i[,1]
        x.data.oil[idx,] <- x.i[,2]
        x.data.clay[idx,] <- x.i[,3]
        l.data[idx,] <- l.i
        accepted[idx,] <- FALSE
        idx <- idx + 1
    }
    if(!(i %% 10000)){
      print(i)
    }
  }
  #accepted <- accepted/m
  data <- list("gas" = x.data.gas, "oil" = x.data.oil, "clay" = x.data.clay, 
               "log.posterior" = l.data, "accepted" = accepted)
  return(data)
}




# Prepare ----
#inline.idx = c(1:248)
#xline.idx = c(1:178) 

# 30x30 area of interest 
inmin <- 171 - 15
inmax <- 200 + 15
xmin <- 86 - 15
xmax <- 115 + 15
by <- 2

inline.idx = seq(inmin, inmax, 2)
xline.idx = seq(xmin, xmax, 2)


unique(xline[xline.idx, inline.idx])

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/zeroOffset.pdf", width=6, height=6)
plot.R0(inline.idx, xline.idx)
#dev.off()

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/gradient.pdf", width=6, height=6)
plot.G(inline.idx, xline.idx)
#dev.off()

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/depthArea.pdf", width=6, height=6)
plot.traveltimes(inline.idx, xline.idx)
#dev.off()

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/depth.pdf", width=6, height=6)
plot.depth.with.area(inline.idx, xline.idx)
#dev.off()

N <- length(inline.idx) * length(xline.idx)

set.seed(2)
# rep(0, N)

tt_area <- traveltimes[xline.idx, inline.idx] # Depth of selected area
mu.prior <- match_traveltimes_mean_trend(tt_area) # Expected value of prior

  
x0.g <- mu.prior[,1] 
x0.o <- mu.prior[,2]   
x0.cl <- mu.prior[,3] 
x0 <- cbind(x0.g, x0.o, x0.cl)


# Save a lot of data frames ----

# q1
m <- 300
save.every <- 3
m/save.every
# Dq1 <- list
system.time(
  d005 <- MCMC(x0, inline.idx, xline.idx, h = 0.005, m, save.every, propose.q3, log.dens.q1)
)



d006 <- MCMC(x0, inline.idx, xline.idx, h = 0.006, m, save.every, propose.q1, log.dens.q1)
Dq1$d006 <- d006
result(Dq1, burn.in = 500)
save(Dq1, file = "Dq1.rda")

d007 <- MCMC(x0, inline.idx, xline.idx, h = 0.007, m, save.every, propose.q1, log.dens.q1)
Dq1$d007 <- d007
result(Dq1, burn.in = 500)
save(Dq1, file = "Dq1.rda")

d008 <- MCMC(x0, inline.idx, xline.idx, h = 0.008, m, save.every, propose.q1, log.dens.q1)
Dq1$d008 <- d008
result(Dq1, burn.in = 500)
save(Dq1, file = "Dq1.rda")


save(Dq1, file = "Dq1.rds")

load("Data/Dq1.rda")
result(Dq1, burn.in = 500)
rm(Dq1)


# q2
m <- 100000
save.every <- 100
m/save.every

#Dq2 <- list()

#d02 <- MCMC(x0, inline.idx, xline.idx, h = 0.02, m, save.every, propose.q2, log.dens.q2)
#Dq2$d02 <- d02
#result(Dq2, burn.in = 500)
#save(Dq2, file = "Dq2.rda")

#d025 <- MCMC(x0, inline.idx, xline.idx, h = 0.025, m, save.every, propose.q2, log.dens.q2)
#Dq2$d025 <- d025
#result(Dq2, burn.in = 500)
#save(Dq2, file = "Dq2.rda")

#d035 <- MCMC(x0, inline.idx, xline.idx, h = 0.035, m, save.every, propose.q2, log.dens.q2)
#Dq2$d035 <- d035
#result(Dq2, burn.in = 500)
#save(Dq2, file = "Dq2.rda")

#d03 <- MCMC(x0, inline.idx, xline.idx, h = 0.03, m, save.every, propose.q2, log.dens.q2)
#Dq2$d03 <- d03
#result(Dq2, burn.in = 500)
#save(Dq2, file = "Dq2.rda")

#d04 <- MCMC(x0, inline.idx, xline.idx, h = 0.04, m, save.every, propose.q2, log.dens.q2)
#Dq2$d04 <- d04
#result(Dq2, burn.in = 500)
#save(Dq2, file = "Dq2.rda")

#d045 <- MCMC(x0, inline.idx, xline.idx, h = 0.045, m, save.every, propose.q2, log.dens.q2)
#Dq2$d045 <- d045
#result(Dq2, burn.in = 500)
#save(Dq2, file = "Dq2.rda")

#load("Dq2.rda")
#result(Dq2, burn.in = 500)
#rm(Dq2)



# q3 
m <- 100000
save.every <- 100
m/save.every

#Dq3 <- list()

#d01 <- MCMC(x0, inline.idx, xline.idx, h = 0.01, m, save.every, propose.q3, log.dens.q3)
#Dq3$d01 <- d01
#result(Dq3, burn.in = 500)
#save(Dq3, file = "Dq3.rda")

#d02 <- MCMC(x0, inline.idx, xline.idx, h = 0.02, m, save.every, propose.q3, log.dens.q3)
#Dq3$d02 <- d02
#result(Dq3, burn.in = 500)
#save(Dq3, file = "Dq3.rda")

#d03 <- MCMC(x0, inline.idx, xline.idx, h = 0.03, m, save.every, propose.q3, log.dens.q3)
#Dq3$d03 <- d03
#result(Dq3, burn.in = 500)
#save(Dq3, file = "Dq3.rda")

#d04 <- MCMC(x0, inline.idx, xline.idx, h = 0.04, m, save.every, propose.q3, log.dens.q3)
#Dq3$d04 <- d04
#result(Dq3, burn.in = 500)
#save(Dq3, file = "Dq3.rda")

#d05 <- MCMC(x0, inline.idx, xline.idx, h = 0.05, m, save.every, propose.q3, log.dens.q3)
#Dq3$d05 <- d05
#result(Dq3, burn.in = 500)
#save(Dq3, file = "Dq3.rda")

#load("Dq3.rda")
#result(Dq3, burn.in = 500)
#save(Dq3, file = "Dq3.rda")
#rm(Dq3)

# Computation time ----

ptm <- proc.time()
MCMC(x0, inline.idx, xline.idx, h = , m = 3000000, save.every = 3000000, propose.q1, log.dens.q1)
comp.time1 <- proc.time() - ptm

ptm <- proc.time()
MCMC(x0, inline.idx, xline.idx, h = 0.035, m = 100000, save.every = 100000, propose.q2, log.dens.q2)
comp.time2 <- proc.time() - ptm

ptm <- proc.time()
MCMC(x0, inline.idx, xline.idx, h = 0.02, m = 100000, save.every = 100000, propose.q3, log.dens.q3)
comp.time3 <- proc.time() - ptm







