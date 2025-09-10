

# File contains functions related to the likelihood in the Alvheim case


# aux functions ----
# speed up calculation of : nabla h(x) (Omega (ys - h(x))^T
library(Rcpp)
cppFunction(
  "NumericVector organize_grad(int N, NumericVector holder, NumericVector dR0dx, NumericVector dGdx) {
  NumericVector dR0Gdx(N);

  for (int i=0; i < N+1; ++i) {
    dR0Gdx[i] = (holder[2*i] * dR0dx[i]) + (holder[2*i+1] * dGdx[i]);
  }
  return dR0Gdx;
}")


# Derivative of log likelihood avo with IOmega (y.s - h(x)) precomputed
derivative.log.likelihood.avo.precomputed.holder <- function(X, N, D.hat, holder){
  dR0dxg <- predict(object = D.hat[[1]], newdata = X)
  dGdxg <- predict(object = D.hat[[4]], newdata = X)
  dR0dxo <- predict(object = D.hat[[2]], newdata = X)
  dGdxo <- predict(object = D.hat[[5]], newdata = X)
  dR0dxc <- predict(object = D.hat[[3]], newdata = X)
  dGdxc <- predict(object = D.hat[[6]], newdata = X)
  
  dR0Gdxg = organize_grad(N, as.numeric(holder), dR0dxg, dGdxg)
  dR0Gdxo = organize_grad(N, as.numeric(holder), dR0dxo, dGdxo)
  dR0Gdxc = organize_grad(N, as.numeric(holder), dR0dxc, dGdxc)
  
  return(rbind(dR0Gdxg, dR0Gdxo, dR0Gdxc))
}



# merges two list alternating so that merge.lists.alternating(h.R0.list, h.G.list)
# gives the correct form [h.R0_1, h.G_1, h.R0_2, h.G_2, ...]
merge.lists.alternating <- function(list1, list2){
  #  length(list1) == length(list2)
  n <- length(list1) 
  L <- rep(0, 2*n)
  idx1 <- seq(1,2*n, 2)
  idx2 <- seq(2,2*n, 2)
  L[idx1] <- list1
  L[idx2] <- list2
  return(L)
}


# data y_s ----

# Get AVO data on the right form. 
getR0G <- function(xline.idx, inline.idx){
  R0 <- r0_data[xline.idx, inline.idx]
  G <- g_data[xline.idx, inline.idx]
  
  N <- length(inline.idx)*length(xline.idx)
  R0.G.list <- rep(NaN, 2*N)
  j <- 1
  for (i in seq(1,2*N,2)){
    R0.G.list[i] <- R0[j]
    R0.G.list[i+1] <- G[j]
    j <- j+1
  }
  return(R0.G.list)
}


# Likelihood ----

# Covariance matrix for likelihood 
create.Omega0 <- function(){
  V.R0 <- 0.003
  V.G <- 0.03
  Cov.R0.G <- -0.6*sqrt(V.G*V.R0)
  Omega0 <- matrix(c(V.R0, Cov.R0.G, Cov.R0.G, V.G), nrow = 2)
  return(Omega0)
}


# Change omega to try to make up for the error caused by using h_hat instead of h. 
create.Omega0.tilde <- function(){
  # from test set w 44144 data points 
  
  var.R0.hat <- 5.588405e-05
  var.G.hat <- 4.320285e-05
  cov.R0.G.hat <- -3.700215e-06
    
  V.R0 <- 0.003 + var.R0.hat
  V.G <- 0.03 + var.G.hat
  Cov.R0.G <- -0.6*sqrt(V.G*V.R0) + cov.R0.G.hat
  Omega0 <- matrix(c(V.R0, Cov.R0.G, Cov.R0.G, V.G), nrow = 2)
  return(Omega0)
}



# Covariance matrix for likelihood 
create.Omega <- function(Omega0, N){
  Omega <- bdiag(replicate(N,Omega0,simplify=FALSE))
  return(Omega)
}


# Change omega to try to make up for the error |h_hatt - h|
create.Omega0.tilde <- function(){
  # from test set w 44144 data points 
  # man(abs(h-h_hatt)) = 0.009, 0.007
  # man(abs(h-h_hatt))^2 = 0.009^2, 0.007^2
  # mse(h,h_hatt) = 0.00018, 0.00011
  diff.R0 <- 0.00018 
  diff.G <- 0.00011
    
  V.R0 <- 0.003 + diff.R0
  V.G <- 0.03 + diff.G
  Cov.R0.G <- -0.6*sqrt(V.G*V.R0)
  Omega0 <- matrix(c(V.R0, Cov.R0.G, Cov.R0.G, V.G), nrow = 2)
  return(Omega0)
}



# AVO sfm ----
# log likelihood with true forward model
log.likelihood.sfm <- function(R0.G, x, tt_area, IOmega){
  mu.likelihood <- sim_forward_model(x[,1], x[,2], x[,3], tt_area)
  p <- -0.5 * t(mu.likelihood - R0.G) %*% IOmega %*% (mu.likelihood-R0.G)
  return(p)
}



# AVO MARS ----
# log likelihood with approximated  forward model using MARS
# used by q1, q2, q3, q4
log.likelihood.avo.mars <- function(R0.G, x, depth_area, IOmega, h.hat.R0, h.hat.G, D.hat){
  X <- cbind(x,depth_area)
  
  mu.likelihood.R0 <- predict(object = h.hat.R0, newdata = X) 
  mu.likelihood.G <- predict(object = h.hat.G, newdata = X) 
  mu.likelihood <- merge.lists.alternating(mu.likelihood.R0, mu.likelihood.G)
  
  p <- -0.5 * t(R0.G-mu.likelihood) %*% IOmega %*% (R0.G-mu.likelihood)
  
  return(p)
}


# Log likelihood avo and derivative of log likelihood avo
# used by q4
log.likelihood.avo.and.derivative.mars <- function(R0.G, x, N, depth_area, IOmega, h.hat.R0, h.hat.G, D.hat){
  N <- length(depth_area)
  X <- cbind(x,depth_area)
  
  mu.likelihood.R0 <- predict(object = h.hat.R0, newdata = X) 
  mu.likelihood.G <- predict(object = h.hat.G, newdata = X) 
  mu.likelihood <- merge.lists.alternating(mu.likelihood.R0, mu.likelihood.G)
  
  holder <- IOmega %*% (R0.G-mu.likelihood)
  p <- - 0.5 * (R0.G-mu.likelihood) %*% holder
  
  dR0G <- derivative.log.likelihood.avo.precomputed.holder(X, N, D.hat, holder)

  return(list("p" = p, "d" = dR0G))
}


# Derivative of log likelihood avo 
# used by q4
derivative.log.likelihood.avo <- function(x, N, depth_area, h.hat.R0, h.hat.G, D.hat, IOmega, R0.G){
  X <- cbind(x,depth_area)
  
  mu.likelihood.R0 <- predict(object = h.hat.R0, newdata = X) 
  mu.likelihood.G <- predict(object = h.hat.G, newdata = X) 
  mu.likelihood <- merge.lists.alternating(mu.likelihood.R0, mu.likelihood.G)

  holder <- IOmega %*% (R0.G-mu.likelihood)
  
  #tic("predict")
  dR0dxg <- predict(object = D.hat[[1]], newdata = X)
  dGdxg <- predict(object = D.hat[[4]], newdata = X)
  dR0dxo <- predict(object = D.hat[[2]], newdata = X)
  dGdxo <- predict(object = D.hat[[5]], newdata = X)
  dR0dxc <- predict(object = D.hat[[3]], newdata = X)
  dGdxc <- predict(object = D.hat[[6]], newdata = X)
  #toc()

  dR0Gdxg = organize_grad(N, as.numeric(holder), dR0dxg, dGdxg)
  dR0Gdxo = organize_grad(N, as.numeric(holder), dR0dxo, dGdxo)
  dR0Gdxc = organize_grad(N, as.numeric(holder), dR0dxc, dGdxc)

  return(rbind(dR0Gdxg, dR0Gdxo, dR0Gdxc))
}


# Well data small area ----

# log likelihood of oil well in smaller are "area of interest"; inline = [156, 214], xline = [71, 129]
# used by q1, q2, q3 and q4
log.likelihood.well.4 <- function(x.prop, xline.idx, inline.idx, idx.well.4){
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  x.prop.oil <- x.prop[,2]
  x.prop.oil <- x.prop.oil[idx.well.4] # i oljeområdet 
  
  x.prop.gas <- x.prop[,1]
  x.prop.gas <- x.prop.gas[idx.well.4] # i oljeområdet 
  mu.oil <- rep(4, length(x.prop.oil)) # noe ekstremt
  mu.gas <- rep(-4, length(x.prop.oil))
  
  sd <- 0.1
  p <- - 1/(2*sd^2) * sum((x.prop.oil - mu.oil)^2) - 1/(2*sd^2) * sum((x.prop.gas - mu.gas)^2)

  return(p)
}


# Log likelihood well and derivative of log likelihood well
# used by q4
log.likelihood.well.4.and.derivative <- function(x, xline.idx, inline.idx, idx.well.4){
  x.oil <- x[,2]
  x.oil <- x.oil[idx.well.4]
  
  x.gas <- x[,1]
  x.gas <- x.gas[idx.well.4] 
  mu.oil <- rep(4, length(x.oil)) 
  mu.gas <- rep(-4, length(x.oil))
  
  sd <- 0.1
  p <- - 1/(2*sd^2) * sum((mu.oil - x.oil)^2) - 1/(2*sd^2) * sum((mu.gas - x.gas)^2)

  dxg <- rep(0, length(x[,1]))
  dxg[idx.well.4] <- - (1/sd^2) * (x.gas - mu.gas)
  dxo <- rep(0, length(x[,2]))
  dxo[idx.well.4] <- - (1/sd^2) * (x.oil - mu.oil) 
  dxc <- rep(0, length(x[,3]))

  d <- rbind(dxg, dxo, dxc)

  return(list("p" = p, "d" = d))

}

# Derivative of log likelihood well
# used by q4
derivative.log.likelihood.well.4 <- function(x.prop, inline.idx, xline.idx, idx.well.4){
  x.prop.oil <- x.prop[,2]
  x.prop.oil <- x.prop.oil[idx.well.4] 
  
  x.prop.gas <- x.prop[,1]
  x.prop.gas <- x.prop.gas[idx.well.4] 
  mu.oil <- rep(4, length(x.prop.oil)) 
  mu.gas <- rep(-4, length(x.prop.oil))
  
  sd <- 0.1
  dxg <- rep(0, length(x.prop[,1]))
  dxg[idx.well.4] <- - (1/sd^2) * (mu.gas - x.prop.gas)
  dxo <- rep(0, length(x.prop[,2]))
  dxo[idx.well.4] <- - (1/sd^2) * (mu.oil - x.prop.oil) 
  dxc <- rep(0, length(x.prop[,3]))

  return(rbind(dxg, dxo, dxc))
} 



# Well data total area ----
# log likelihood of all four wells
# used by q1, q2, q3 and q4
log.likelihood.wells <- function(x.prop, xline.idx, inline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4){
  x.prop.gas <- x.prop[,1]
  x.prop.oil <- x.prop[,2]
  
  
  # well 1 (gas) (484,4798)
  x.prop.gas.well.1 <- x.prop.gas[idx.well.1] # in well 1 area
  x.prop.oil.well.1 <- x.prop.oil[idx.well.1] # in well 1 area
  
  
  # well 2 (gas) (716,4914)
  x.prop.gas.well.2 <- x.prop.gas[idx.well.2] # in well 2 area
  x.prop.oil.well.2 <- x.prop.oil[idx.well.2] # in well 2 area
  
  
  # well 3 (gas) (1300,5150)
  x.prop.gas.well.3 <- x.prop.gas[idx.well.3] # in well 3 area
  x.prop.oil.well.3 <- x.prop.oil[idx.well.3] # in well 3 area
  
  
  # well 4 (oil) (1026,4924)
  x.prop.gas.well.4 <- x.prop.gas[idx.well.4] # in well 4 area
  x.prop.oil.well.4 <- x.prop.oil[idx.well.4] # in well 4 area
  
  # likelihood
  x.prop.gas.wells <- c(x.prop.gas.well.1, x.prop.gas.well.2, x.prop.gas.well.3, x.prop.gas.well.4)
  x.prop.oil.wells <- c(x.prop.oil.well.1, x.prop.oil.well.2, x.prop.oil.well.3, x.prop.oil.well.4)
  
  high <- 4
  low <- -4
  
  mu.gas <- c(rep(high, length(x.prop.gas.well.1)), rep(high, length(x.prop.gas.well.2)), 
               rep(high, length(x.prop.gas.well.3)), rep(low, length(x.prop.gas.well.4)))
  mu.oil <- c(rep(low, length(x.prop.oil.well.1)), rep(low, length(x.prop.oil.well.2)), 
               rep(low, length(x.prop.oil.well.3)), rep(high, length(x.prop.oil.well.4)) ) 
  
  sd <- 0.1
  p <- - 1/(2*sd^2) *( sum((mu.oil - x.prop.oil.wells)^2) + sum((mu.gas - x.prop.gas.wells)^2) )
  return(p)
}


# log likelihood and derivative of log likelihood
# used q4
log.likelihood.and.deriavtive.wells <- function(x.prop, xline.idx, inline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4){
  x.prop.gas <- x.prop[,1]
  x.prop.oil <- x.prop[,2]
  
  # well 1 (gas) (484,4798)
  x.prop.gas.well.1 <- x.prop.gas[idx.well.1] # in well 1 area
  x.prop.oil.well.1 <- x.prop.oil[idx.well.1] # in well 1 area
  
  
  # well 2 (gas) (716,4914)
  x.prop.gas.well.2 <- x.prop.gas[idx.well.2] # in well 2 area
  x.prop.oil.well.2 <- x.prop.oil[idx.well.2] # in well 2 area
  
  
  # well 3 (gas) (1300,5150)
  x.prop.gas.well.3 <- x.prop.gas[idx.well.3] # in well 3 area
  x.prop.oil.well.3 <- x.prop.oil[idx.well.3] # in well 3 area
  
  
  # well 4 (oil) (1026,4924)
  x.prop.gas.well.4 <- x.prop.gas[idx.well.4] # in well 4 area
  x.prop.oil.well.4 <- x.prop.oil[idx.well.4] # in well 4 area
  
  # likelihood
  x.prop.gas.wells <- c(x.prop.gas.well.1, x.prop.gas.well.2, x.prop.gas.well.3, x.prop.gas.well.4)
  x.prop.oil.wells <- c(x.prop.oil.well.1, x.prop.oil.well.2, x.prop.oil.well.3, x.prop.oil.well.4)
  
  high <- 4
  low <- -4
  
  mu.gas <- c(rep(high, length(x.prop.gas.well.1)), rep(high, length(x.prop.gas.well.2)), 
              rep(high, length(x.prop.gas.well.3)), rep(low, length(x.prop.gas.well.4)))
  mu.oil <- c(rep(low, length(x.prop.oil.well.1)), rep(low, length(x.prop.oil.well.2)), 
              rep(low, length(x.prop.oil.well.3)), rep(high, length(x.prop.oil.well.4)) ) 
  
  
  sd <- 0.1
  dxg <- rep(0, length(x.prop[,1]))
  dxo <- rep(0, length(x.prop[,2]))
  dxc <- rep(0, length(x.prop[,3]))
  
  dxg[idx.well.1] <- (1/sd^2) * (high - x.prop.gas.well.1) 
  dxo[idx.well.1] <- (1/sd^2) * (low - x.prop.oil.well.1)
  
  dxg[idx.well.2] <- (1/sd^2) * (high - x.prop.gas.well.2) 
  dxo[idx.well.2] <- (1/sd^2) * (low - x.prop.oil.well.2)
  
  dxg[idx.well.3] <- (1/sd^2) * (high - x.prop.gas.well.3) 
  dxo[idx.well.3] <- (1/sd^2) * (low - x.prop.oil.well.3)
  
  dxg[idx.well.4] <- (1/sd^2) * (low - x.prop.gas.well.4) 
  dxo[idx.well.4] <- (1/sd^2) * (high - x.prop.oil.well.4)

  
  p <- - 1/(2*sd^2) *( sum((mu.oil - x.prop.oil.wells)^2) + sum((mu.gas - x.prop.gas.wells)^2) )
  
  d <- rbind(dxg, dxo, dxc)
  
  return(list("p"=p, "d"=d))
}


# derivative of log likelihood
# used q4
derivative.log.likelihood.well <- function(x, inline.idx, xline.idx, idx.well.1, idx.well.2, idx.well.3, idx.well.4){
  x.prop.gas <- x[,1]
  x.prop.oil <- x[,2]
  
  # well 1 (gas) (484,4798)
  x.prop.gas.well.1 <- x.prop.gas[idx.well.1] # in well 1 area
  x.prop.oil.well.1 <- x.prop.oil[idx.well.1] # in well 1 area
  
  # well 2 (gas) (716,4914)
  x.prop.gas.well.2 <- x.prop.gas[idx.well.2] # in well 2 area
  x.prop.oil.well.2 <- x.prop.oil[idx.well.2] # in well 2 area
  
  # well 3 (gas) (1300,5150)
  x.prop.gas.well.3 <- x.prop.gas[idx.well.3] # in well 3 area
  x.prop.oil.well.3 <- x.prop.oil[idx.well.3] # in well 3 area
  
  # well 4 (oil) (1026,4924)
  x.prop.gas.well.4 <- x.prop.gas[idx.well.4] # in well 4 area
  x.prop.oil.well.4 <- x.prop.oil[idx.well.4] # in well 4 area

  # likelihood
  x.prop.gas.wells <- c(x.prop.gas.well.1, x.prop.gas.well.2, x.prop.gas.well.3, x.prop.gas.well.4)
  x.prop.oil.wells <- c(x.prop.oil.well.1, x.prop.oil.well.2, x.prop.oil.well.3, x.prop.oil.well.4)

  high <- 4
  low <- -4
  
  mu.gas <- c(rep(low, length(x.prop.gas.well.1)), rep(high, length(x.prop.gas.well.2)), 
              rep(high, length(x.prop.gas.well.3)), rep(low, length(x.prop.gas.well.4)))
  mu.oil <- c(rep(high, length(x.prop.oil.well.1)), rep(low, length(x.prop.oil.well.2)), 
              rep(low, length(x.prop.oil.well.3)), rep(high, length(x.prop.oil.well.4)) ) 
  
  sd <- 0.1
  dxg <- rep(0, length(x[,1]))
  dxo <- rep(0, length(x[,2]))
  dxc <- rep(0, length(x[,3]))
  
  dxg[idx.well.1] <- (1/sd^2) * (x.prop.gas.well.1 - 4) 
  dxo[idx.well.1] <- (1/sd^2) * (x.prop.oil.well.1 + 4)
  
  dxg[idx.well.2] <- (1/sd^2) * (x.prop.gas.well.2 - 4) 
  dxo[idx.well.2] <- (1/sd^2) * (x.prop.oil.well.2 + 4)
  
  dxg[idx.well.3] <- (1/sd^2) * (x.prop.gas.well.3 - 4) 
  dxo[idx.well.3] <- (1/sd^2) * (x.prop.oil.well.3 + 4)
  
  dxg[idx.well.4] <- (1/sd^2) * (x.prop.gas.well.4 + 4) 
  dxo[idx.well.4] <- (1/sd^2) * (x.prop.oil.well.4 - 4)
  
  
  d <- rbind(dxg, dxo, dxc)
  
  return(d)
} 

