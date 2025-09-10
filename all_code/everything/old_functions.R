

derivative.log.likelihood.well.4 <- function(x.prop, inline.idx, xline.idx){
  i <- which(xline.idx == 105)
  j <- which(inline.idx == 180)
  
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  x.prop.oil <- x.prop[,2]
  x.prop.oil <- x.prop.oil[c(((j-1)*r+i),((j)*r+i))] # i oljeområdet 
  
  x.prop.gas <- x.prop[,1]
  x.prop.gas <- x.prop.gas[c(((j-1)*r+i),((j)*r+i))] # i oljeområdet 
  mu.oil <- rep(4, length(x.prop.oil)) # noe ekstremt
  mu.gas <- rep(-4, length(x.prop.oil))
  
  sd <- 0.1
  dxg <- - (1/sd^2) * sum(mu.gas - x.prop.gas)
  dxo <- - (1/sd^2) * sum(mu.oil - x.prop.oil) 
  
  return(c(dxg, dxo))
} 


# Derivative of log likelihood avo 
derivative.log.likelihood.avo_old <- function(x, depth_area, h.hat.R0, h.hat.G, D.hat, IOmega, R0.G){
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
  
  dR0Gdxg <- rep(0, N)
  dR0Gdxo <- rep(0, N)
  dR0Gdxc <- rep(0, N)
  
  for (i in 1:N){
    dR0Gdxg[i] <- holder[2*i-1] * dR0dxg[i] + holder[2*i] * dGdxg[i]
    dR0Gdxo[i] <- holder[2*i-1] * dR0dxo[i] + holder[2*i] * dGdxo[i]
    dR0Gdxc[i] <- holder[2*i-1] * dR0dxc[i] + holder[2*i] * dGdxc[i]
  }
  
  return(rbind(dR0Gdxg, dR0Gdxo, dR0Gdxc))
}



# Derivative of log likelihood avo with IOmega (y.s - h(x)) precomputed
derivative.log.likelihood.avo.precomputed.holder_old <- function(X, N, D.hat, holder){
  dR0dxg <- predict(object = D.hat[[1]], newdata = X)
  dGdxg <- predict(object = D.hat[[4]], newdata = X)
  dR0dxo <- predict(object = D.hat[[2]], newdata = X)
  dGdxo <- predict(object = D.hat[[5]], newdata = X)
  dR0dxc <- predict(object = D.hat[[3]], newdata = X)
  dGdxc <- predict(object = D.hat[[6]], newdata = X)
  
  dR0Gdxg <- rep(0, N)
  dR0Gdxo <- rep(0, N)
  dR0Gdxc <- rep(0, N)
  
  j <- 1
  for (i in 1:N){
    dR0Gdxg[i] <- holder[j] * dR0dxg[i] + holder[j+1] * dGdxg[i]
    dR0Gdxo[i] <- holder[j] * dR0dxo[i] + holder[j+1] * dGdxo[i]
    dR0Gdxc[i] <- holder[j] * dR0dxc[i] + holder[j+1] * dGdxc[i]
    j <- j+2
  }
  
  return(rbind(dR0Gdxg, dR0Gdxo, dR0Gdxc))
}



derivative.log.likelihood.avo.precomputed.holder.old <- function(x, depth_area, holder, D.hat){
  X <- cbind(x,depth_area)

  dR0G_mat <- gradient.h.matrix(X, D.hat)

  log.derivative.likelihood.X <-  t(holder) %*% dR0G_mat 
  log.derivative.likelihood <- matrix(log.derivative.likelihood.X, nrow=3)

  return(log.derivative.likelihood)
}

# Log likelihood avo and derivative of log likelihood avo
# old
log.likelihood.avo.and.derivative.mars.old <- function(R0.G, x, depth_area, IOmega, h.hat.R0, h.hat.G, D.hat){
  N <- length(depth_area)
  X <- cbind(x,depth_area)
  
  mu.likelihood.R0 <- predict(object = h.hat.R0, newdata = X) 
  mu.likelihood.G <- predict(object = h.hat.G, newdata = X) 
  mu.likelihood <- merge.lists.alternating(mu.likelihood.R0, mu.likelihood.G)
  
  holder <- IOmega %*% (R0.G-mu.likelihood)
  p <- - 0.5 * (R0.G-mu.likelihood) %*% holder
  dR0G <- derivative.log.likelihood.avo.precomputed.holder(x, depth_area, holder, D.hat)

  return(list("p" = p, "d" = dR0G))
}



# Aux function to create gradient of h(x)
# old
gradient.h.matrix <- function(X, D.hat){
  dR0dxg <- predict(object = D.hat[[1]], newdata = X)
  dGdxg <- predict(object = D.hat[[4]], newdata = X)
  dR0dxo <- predict(object = D.hat[[2]], newdata = X)
  dGdxo <- predict(object = D.hat[[5]], newdata = X)
  dR0dxc <- predict(object = D.hat[[3]], newdata = X)
  dGdxgc <- predict(object = D.hat[[6]], newdata = X)

  dR0G_mat <- matrix(0, ncol=3*N, nrow=2*N)

  for (i in 1:N){
    row.idx <- (2*i-1):(2*i)
    col.idx <- (3*i-2):(3*i)
    dR0G_mat[row.idx,col.idx] <- matrix(rbind(c(dR0dxg[i], dR0dxo[i], dR0dxc[i]), c(dR0dxg[i], dR0dxo[i], dR0dxc[i])), ncol=3) 
  }

  return(dR0G_mat)
}