### Alvheim
# Karen

# Prior ----
# Calculate Manhattan distance between two points 
distance <- function(inline.idx1, xline.idx1, inline.idx2, xline.idx2){
  inline1 <- inline[xline.idx1, inline.idx1]
  xline1 <- xline[xline.idx1, inline.idx1]
  inline2 <- inline[xline.idx2, inline.idx2]
  xline2 <- xline[xline.idx2, inline.idx2]
  
  d <- sqrt((inline1 - inline2)^2 + (xline1 - xline2)^2)
  return(d)
}


i.to.inline.xline <- function(i, inline.idx, xline.idx){
  r <- length(xline.idx)
  k <- length(inline.idx)
  #e.i <- (i-1) %% k # mod
  #n.i <- (i-1) %/% r # heltallsdivisjon
  e.i <- (i-1) %/% k # heltallsdivisjon
  n.i <- (i-1) %% r # mod
  idxs <- c(inline.idx[e.i+1], xline.idx[n.i+1])
  return(idxs)
}


create.Sigma.g <- function(inline.idx, xline.idx){
  sd.g <- 2.83
  l <- 15
  N <- length(inline.idx) * length(xline.idx)
  Sig <- matrix(data=NaN,nrow=N,ncol=N)
  for (r in 1:N){
    r.idx <- i.to.inline.xline(r, inline.idx, xline.idx)
    for (k in 1:N){
      k.idx <- i.to.inline.xline(k, inline.idx, xline.idx)
      d <- distance(r.idx[1], r.idx[2], k.idx[1], k.idx[2])
      Sig[r,k] <- sd.g^2 * exp(-0.5 * (d/l)^2)
    }
  }
  Sig = Sig + 0.01 * diag(N)
  return(Sig)
}


create.Sigma.o <- function(inline.idx, xline.idx){
  sd.o <- 2.68
  l <- 15
  N <- length(inline.idx) * length(xline.idx)
  Sig <- matrix(data=NaN,nrow=N,ncol=N)
  for (r in 1:N){
    r.idx <- i.to.inline.xline(r, inline.idx, xline.idx)
    for (k in 1:N){
      k.idx <- i.to.inline.xline(k, inline.idx, xline.idx)
      d <- distance(r.idx[1], r.idx[2], k.idx[1], k.idx[2])
      Sig[r,k] <- sd.o^2 * exp(-0.5 * (d/l)^2)
    }
  }
  Sig = Sig + 0.01 * diag(N)
  return(Sig)
}


create.Sigma.cl <- function(inline.idx, xline.idx){
  sd.cl <- 1.64
  l <- 15
  N <- length(inline.idx) * length(xline.idx)
  Sig <- matrix(data=NaN,nrow=N,ncol=N)
  for (r in 1:N){
    r.idx <- i.to.inline.xline(r, inline.idx, xline.idx)
    for (k in 1:N){
      k.idx <- i.to.inline.xline(k, inline.idx, xline.idx)
      d <- distance(r.idx[1], r.idx[2], k.idx[1], k.idx[2])
      Sig[r,k] <- sd.cl^2 * exp(-0.5 * (d/l)^2)
    }
  }
  Sig = Sig + 0.01 * diag(N)
  return(Sig)
}


# Log prior
log.prior <- function(x.prop, mu.prior, ISig.g, ISig.o, ISig.cl){
  p <- -0.5 * ( t(x.prop[,1]-mu.prior[,1]) %*% ISig.g %*% (x.prop[,1]-mu.prior[,1]) 
                + t(x.prop[,2]-mu.prior[,2]) %*% ISig.o %*% (x.prop[,2]-mu.prior[,2])
                + t(x.prop[,3]-mu.prior[,3]) %*% ISig.cl %*% (x.prop[,3]-mu.prior[,3]) )
  return(p)
}
