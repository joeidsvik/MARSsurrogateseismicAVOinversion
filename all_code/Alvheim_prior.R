

# File contains functions related to the likelihood in the Alvheim case


library(compiler)
source("functions_fft.R")


# Prior total area ----

# log prior total area
# used by q1, q2, q3 and q4
log.prior.total.area <- function(x, mu.prior, nxy, cm.g, cm.o, cm.c){
  p.g <- - pdf_fast_eval(x[,1], mu = mu.prior[,1], nxy = nxy, cm.g) # 0.5(x-mu)Sigma.inv(x-mu)
  p.o <- - pdf_fast_eval(x[,2], mu = mu.prior[,2], nxy = nxy, cm.o) # 0.5(x-mu)Sigma.inv(x-mu)
  p.c <- - pdf_fast_eval(x[,3], mu = mu.prior[,3], nxy = nxy, cm.c) # 0.5(x-mu)Sigma.inv(x-mu)
  p <- p.g + p.o + p.c
  return(p)
}


# Derivative of log prior
# used q4
derivative.log.prior.total.area <- function(x, mu.prior, nxy, cm.g, cm.o, cm.c){
  d.p.g <- deriv_fast_eval(x[,1], mu = mu.prior[,1], nxy = nxy, cm.g)
  d.p.o <- deriv_fast_eval(x[,2], mu = mu.prior[,2], nxy = nxy, cm.o)
  d.p.c <- deriv_fast_eval(x[,3], mu = mu.prior[,3], nxy = nxy, cm.c)
  d <- rbind(d.p.g, d.p.o, d.p.c)
  return(d)
}

log.prior.and.derivative.total.area <- function(x, mu.prior, nxy, cm.g, cm.o, cm.c){
  p.g <- - pdf_fast_eval(x[,1], mu = mu.prior[,1], nxy = nxy, cm.g) # 0.5(x-mu)Sigma.inv(x-mu)
  p.o <- - pdf_fast_eval(x[,2], mu = mu.prior[,2], nxy = nxy, cm.o) # 0.5(x-mu)Sigma.inv(x-mu)
  p.c <- - pdf_fast_eval(x[,3], mu = mu.prior[,3], nxy = nxy, cm.c) # 0.5(x-mu)Sigma.inv(x-mu)
  p <- p.g + p.o + p.c

  d.p.g <- deriv_fast_eval(x[,1], mu = mu.prior[,1], nxy = nxy, cm.g)
  d.p.o <- deriv_fast_eval(x[,2], mu = mu.prior[,2], nxy = nxy, cm.o)
  d.p.c <- deriv_fast_eval(x[,3], mu = mu.prior[,3], nxy = nxy, cm.c)
  d <- rbind(d.p.g, d.p.o, d.p.c)
  return(list("p" = p, "d" = d))
}


# Prior small area ----

# Log prior small area
# used by q1, q2, q3 and q4
log.prior <- function(x, mu.prior, ISig.g, ISig.o, ISig.c){
  p <- -0.5 * (t(x[,1]-mu.prior[,1]) %*% ISig.g  %*% (x[,1]-mu.prior[,1]) 
               + t(x[,2]-mu.prior[,2]) %*% ISig.o  %*% (x[,2]-mu.prior[,2]) 
               + t(x[,3]-mu.prior[,3]) %*% ISig.c  %*%(x[,3]-mu.prior[,3]))
  return(p)
}

# Log prior and derivative of log prior
# used q4
log.prior.and.derivative <- function(x, mu.prior, ISig.g, ISig.o, ISig.c){
  dxg <- c(ISig.g %*% (x[,1]-mu.prior[,1]))
  dxo <- c(ISig.o %*% (x[,2]-mu.prior[,2]))
  dxc <- c(ISig.c %*%(x[,3]-mu.prior[,3]))
  d <- - rbind(dxg, dxo, dxc)
  p <- - 0.5 * (t(x[,1]-mu.prior[,1]) %*% dxg 
                + t(x[,2]-mu.prior[,2]) %*% dxo 
                + t(x[,3]-mu.prior[,3]) %*% dxc )
  return(list("p" = p, "d" = d))
}


# Derivative of log prior
# used q4
derivative.log.prior <- function(x, mu.prior, ISig.g, ISig.o, ISig.c){
  dxg <- c(ISig.g %*% (x[,1]-mu.prior[,1]))
  dxo <- c(ISig.o %*% (x[,2]-mu.prior[,2]))
  dxc <- c(ISig.c %*%(x[,3]-mu.prior[,3]))
  d <- - rbind(dxg, dxo, dxc)
  return(d)
}


# aux functions for prior small area ----

# Calculate Manhattan distance between two points 
distance <- function(inline.idx1, xline.idx1, inline.idx2, xline.idx2){
  inline1 <- inline[xline.idx1, inline.idx1]
  xline1 <- xline[xline.idx1, inline.idx1]
  inline2 <- inline[xline.idx2, inline.idx2]
  xline2 <- xline[xline.idx2, inline.idx2]
  
  d <- sqrt((inline1 - inline2)^2 + (xline1 - xline2)^2)
  return(d)
}


i.to.inline.xline <- function(i, xline.idx, inline.idx){
  l.x <- length(xline.idx)
  #l.in <- length(inline.idx)

  e.i <- (i-1) %/% l.x # l.in # heltallsdivisjon
  n.i <- (i-1) %% l.x # mod
  idxs <- c(inline.idx[e.i+1], xline.idx[n.i+1])
  return(idxs)
}


create.Sigma.g <- function(xline.idx, inline.idx){
  sd.g <- 2.83
  l <- 15
  l.x <- length(xline.idx)
  l.in <- length(inline.idx)
  N <- l.x*l.in
  Sig <- matrix(data=NaN,nrow=N,ncol=N)
  for (r in 1:N){
    r.idx <- i.to.inline.xline(r, xline.idx, inline.idx)
    for (k in 1:N){
      k.idx <- i.to.inline.xline(k, xline.idx, inline.idx)
      d <- distance(r.idx[1], r.idx[2], k.idx[1], k.idx[2])
      Sig[r,k] <- sd.g^2 * exp(-0.5 * (d/l)^2)
    }
  }
  Sig = Sig + 0.01 * diag(N)
  return(Sig)
}


create.Sigma.o <- function(xline.idx, inline.idx){
  sd.o <- 2.68
  l <- 15
  l.x <- length(xline.idx)
  l.in <- length(inline.idx)
  N <- l.x*l.in
  Sig <- matrix(data=NaN,nrow=N,ncol=N)
  for (r in 1:N){
    r.idx <- i.to.inline.xline(r, xline.idx, inline.idx)
    for (k in 1:N){
      k.idx <- i.to.inline.xline(k, xline.idx, inline.idx)
      d <- distance(r.idx[1], r.idx[2], k.idx[1], k.idx[2])
      Sig[r,k] <- sd.o^2 * exp(-0.5 * (d/l)^2)
    }
  }
  Sig = Sig + 0.01 * diag(N)
  return(Sig)
}
create.Sigma.g <- cmpfun(create.Sigma.g)


create.Sigma.c <- function(xline.idx, inline.idx){
  sd.c <- 1.64
  l <- 15
  l.x <- length(xline.idx)
  l.in <- length(inline.idx)
  N <- l.x*l.in
  Sig <- matrix(data=NaN,nrow=N,ncol=N)
  for (r in 1:N){
    r.idx <- i.to.inline.xline(r, xline.idx, inline.idx)
    for (k in 1:N){
      k.idx <- i.to.inline.xline(k, xline.idx, inline.idx)
      d <- distance(r.idx[1], r.idx[2], k.idx[1], k.idx[2])
      Sig[r,k] <- sd.c^2 * exp(-0.5 * (d/l)^2)
    }
  }
  Sig = Sig + 0.01 * diag(N)
  return(Sig)
}


# Covariance matrix for prior and proposal
getSigmas <- function(xline.idx, inline.idx, LOAD=FALSE){
  if (LOAD){
    load("Sigma_obj/ISig_g.Rda")
    load("Sigma_obj/cholSig_g.Rda")
    
    load("Sigma_obj/ISig_o.Rda")
    load("Sigma_obj/cholSig_o.Rda")
    
    load("Sigma_obj/ISig_c.Rda")
    load("Sigma_obj/cholSig_c.Rda")
  } else{
    Sig.g <- create.Sigma.g(xline.idx, inline.idx)
    ISig.g <- solve(Sig.g)
    chol.Sig.g <- chol(Sig.g)
    save(ISig.g, file="Sigma_obj/ISig_g.Rda")
    save(chol.Sig.g, file="Sigma_obj/cholSig_g.Rda")
    
    Sig.o <- create.Sigma.o(xline.idx, inline.idx)
    ISig.o <- solve(Sig.o)
    chol.Sig.o <- chol(Sig.o)
    save(ISig.o, file="Sigma_obj/ISig_o.Rda")
    save(chol.Sig.o, file="Sigma_obj/cholSig_o.Rda")
    
    Sig.c <- create.Sigma.c(xline.idx, inline.idx)
    ISig.c <- solve(Sig.c)
    chol.Sig.c <- chol(Sig.c)
    save(ISig.c, file="Sigma_obj/ISig_c.Rda")
    save(chol.Sig.c, file="Sigma_obj/cholSig_c.Rda")
  }
  
  Sigmas=list("Ig" = ISig.g, "Io" = ISig.o, "Ic" = ISig.c, 
              "chol.g" = chol.Sig.g, "chol.o" = chol.Sig.o, "chol.c"= chol.Sig.c)
  return(Sigmas)
}



















