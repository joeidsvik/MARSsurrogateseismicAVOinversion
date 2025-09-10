### Alvheim
# Karen

source("main_fm.R")

# Likelihood ----
create.Omega0 <- function(){
  V.R0 <- 0.003
  V.G <- 0.03
  Cov.R0.G <- -0.6*sqrt(V.G*V.R0)
  Omega0 <- matrix(c(V.R0, Cov.R0.G, Cov.R0.G, V.G), nrow = 2)
  return(Omega0)
}


create.Omega <- function(Omega0, N){
  Omega <- bdiag(replicate(N,Omega0,simplify=FALSE))
  return(Omega)
}


getR0G <- function(inline.idx, xline.idx){
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


log.likelihood <- function(R0.G, x.prop, tt_area, IOmega){
  # R0.G skal være på samme form som simdata. Ønsker at getRG skal returnere dataen vi har i 2N vektor.
  mu.likelihood <- sim_forward_model(x.prop[,1], x.prop[,2], x.prop[,3], tt_area)
  p <- -0.5 * t(mu.likelihood - R0.G) %*% IOmega %*% (mu.likelihood-R0.G)
  return(p)
}


log.likelihood.well <- function(x.prop, inline.idx, xline.idx){
  #inline.idx = c(171:200)
  #xline.idx = c(86:115) 
  # x[,2] = x.oil
  # well <- c(1026, 4924), inline.idx=180:181, xline.idx=105
  #inline.idx = c(1:248)
  #xline.idx = c(1:178) 
  
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
  p <- - 1/(2*sd^2) * sum((mu.oil - x.prop.oil)^2) - 1/(2*sd^2) * sum((mu.gas - x.prop.gas)^2)
  return(p)
}
