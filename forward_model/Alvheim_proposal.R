### Alvheim
# Karen

# Proposal ----

# Normal random walk, q1
propose.q1 <- function(x.i, mu.prior, h, N, Sig.g, Sig.o, Sig.cl){
  N <- length(x.i[,1])
  x.prop <- cbind(x.i[,1] + h * rnorm(N), x.i[,2] + h * rnorm(N), x.i[,3] + h * rnorm(N))
  return(x.prop)
}


log.dens.q1 <- function(x.prop=0, x.i=0, h=0, ISig.g=0, ISig.o=0, ISig.cl=0, mu.prior=0){
  return(0) # symmetrisk
}


# Normal random walk with covariance, q2
propose.q2 <- function(x.i, mu.prior, h, N, Sig.g, Sig.o, Sig.cl){
  x.i.g <- x.i[,1]
  x.i.o <- x.i[,2]
  x.i.cl <- x.i[,3]
  
  x.prop.g <- c(rmvnorm(n = 1, mean = x.i.g, sigma = h^2 * Sig.g))
  x.prop.o <- c(rmvnorm(n = 1, mean = x.i.o, sigma = h^2 * Sig.o)) 
  x.prop.cl <- c(rmvnorm(n = 1, mean = x.i.cl, sigma = h^2 * Sig.cl))
  
  x.prop <- cbind(x.prop.g, x.prop.o, x.prop.cl)
  return(x.prop)
}


log.dens.q2 <- function(x.prop=0, x.i=0, h=0, ISig.g=0, ISig.o=0, ISig.cl=0, mu.prior=0){
  return(0) # symmetrisk
}


# q3
propose.q3 <- function(x.i, mu.prior, h, N, Sig.g, Sig.o, Sig.cl){
  mu.g <- sqrt(1-h^2)*x.i[,1] #- mu.prior[,1])
  mu.o <- sqrt(1-h^2)*x.i[,2] #- mu.prior[,2])
  mu.cl <- sqrt(1-h^2)*x.i[,3] #- mu.prior[,3])
  
  x.prop.g <- c(rmvnorm(n = 1, mean = mu.g, sigma = h^2 * Sig.g))
  x.prop.o <- c(rmvnorm(n = 1, mean = mu.o, sigma = h^2 * Sig.o)) 
  x.prop.cl <- c(rmvnorm(n = 1, mean = mu.cl, sigma = h^2 * Sig.cl))
  
  x.prop <- cbind(x.prop.g, x.prop.o, x.prop.cl)
  #print(x.prop.g[1])
  return(x.prop)
}


log.dens.q3 <- function(x.prop, x.i, h, ISig.g, ISig.o, ISig.cl, mu.prior){
  x.prop.g <- x.prop[,1]
  x.prop.o <- x.prop[,2]
  x.prop.cl <- x.prop[,3]
  
  mu.g <- sqrt(1-h^2)*x.i[,1] # - mu.prior[,1])
  mu.o <- sqrt(1-h^2)*x.i[,2] # - mu.prior[,2])
  mu.cl <- sqrt(1-h^2)*x.i[,3] # - mu.prior[,3])
  
  p <- -(0.5 / h^2) * ( t(x.prop.g-mu.g) %*% ISig.g %*% (x.prop.g-mu.g) 
                      + t(x.prop.o-mu.o) %*% ISig.o %*% (x.prop.o-mu.o)
                      + t(x.prop.cl-mu.cl) %*% ISig.cl %*% (x.prop.cl-mu.cl) )
  return(p)
}
