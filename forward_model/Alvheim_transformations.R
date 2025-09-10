### Alvheim
# Karen

# Transformations ----
S.g <- function(x){
  S <- exp(x[,1]) / (1 + exp(x[,1]) + exp(x[,2]))
  return(S)
}


S.o <- function(x){
  S <- exp(x[,2]) / (1 + exp(x[,1]) + exp(x[,2]))
  return(S)
}


S.b <- function(x){
  S <- 1 / (1 + exp(x[,1]) + exp(x[,2]))
  return(S)
}


V.cl <- function(x){
  V <- exp(x[,3]) / (1 + exp(x[,3]))
  return(V)
}