

# Contains
# create data frame drawn from prior
# calculate numerical derivative sim_forward_model



source("/Users/karenauestad/Dokumenter/master/Alvheim_kode/Alvheim_data.R")


# create data frame ----
create.df <- function(n){
  
  # Convert between travel times and depth
  convrt <- 1.055
  
  
  # sigma values
  sd.g <- 2.83
  sd.o <- 2.68
  sd.c <- 1.64
  
  
  # Lower and upper boundary for travel times
  min.traveltimes <- min(traveltimes)
  max.traveltimes <- max(traveltimes)
  
  
  # Empty data frame
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(df) <- c('x.gas', 'x.oil', 'x.clay', 'x.depth', 'R0', 'G')
  
  
  for (i in 1:n){
    # Draw travel time
    traveltime <- runif(1, min.traveltimes, max.traveltimes)
    
    # Depth
    depth <- traveltime*convrt
    
    # Prior mean
    mu.prior <- match_traveltimes_mean_trend(traveltime)
    
    # Draw x values
    x.g <- rnorm(1, mu.prior[1], sd.g)
    x.o <- rnorm(1, mu.prior[2], sd.o)
    x.c <- rnorm(1, mu.prior[3], sd.c)
    
    d <- length(x.g)
    
    # true gi
    sim_data <- sim_forward_model(x.g, x.o, x.c, traveltime)
    R0 <- sim_data[1]
    G <- sim_data[2]
    
    
    # Add row to data frame 
    
    df2 <- data.frame(matrix(c(x.g, x.o, x.c, depth, R0, G), ncol=6, nrow=1))
    colnames(df2) <- c('x.gas', 'x.oil', 'x.clay', 'x.depth', 'R0', 'G')
    df <- rbind(df,df2)
    
  }
  return(df)
}


# numerical derivative of sim_forward_model ----
sfm.derivative <- function(x){
  convrt <- 1.055
  x.g <- x[,1]
  x.o <- x[,2]
  x.c <- x[,3]
  tt <- x[,4]/convrt
  eps <- 0.00001
  d.smf.dxg <- (sim_forward_model(x.g + eps, x.o, x.c, tt) - sim_forward_model(x.g, x.o, x.c, tt)) / eps
  d.smf.dxo <- (sim_forward_model(x.g, x.o + eps, x.c, tt) - sim_forward_model(x.g, x.o, x.c, tt)) / eps 
  d.smf.dxc <- (sim_forward_model(x.g, x.o, x.c + eps, tt) - sim_forward_model(x.g, x.o, x.c, tt)) / eps
  dR0.SMF <- c(d.smf.dxg[1], d.smf.dxo[1], d.smf.dxc[1])
  dG.SMF <- c(d.smf.dxg[2], d.smf.dxo[2], d.smf.dxc[2])
  return(c(dR0.SMF, dG.SMF))
}



# performance metrics ----

# sample correlation
corr <- function(y, y.hat){
  mean.y <- mean(y)
  mean.y.hat <- mean(y.hat)
  
  # factor 1/(n-1) cancels everywhere
  cov.y.y.hat <- sum((y-mean.y)*(y.hat-mean.y.hat))
  sd.y <- sqrt(sum((y-mean.y)^2))
  sd.y.hat <- sqrt(sum((y.hat-mean.y.hat)^2)) 
  
  corr.y.y.hat <- cov.y.y.hat / (sd.y * sd.y.hat)
  return(corr.y.y.hat)
}



# Create and save training and test set ----

# df train
#df.train <- create.df(n=20000)
#save(df.train, file="h_hat_models/df_train.Rda")
#load("h_hat_models/df_train.Rda")


# df test
#df.test <- create.df(n=44144)
#save(df.test, file="h_hat_models/df_test.Rda")
#load("h_hat_models/df_test.Rda")


# create numerical derivatives for df_train and df_test
# library(tictoc)
# tic()
# n <- length(df.train[,1])
# columns = c("dR0dxg", "dR0dxo", "dR0dxc", "dGdxg", "dGdxo", "dGdxc")
# num.deriv.sfm.train = data.frame(matrix(nrow = n, ncol = length(columns)))
# colnames(num.deriv.sfm.train) = columns
# 
# for (i in 1:n){
#   x <- as.matrix(df.train[i,1:4])
# 
#   dsmf <- sfm.derivative(x)
# 
#   num.deriv.sfm.train$dR0dxg[i] <- dsmf[1]
#   num.deriv.sfm.train$dR0dxo[i] <- dsmf[2]
#   num.deriv.sfm.train$dR0dxc[i] <- dsmf[3]
# 
#   num.deriv.sfm.train$dGdxg[i] <- dsmf[4]
#   num.deriv.sfm.train$dGdxo[i] <- dsmf[5]
#   num.deriv.sfm.train$dGdxc[i] <- dsmf[6]
# }
# save(num.deriv.sfm.train, file="h_hat_models/derivative_sfm_df_train.Rda")
# toc()

# tic()
# n <- length(df.test[,1])
# columns = c("dR0dxg", "dR0dxo", "dR0dxc", "dGdxg", "dGdxo", "dGdxc")
# num.deriv.sfm = data.frame(matrix(nrow = n, ncol = length(columns)))
# colnames(num.deriv.sfm) = columns
# 
# for (i in 1:n){
#   x <- as.matrix(df.test[i,1:4])
# 
#   dsmf <- sfm.derivative(x)
# 
#   num.deriv.sfm$dR0dxg[i] <- dsmf[1]
#   num.deriv.sfm$dR0dxo[i] <- dsmf[2]
#   num.deriv.sfm$dR0dxc[i] <- dsmf[3]
# 
#   num.deriv.sfm$dGdxg[i] <- dsmf[4]
#   num.deriv.sfm$dGdxo[i] <- dsmf[5]
#   num.deriv.sfm$dGdxc[i] <- dsmf[6]
# }
# save(num.deriv.sfm, file="h_hat_models/derivative_sfm_df_test.Rda")
# toc()












