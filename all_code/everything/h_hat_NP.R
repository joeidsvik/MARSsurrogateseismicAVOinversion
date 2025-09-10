

# File with NP models


library(np)
library(stats)
library(plotly)

# Data 
load(file="/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/df_train.Rda")
load(file="/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/df_test.Rda")


# models
load("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/NPR_R0_1000.Rda")
load("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/NPR_G_1000.Rda")

load("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/NPR_R0_4000.Rda")
load("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/NPR_G_4000.Rda")



fit.NPR.models <- function(){
  npr.R0.1000 <- fit.NP.R0(1000)
  npr.R0.4000 <- fit.NP.R0(4000)
  npr.G.1000 <- fit.NP.G(1000)
  npr.G.4000 <- fit.NP.G(4000)
  
  save(npr.R0.1000, file="h_hat_models/NPR_R0_1000.Rda")
  save(npr.R0.4000, file="h_hat_models/NPR_R0_4000.Rda")
  save(npr.G.1000, file="h_hat_models/NPR_G_1000.Rda")
  save(npr.G.4000, file="h_hat_models/NPR_G_4000.Rda")
}

fit.NP.2D <- function(df){
  #df must contain columns x1, x2 and y
  npr.model <- npreg(y ~ x1 + x2,
                     gradients = TRUE,
                     data = df)
  return(npr.model)
}


fit.NP.R0 <- function(n){
  npr.model <- npreg(R0 ~ x.gas + x.oil + x.clay + x.depth,
                     gradients = TRUE,
                     data = df.train[1:n,])
  return(npr.model)
}


fit.NP.G <- function(n){
  npr.model <- npreg(G ~ x.gas + x.oil + x.clay + x.depth,
                    gradients = TRUE,
                    data = df.train[1:n,])
  
  return(npr.model)
}


min.distance <- function(x, npr.model){
  dist <- distances(x, npr.model)
  row.idx <- which(dist==min(dist))
  return(row.idx)
}



distances <- function(x, datapoints){
  distances <- abs(x[,1] - datapoints[,1]) + abs(x[,2] - datapoints[,2]) + abs(x[,3] - datapoints[,3]) + abs(x[,4] - datapoints[,4])
  return(cbind(distances, datapoints[,6], row.names = NULL))
}

distances.2D <- function(x, datapoints){
  distances <- abs(x[,1] - datapoints[,1]) + abs(x[,2] - datapoints[,2]) 
  return(cbind(distances, datapoints[,4], row.names = NULL))
}




binary.search.min.distance <- function(x, npr.model){
  dist <- distances(x, npr.model)
  min.dist <-  min(dist)
  return(row.idx)
}



library(evmix)
library(mvtnorm)



getKh <- function(x.new, x, h, n){
  # Kg = kdgaussian(x = (x.new[1] - x[,1]), bw = h[1])
  # Ko = kdgaussian(x = (x.new[2] - x[,2]), bw = h[2])
  # Kc = kdgaussian(x = (x.new[3] - x[,3]), bw = h[3])
  # Kd = kdgaussian(x = (x.new[4] - x[,4]), bw = h[4])
  
  Kg = 1 / (4*pi^2*h[1]) * exp(-(x.new[1] - x[,1])^2/(2*h[1]^2))
  Ko = 1 / (4*pi^2*h[2]) * exp(-(x.new[2] - x[,2])^2/(2*h[2]^2))
  Kc = 1 / (4*pi^2*h[3]) * exp(-(x.new[3] - x[,3])^2/(2*h[3]^2))
  Kd = 1 / (4*pi^2*h[4]) * exp(-(x.new[4] - x[,4])^2/(2*h[4]^2))
  
  return(Kg*Ko*Kc*Kd)
}

getKh.2D <- function(x.new, x, h, n){
  
  Kx = 1 / (4*pi^2*h[1]) * exp(-(x.new[1] - x[,1])^2/(2*h[1]^2))
  Ky = 1 / (4*pi^2*h[2]) * exp(-(x.new[2] - x[,2])^2/(2*h[2]^2))

  return(Kx*Ky)
}



derivative.npreg.anal <- function(x.new, npr.model){
  x <- npr.model$eval
  n <- length(npr.model$eval[,1])
  h <- npr.model$bw
  
  D <- apply(X=x.new, MARGIN = 1, function(x.new.row){
    
    Kh <- getKh(as.numeric(x.new.row), x, h, n)
    
    dKg = (rep(x.new.row[1], n) - x[,1])/(h[1]^2)
    dKo = (rep(x.new.row[2], n) - x[,2])/(h[2]^2)
    dKc = (rep(x.new.row[3], n) - x[,3])/(h[3]^2)
    
    
    Y <- npr.model$mean
    Dg <- sum( (Y * Kh * (sum(dKg*Kh) - dKg*sum(Kh))) / sum(Kh)^2 )
    Do <- sum( (Y * Kh * (sum(dKo*Kh) - dKo*sum(Kh))) / sum(Kh)^2 )
    Dc <- sum( (Y * Kh * (sum(dKc*Kh) - dKc*sum(Kh))) / sum(Kh)^2 )
    return(c(Dg, Do, Dc))
  })
  return(D)
}


derivative.npreg.anal.2D <- function(x.new, npr.model){
  x <- npr.model$eval
  n <- length(npr.model$eval[,1])
  h <- npr.model$bw
  
  D <- apply(X=x.new, MARGIN = 1, function(x.new.row){
    
    Kh <- getKh.2D(as.numeric(x.new.row), x, h, n)
    
    dKx = (rep(x.new.row[1], n) - x[,1])/(h[1]^2)
    dKy = (rep(x.new.row[2], n) - x[,2])/(h[2]^2)
    
    
    Y <- npr.model$mean
    Dx <- sum( (Y * Kh * (sum(dKx*Kh) - dKx*sum(Kh))) / sum(Kh)^2 )
    Dy <- sum( (Y * Kh * (sum(dKy*Kh) - dKy*sum(Kh))) / sum(Kh)^2 )
    
    return(c(Dx, Dy))
  })
  return(D)
}



derivative.npreg <- function(x.new, X.cl, s, km.classif, npr.model){
  # X.cl : cluster of trainig data
  # s : neighbors to use in means of neighbors 
  # km.classif : kmeans obj
  
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x.new)), function(i) apply(km.classif$centers, 1, function(v) sum((x.new[i, ]-v)^2)))
  
  # find index of min distance (find cluster)
  cl <- max.col(-t(tmp))  
  
  # find min distance from x to points in that cluster
  d <- sapply(X=c(1:nrow(x.new)), function(i) {distances(x.new[i,], X.cl[[cl[i]]])})
  
  # find s smallest distances
  I <- sapply(X=d, function(D){D[order(D[,1]),][1:s, 2]})
  
  
  # find mean of gradient at the s closest data points
  G <- apply(X=I, MAR=2, function(i){gradients(npr.model)[i,]})
  
  
  gr <- apply(X=G, MAR=2, function(g){c(mean(g[1:s]), mean(g[(s+1):(2*s)]), mean(g[(2*s+1):(3*s)]))})

  return(gr)
}



derivative.npreg.2D.one.data.point <- function(x.new, X.cl, s, km.classif, npr.model){
  # X.cl : cluster of trainig data
  # s : neighbors to use in means of neighbors 
  # km.classif : kmeans obj
  
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- apply(km.classif$centers, 1, function(v) sum((x.new-v)^2))
  
  # find index of min distance (find cluster)
  cl <- max.col(-t(tmp))  
  
  # find min distance from x to points in that cluster
  d <- distances.2D(x.new, as.data.frame(X.cl[cl]))
  
  # find s smallest distances
  I <- d[order(d[,1]),][1:s, 2]
  
  
  # find mean of gradient at the s closest data points
  G <- gradients(npr.model)[I,]
  
  
  gr <- c(mean(G[,1]), mean(G[,2]))
  
  return(gr)
}



clusterX <- function(npr.model, km.classif){
  X <- cbind(npr.model$eval, km.classif$cluster,  c(1:nrow(npr.model$eval)))
  colnames(X) <- c("x.gas","x.oil","x.clay", "x.depth", "cl", "idx")
  
  X.cl = list()
  for (k in sort(unique(km.classif$cluster))){
    X.cl[[k]] <-  X[X$cl==k,]
  }
  
  return(X.cl)
}


clusterX.2D <- function(npr.model, km.classif){
  X <- cbind(npr.model$eval, km.classif$cluster,  c(1:nrow(npr.model$eval)))
  colnames(X) <- c("x","y", "cl", "idx")
  
  X.cl = list()
  for (k in sort(unique(km.classif$cluster))){
    X.cl[[k]] <-  X[X$cl==k,]
  }
  
  return(X.cl)
}


Kmeans.visualize <- function(npr.model){
  
  X <- npr.model$eval
  km.classif <- kmeans(x=X, centers=3)
  
  
  # Plot kmeans

  col <- c("dodgerblue3", "goldenrod1", "mediumslateblue", "forestgreen", 
           "brown2", "palevioletred1", "lightskyblue", "maroon", "sienna2",  "violet" )
  
  
  
  fig1 <- plot_ly(x = X[,1], y = X[,2], z = X[,3], color = km.classif$cluster, colors=col)
  fig1 <- fig1 %>% add_markers()
  fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = 'x.gas'),
                                       yaxis = list(title = 'x.oil'),
                                       zaxis = list(title = 'x.clay')))
  fig1
  
  
  fig2 <- plot_ly(x = X[,4], y = X[,2], z = X[,3], color = km.classif$cluster, colors=col)
  fig2 <- fig2 %>% add_markers()
  fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = 'x.depth'),
                                       yaxis = list(title = 'x.oil'),
                                       zaxis = list(title = 'x.clay')))
  fig2
  
  
  fig3 <- plot_ly(x = X[,4], y = X[,1], z = X[,3], color = km.classif$cluster, colors=col)
  fig3 <- fig3 %>% add_markers()
  fig3 <- fig3 %>% layout(scene = list(xaxis = list(title = 'x.depth'),
                                       yaxis = list(title = 'x.gas'),
                                       zaxis = list(title = 'x.clay')))
  fig3
  
  
  fig4 <- plot_ly(x = X[,4], y = X[,2], z = X[,1], color = km.classif$cluster, colors=col)
  fig4 <- fig4 %>% add_markers()
  fig4 <- fig4 %>% layout(scene = list(xaxis = list(title = 'x.depth'),
                                       yaxis = list(title = 'x.oil'),
                                       zaxis = list(title = 'x.gas')))
  fig4
}
