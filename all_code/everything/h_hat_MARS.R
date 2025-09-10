

# File with MARS models


library(earth)     
library(caret) 
library(fastmatch)

load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/h_hat_models/df_train.Rda")
load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/h_hat_models/df_test.Rda")


load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/h_hat_models/mars_G_bf.Rda")
load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/h_hat_models/mars_R0_bf.Rda")
load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/h_hat_models/D_mars.Rda")

load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/h_hat_models/mars_R0_of.Rda")
load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/h_hat_models/mars_G_of.Rda")


# fit models ----

fit.mars.2D <- function(df){
  #df must contain columns x1, x2 and y
  mars.model <- earth(
    x = subset(df, select = c(x1, x2)),
    y = df$y,
    penalty=3,
    degree=2,
    nk=50,
    pmethod="exhaustive"
  )
  return(mars.model)
}


fit.mars.2D.large <- function(df){
  mars <- earth(x = subset(df, select = c(x1, x2)),
                y = df$y,
                penalty = 0,
                pmethod="none",
                nk=200,
                thresh= 0.0000000001)
  return(mars)
}


# MARS.R0
fit.mars.R0.large <- function(data.set){
  mars.R0 <- earth(
    x = subset(data.set, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = data.set$R0,
    penalty = 0,
    pmethod="exhaustive",
    nk=80,
    thresh= 0.00001)
  return(mars.R0)
}


fit.mars.R0 <- function(data.set){
  mars.R0 <- earth(
    x = subset(data.set, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = data.set$R0,
    penalty=3,
    degree=4,
    nk=50,
    pmethod="exhaustive"
  )
  return(mars.R0)
}



# MARS.G
fit.mars.G.large <- function(data.set){
  mars.G <- earth(
    x = subset(data.set, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = data.set$G,
    penalty = 0,
    pmethod="exhaustive",
    nk=80,
    thresh= 0.00001)
  return(mars.G)
}


# MARS.G
fit.mars.G <- function(data.set){
  mars.G <- earth(
    x = subset(data.set, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = data.set$G,
    degree=4,
    penalty=3,
    nk=50,
    pmethod="exhaustive"
  )
  return(mars.G)
}



# create MARS models ----


create.mars.models <- function(){
  load("h_hat_models/df_train.Rda")
  load("h_hat_models/df_test.Rda")
  
  
  # mars.R0.bf_2 <- fit.mars.R0(df.train)
  # summary(mars.R0.bf_2)
  # system.time(mars.R0.pred <- predict(object= mars.R0.bf_2, newdata = df.test))
  # corr(df.test$R0, mars.R0.pred)
  # mean((df.test$R0 - mars.R0.pred)^2)
  # 
  # mars.G.bf_2 <- fit.mars.G(df.train)
  # summary(mars.G.bf_2)
  # system.time(mars.G.pred <- predict(object= mars.G.bf_2, newdata = df.test))
  # corr(df.test$G, mars.G.pred)
  # mean((df.test$G - mars.G.pred)^2)
  
  
  # best fit
  #mars.R0.bf <- fit.mars.R0(df.train)
  #save(mars.R0.bf, file="h_hat_models/mars_R0_bf.Rda")
  load("h_hat_models/mars_R0_bf.Rda")
  
  #mars.G.bf <- fit.mars.G(df.train)
  #save(mars.G.bf, file="h_hat_models/mars_G_bf.Rda")
  load("h_hat_models/mars_G_bf.Rda")
  
  # over fit
  #mars.R0.of <- fit.mars.R0.large(df.train)
  #save(mars.R0.of, file="h_hat_models/mars_R0_of.Rda")
  load("h_hat_models/mars_R0_of.Rda")
  
  #mars.G.of <- fit.mars.G.large(df.train)
  #save(mars.G.of, file="h_hat_models/mars_G_of.Rda")
  load("h_hat_models/mars_G_of.Rda")
  
}







# Derivatives ----

## Mars model for derivatives ----


fit.mars.derivative <- function(){
  
  load("h_hat_models/df_train.Rda")
  load("h_hat_models/derivative_sfm_df_train.Rda")
  df.train.w.d <- cbind(df.train, num.deriv.sfm.train)
  
  
  param_grid <- expand.grid(degree = seq(1,10,2), nprune = floor(seq(5, 50, length.out = 20)))
  
  
  # dR0dxg
  cv_marsdR0dxg <- train(
    x = subset(df.train.w.d, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = df.train.w.d$dR0dxg,
    method = "earth",
    penalty=0,
    metric = "RMSE",
    trControl = trainControl(method = "cv", number = 5),
    tuneGrid = param_grid)
  mars.dR0dxg <- cv_marsdR0dxg$finalModel
  
  
  # dR0dxo
  cv_marsdR0dxo <- train(
    x = subset(df.train.w.d, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = df.train.w.d$dR0dxo,
    method = "earth",
    penalty=0,
    metric = "RMSE",
    trControl = trainControl(method = "cv", number = 5),
    tuneGrid = param_grid)
  mars.dR0dxo <- cv_marsdR0dxo$finalModel
  
  
  # dR0dxc
  cv_marsdR0dxc <- train(
    x = subset(df.train.w.d, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = df.train.w.d$dR0dxc,
    method = "earth",
    penalty=0,
    metric = "RMSE",
    trControl = trainControl(method = "cv", number = 5),
    tuneGrid = param_grid)
  mars.dR0dxc <- cv_marsdR0dxc$finalModel
  
  
  # dGdxg
  cv_marsdGdxg <- train(
    x = subset(df.train.w.d, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = df.train.w.d$dGdxg,
    method = "earth",
    penalty=0,
    metric = "RMSE",
    trControl = trainControl(method = "cv", number = 5),
    tuneGrid = param_grid)
  mars.dGdxg <- cv_marsdGdxg$finalModel
  
  
  # dGdxo
  cv_marsdGdxo <- train(
    x = subset(df.train.w.d, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = df.train.w.d$dGdxo,
    method = "earth",
    penalty=0,
    metric = "RMSE",
    trControl = trainControl(method = "cv", number = 5),
    tuneGrid = param_grid)
  mars.dGdxo <- cv_marsdGdxo$finalModel
  
  
  # dGdxc
  cv_marsdGdxc <- train(
    x = subset(df.train.w.d, select = c(x.gas, x.oil, x.clay, x.depth)),
    y = df.train.w.d$dGdxc,
    method = "earth",
    penalty=0,
    metric = "RMSE",
    trControl = trainControl(method = "cv", number = 5),
    tuneGrid = param_grid)
  mars.dGdxc <- cv_marsdGdxc$finalModel
  
  return(list(mars.dR0dxg, mars.dR0dxo, mars.dR0dxc, mars.dGdxg, mars.dGdxo, mars.dGdxc))
}



mars.derivatives.pred.R0 <- function(D.mars, data){
  derivatives <- matrix(, nrow=dim(data)[1], ncol=3)
  derivatives[,1] <- predict(object = D.mars[[1]], newdata = data)
  derivatives[,2] <- predict(object = D.mars[[2]], newdata = data)
  derivatives[,3] <- predict(object = D.mars[[3]], newdata = data)
  return(derivatives)
}

mars.derivatives.pred.G <- function(D.mars, data){
  derivatives <- matrix(, nrow=dim(data)[1], ncol=3)
  derivatives[,1] <- predict(object = D.mars[[4]], newdata = data)
  derivatives[,2] <- predict(object = D.mars[[5]], newdata = data)
  derivatives[,3] <- predict(object = D.mars[[6]], newdata = data)
  return(derivatives)
}

mars.derivatives.pred <- function(D.mars, data){
  derivatives <- matrix(, nrow=dim(data)[1], ncol=6)
  derivatives[,1] <- predict(object = D.mars[[1]], newdata = data)
  derivatives[,2] <- predict(object = D.mars[[2]], newdata = data)
  derivatives[,3] <- predict(object = D.mars[[3]], newdata = data)
  derivatives[,4] <- predict(object = D.mars[[4]], newdata = data)
  derivatives[,5] <- predict(object = D.mars[[5]], newdata = data)
  derivatives[,6] <- predict(object = D.mars[[6]], newdata = data)
  return(derivatives)
}




## Aux functions ----
# Get cuts and dirs
getCutsDirs <- function(mars.model){
  mars.cuts <- as.matrix(mars.model$cuts[mars.model$selected.terms,])
  mars.dirs <- as.matrix(mars.model$dirs[mars.model$selected.terms,])
  return(list(cuts = mars.cuts, dirs = mars.dirs))
}



get.cu <- function(mars.model){
  mars.cuts <- as.matrix(mars.model$cuts[mars.model$selected.terms,])
  mars.dirs <- as.matrix(mars.model$dirs[mars.model$selected.terms,])
  
  M <-  length(mars.model$coefficients)
  
  # (Mx4) which variable has cuts at term m 
  non.zero <- t(apply(X=mars.dirs, MARGIN=1, function(md) {md != 0} ))
  
  idx.non.interactions <- which(rowSums(non.zero) == 1)
  idx.interactions <- which(rowSums(non.zero) > 1)
  
  # non interactions
  cu <- as.data.frame(mars.cuts[idx.non.interactions,])
  cu <- cbind(cu, idx.non.interactions)
  
  cu.gas <- cu[which(cu[,1] != 0), ]
  cu.gas.sorted <- cu.gas[order(cu.gas[,1]), ]
  cu.gas.sorted.pos <- cu.gas.sorted[which(mars.dirs[cu.gas.sorted[,5], 1]==1), ]
  cu.gas.neg <- cu.gas.sorted[which(mars.dirs[cu.gas.sorted[,5], 1]==-1), ]
  cu.gas.sorted.neg <- cu.gas.neg[order(cu.gas.neg[,3], decreasing = TRUE), ]
  
  cu.oil <- cu[which(cu[,2] != 0), ]
  cu.oil.sorted <- cu.oil[order(cu.oil[,2]), ]
  cu.oil.sorted.pos <- cu.oil.sorted[which(mars.dirs[cu.oil.sorted[,5], 2]==1), ]
  cu.oil.neg <- cu.oil.sorted[which(mars.dirs[cu.oil.sorted[,5], 2]==-1), ]
  cu.oil.sorted.neg <- cu.oil.neg[order(cu.oil.neg[,3], decreasing = TRUE), ]
  
  cu.clay <- cu[which(cu[,3] != 0), ]
  cu.clay.sorted <- cu.clay[order(cu.clay[,3]), ]
  cu.clay.sorted.pos <- cu.clay.sorted[which(mars.dirs[cu.clay.sorted[,5], 3]==1), ]
  cu.clay.neg <- cu.clay.sorted[which(mars.dirs[cu.clay.sorted[,5], 3]==-1), ]
  cu.clay.sorted.neg <- cu.clay.neg[order(cu.clay.neg[,3], decreasing = TRUE), ]
  
  # interactions
  cu.int <- as.data.frame(mars.cuts[idx.interactions,])
  cu.int <- cbind(cu.int, idx.interactions)
  
  cu.gas.int <- cu.int[which(cu.int[,1] != 0), ]
  cu.gas.sorted.int <- cu.gas.int[order(cu.gas.int[,1]), ]
  cu.gas.sorted.pos.int <- cu.gas.sorted.int[which(mars.dirs[cu.gas.sorted.int[,5], 1]==1), ]
  cu.gas.neg.int <- cu.gas.sorted.int[which(mars.dirs[cu.gas.sorted.int[,5], 1]==-1), ]
  cu.gas.sorted.neg.int <- cu.gas.neg.int[order(cu.gas.neg.int[,1], decreasing = TRUE), ] # 3
  
  cu.oil.int <- cu.int[which(cu.int[,2] != 0), ]
  cu.oil.sorted.int <- cu.oil.int[order(cu.oil.int[,2]), ]
  cu.oil.sorted.pos.int <- cu.oil.sorted.int[which(mars.dirs[cu.oil.sorted.int[,5], 2]==1), ]
  cu.oil.neg.int <- cu.oil.sorted.int[which(mars.dirs[cu.oil.sorted.int[,5], 2]==-1), ]
  cu.oil.sorted.neg.int <- cu.oil.neg.int[order(cu.oil.neg.int[,2], decreasing = TRUE), ] # 3
  
  cu.clay.int <- cu.int[which(cu.int[,3] != 0), ]
  cu.clay.sorted.int <- cu.clay.int[order(cu.clay.int[,3]), ]
  cu.clay.sorted.pos.int <- cu.clay.sorted.int[which(mars.dirs[cu.clay.sorted.int[,5], 3]==1), ]
  cu.clay.neg.int <- cu.clay.sorted.int[which(mars.dirs[cu.clay.sorted.int[,5], 3]==-1), ]
  cu.clay.sorted.neg.int <- cu.clay.neg.int[order(cu.clay.neg.int[,3], decreasing = TRUE), ]
  
  cu.depth.int <- cu.int[which(cu.int[,4] != 0), ]
  cu.depth.sorted.int <- cu.depth.int[order(cu.depth.int[,4]), ]
  cu.depth.sorted.pos.int <- cu.depth.sorted.int[which(mars.dirs[cu.depth.sorted.int[,5], 4]==1), ]
  cu.depth.neg.int <- cu.depth.sorted.int[which(mars.dirs[cu.depth.sorted.int[,5], 4]==-1), ]
  cu.depth.sorted.neg.int <- cu.depth.neg.int[order(cu.depth.neg.int[,4], decreasing = TRUE), ]
  
  
  cu.list <- list("cu.gas.sorted.pos"=cu.gas.sorted.pos, "cu.gas.sorted.neg"=cu.gas.sorted.neg, 
                  "cu.gas.sorted.pos.int"=cu.gas.sorted.pos.int, "cu.gas.sorted.neg.int"=cu.gas.sorted.neg.int, 
                  "cu.oil.sorted.pos"=cu.oil.sorted.pos, "cu.oil.sorted.neg"=cu.oil.sorted.neg,
                  "cu.oil.sorted.pos.int"=cu.oil.sorted.pos.int, "cu.oil.sorted.neg.int"=cu.oil.sorted.neg.int,
                  "cu.clay.sorted.pos"=cu.clay.sorted.pos, "cu.clay.sorted.neg"=cu.clay.sorted.neg,
                  "cu.clay.sorted.pos.int"=cu.clay.sorted.pos.int, "cu.clay.sorted.neg.int"=cu.clay.sorted.neg.int,
                  "cu.depth.sorted.pos.int"=cu.depth.sorted.pos.int, "cu.depth.sorted.neg.int"=cu.depth.sorted.neg.int,
                  "idx.non.interactions"= idx.non.interactions, "idx.interactions"=idx.interactions, "M"=M)
  return(cu.list)
}



## Derivatives for-loops (slow) ----
mars.derivative <- function(x, mars.model, mars.cuts, mars.dirs){
  # Easy to interperet but slow
  
  # Empty vector for derivatives
  derivative <- rep(0,length(x))
  
  for (m in 2:length(mars.model$coefficients)){
    
    # which variable has cuts at term m
    non.zero <- which(mars.dirs[m,] != 0) # which(mars.cuts[m,] != 0) 
    
    # idx where the conditions x_j > t or -x_j > -t  are met 
    idx <- which((x*mars.dirs[m,]) > (mars.cuts[m,] * mars.dirs[m,]))
    
    # all conditions are met 
    if (length(idx) == length(non.zero)){ 
      
      # interactions
      if (length(idx) > 1){ 
        
        
        for (j in idx){
          dirs.aux <- mars.dirs[m,]
          dirs.aux[j] <- 0
          
          # [abs(x_k-t_k)*dirs] k  = 1, 2, 3, 4 and k ≠ j
          diff <- (x - mars.cuts[m,])*dirs.aux # abs(x - mars.cuts[m,])*dirs.aux
          
          # derivative += beta_m * (1 or -1) * (x_k-t_k) *  ..., k  = 1, 2, 3, 4 and k ≠ j
          derivative[j] <- derivative[j] + mars.model$coefficients[m]*mars.dirs[m,j]*prod(diff[which(diff != 0)])
        }
      }
      
      # no interactions
      else if (length(idx) == 1){
        derivative[idx] <- derivative[idx] + mars.model$coefficients[m]*mars.dirs[m,idx]
      }
    }
  }
  return(derivative[1:3]) 
}

mars.derivative.X <- function(X, mars.model, mars.cuts, mars.dirs){
  d <- apply(X, MARGIN = 1, FUN = function(x) mars.derivative(x, mars.model, mars.cuts, mars.dirs)) 
  return(d)
}

mars.derivative.2D <- function(x, mars.model){
  mars.cuts <- mars.model$cuts[mars.model$selected.terms,]
  mars.dirs <- mars.model$dirs[mars.model$selected.terms,]
  
  # Empty vector for derivatives
  derivative <- rep(0,length(x))
  
  for (m in 2:length(mars.model$coefficients)){
    # m <- 17
    # which variable has cuts at term m
    non.zero <- which(mars.cuts[m,] != 0) 
    # are the conditions x_j > t or -x_j > -t met
    idx <- which((x*mars.dirs[m,]) > (mars.cuts[m,] * mars.dirs[m,]))
    
    # all conditions are met 
    if (length(idx) == length(non.zero)){ 
      
      # interactions
      if (length(idx) > 1){ 
        
        for (j in idx){
          dirs.aux <- (mars.dirs[m,]) # abs(mars.dirs[m,])
          dirs.aux[j] <- 0
          
          # [(x_k-t_k)]*dirs  k  = 1, 2, 3, 4 and k \neq j
          diff <- (x - mars.cuts[m,])*dirs.aux # abs(x - mars.cuts[m,])*dirs.aux
          
          #b <- 13.423401
          #c(-b*())
          
          # derivative += beta_m * (1 or -1) * (x_k-t_k) *  ..., 
          # k  = 1, 2, 3, 4 and k \neq j
          derivative[j] <- derivative[j] + mars.model$coefficients[m]*mars.dirs[m,j]*prod(diff[which(diff != 0)])
        }
      }
      
      # no interactions
      else if (length(idx) == 1){
        derivative[idx] <- derivative[idx] + 
          mars.model$coefficients[m]*mars.dirs[m,idx]
      }
      
    }
  }
  return(derivative) 
}


## Derivatives matrix multiplication ----
# Faster version of mars.derivative
mars.derivative.matrixmult <- function(x, mars.model, mars.cuts, mars.dirs){
  
  # Empty vector for derivatives
  derivative <- rep(0,4)
  
  M <-  length(mars.model$coefficients)
  
  # (Mx4) which variable has cuts at term m 
  non.zero <- t(apply(X=mars.dirs, MARGIN=1, function(md) {md != 0} ))
  
  # (Mx4) idx where the conditions x_j > t or -x_j > -t  are met
  x.dirs <- mars.dirs * outer(rep(1,M), x)#[,1,]
  idx <- as.matrix((x.dirs > (mars.cuts * mars.dirs)))
  
  # all conditions are met for index m
  m <- which(rowSums(non.zero == idx) == 4)
  
  # index for non-interactions terms 
  ni <- which(rowSums(idx) == 1)
  m.ni <- intersect(m, ni)
  
  # non-interactions terms 
  derivative <- derivative + colSums(mars.model$coefficients[m.ni]*mars.dirs[m.ni,])
  
  # index for interaction terms
  i <- which(rowSums(idx) > 1)
  m.i <- intersect(m, i)
  
  # interaction terms This part might be improved in speed
  for(i in 1:length(m.i)){
    J <- which(idx[m.i[i],])
    for (j in J){
      dirs.aux <- mars.dirs[m.i[i],]
      dirs.aux[j] <- 0
      diff <- (x - mars.cuts[m.i[i],])*dirs.aux
      derivative[j] <- derivative[j] + mars.model$coefficients[m.i[i]]*mars.dirs[m.i[i],j]*prod(diff[which(diff != 0)])
    }
  }
  
  return(derivative[1:3]) 
}


mars.derivative.matrixmult.X <- function(X, mars.model, mars.cuts, mars.dirs){
  d <- apply(X, MARGIN = 1, function(x) mars.derivative.matrixmult(x, mars.model, mars.cuts, mars.dirs)) 
  return(d)
}



## Derivatives binary search interactions ----
mars.derivative.binary <- function(x, mars.model, cu.list){

  n <- length(x[,1])

  # Empty vector for derivatives
  derivatives <- matrix(0, nrow=n, ncol=4)
  
  # number of terms in the model
  M <- cu.list$M
  
  cuts <- cbind(getCutsDirs(mars.model)$cuts, 1:M)
  dirs <- cbind(getCutsDirs(mars.model)$dirs, 1:M)
  m.interactions <- cu.list$idx.interactions
  
  # sorted cuts, pos means (x-t)+ and neg mens (t-x)+
  cu.gas.sorted.pos <- cu.list$cu.gas.sorted.pos
  cu.gas.sorted.neg <- cu.list$cu.gas.sorted.neg
  cu.oil.sorted.pos <- cu.list$cu.oil.sorted.pos
  cu.oil.sorted.neg <- cu.list$cu.oil.sorted.neg
  cu.clay.sorted.pos <- cu.list$cu.clay.sorted.pos
  cu.clay.sorted.neg <- cu.list$cu.clay.sorted.neg
  
  cu.gas.sorted.pos.int <- cu.list$cu.gas.sorted.pos.int
  cu.gas.sorted.neg.int <- cu.list$cu.gas.sorted.neg.int
  cu.oil.sorted.pos.int <- cu.list$cu.oil.sorted.pos.int
  cu.oil.sorted.neg.int <- cu.list$cu.oil.sorted.neg.int
  cu.clay.sorted.pos.int <- cu.list$cu.clay.sorted.pos.int
  cu.clay.sorted.neg.int <- cu.list$cu.clay.sorted.neg.int
  cu.depth.sorted.pos.int <- cu.list$cu.depth.sorted.pos.int
  cu.depth.sorted.neg.int <- cu.list$cu.depth.sorted.neg.int
  
  
  # non-interactions gas
  m.gas.pos <- findInterval(x[,1], cu.gas.sorted.pos[,1], left.open = TRUE)
  m.gas.neg <- findInterval(-x[,1], -cu.gas.sorted.neg[,1], left.open = TRUE)
  
  derivatives[,1] = derivatives[,1] + sapply(m.gas.pos, function(m){
    if (m > 0){
      return(sum(mars.model$coefficients[cu.gas.sorted.pos[c(1:m),5]]))
    }
    else{return(0)}})
  
  derivatives[,1] = derivatives[,1] + sapply(m.gas.neg, function(m){
    if (m > 0){
      return(sum( - mars.model$coefficients[cu.gas.sorted.neg[c(1:m),5]]))
    }
    else{return(0)}})

  
  # non-interactions oil
  m.oil.pos <- findInterval(x[,2], cu.oil.sorted.pos[,2], left.open = TRUE)
  m.oil.neg <- findInterval(-x[,2], -cu.oil.sorted.neg[,2], left.open = TRUE)

  derivatives[,2] = derivatives[,2] + sapply(m.oil.pos, function(m){
    if (m > 0){
      return(sum(mars.model$coefficients[cu.oil.sorted.pos[c(1:m),5]]))
    }
    else{return(0)}})
  
  derivatives[,2] = derivatives[,2] + sapply(m.oil.neg, function(m){
    if (m > 0){
      return(sum(-mars.model$coefficients[cu.oil.sorted.neg[c(1:m),5]]))
    }
    else{return(0)}})
  
  
  # non-interactions clay
  m.clay.pos <- findInterval(x[,3], cu.clay.sorted.pos[,3], left.open = TRUE)
  m.clay.neg <- findInterval(-x[,3], -cu.clay.sorted.neg[,3], left.open = TRUE)
  
  derivatives[,3] = derivatives[,3] + sapply(m.clay.pos, function(m){
    if (m > 0){
      return(sum(mars.model$coefficients[cu.clay.sorted.pos[c(1:m),5]]))
    }
    else{return(0)}})
  
  derivatives[,3] = derivatives[,3] + sapply(m.clay.neg, function(m){
    if (m > 0){
      return(sum(-mars.model$coefficients[cu.clay.sorted.neg[c(1:m),5]]))
    }
    else{return(0)}})
  
  
  # interactions
  m.gas.pos.int.idx=0
  m.oil.pos.int.idx=0
  m.clay.pos.int.idx=0
  m.depth.pos.int.idx=0
  m.gas.neg.int.idx=0
  m.oil.neg.int.idx=0
  m.clay.neg.int.idx=0 
  m.depth.neg.int.idx=0
  
  if (1 %in% dirs[m.interactions,1]){
    m.gas.pos.int <- findInterval(x[,1], cu.gas.sorted.pos.int[,1])
    m.gas.pos.int.idx <- sapply(X=m.gas.pos.int, function(m){
      if (m > 0){
        cu.gas.sorted.pos.int[c(1:m),5]}
    })}

  if (-1 %in% dirs[m.interactions,1]){
  m.gas.neg.int <- findInterval(-x[,1], -cu.gas.sorted.neg.int[,1])
  m.gas.neg.int.idx <- sapply(X=m.gas.neg.int, function(m){
    if (m > 0){
    cu.gas.sorted.neg.int[c(1:m),5]}
  })}
  
  if (1 %in% dirs[m.interactions,2]){
  m.oil.pos.int <- findInterval(x[,2], cu.oil.sorted.pos.int[,2])
  m.oil.pos.int.idx <- sapply(X=m.oil.pos.int, function(m){
    if (m > 0){
    cu.oil.sorted.pos.int[c(1:m),5]}
  })}

  if (-1 %in% dirs[m.interactions,2]){
  m.oil.neg.int <- findInterval(-x[,2], -cu.oil.sorted.neg.int[,2])
  m.oil.neg.int.idx <- sapply(X=m.oil.neg.int, function(m){
    if (m > 0){
    cu.oil.sorted.neg.int[c(1:m),5]}
  })}

  if (1 %in% dirs[m.interactions,3]){
  m.clay.pos.int <- findInterval(x[,3], cu.clay.sorted.pos.int[,3])
  m.clay.pos.int.idx <- sapply(X=m.clay.pos.int, function(m){
    if (m > 0){
    cu.clay.sorted.pos.int[c(1:m),5]}
  })}

  if (-1 %in% dirs[m.interactions,3]){
  m.clay.neg.int <- findInterval(-x[,3], -cu.clay.sorted.neg.int[,3])
  m.clay.neg.int.idx <- sapply(X=m.clay.neg.int, function(m){
    if (m > 0){
    cu.clay.sorted.neg.int[c(1:m),5]}
  })}

  if (1 %in% dirs[m.interactions,4]){
  m.depth.pos.int <- findInterval(x[,4], cu.depth.sorted.pos.int[,4])
  m.depth.pos.int.idx <- sapply(X=m.depth.pos.int, function(m){
    if (m > 0){
    cu.depth.sorted.pos.int[c(1:m),5]}
  })}

  if (-1 %in% dirs[m.interactions,4]){
  m.depth.neg.int <- findInterval(-x[,4], -cu.depth.sorted.neg.int[,4])
  m.depth.neg.int.idx <- sapply(X=m.depth.neg.int, function(m){
    if (m > 0){
    cu.depth.sorted.neg.int[c(1:m),5]}
  })}


  J.pos <- list(m.gas.pos.int.idx, m.oil.pos.int.idx, m.clay.pos.int.idx, m.depth.pos.int.idx)
  J.neg <- list(m.gas.neg.int.idx, m.oil.neg.int.idx, m.clay.neg.int.idx, m.depth.neg.int.idx)
  
  
  for (m in m.interactions){
    D <- dirs[m,1:4]
    C <- cuts[m,1:4]
    j.pos <- which(D==1)
    j.neg <- which(D==-1)
    
    if(length(j.pos) > 0){
      L.pos <- sapply(j.pos, function(j){
        sapply(J.pos[[j]], function(J){
          fmatch(m,J)})
      })}
    
    if(length(j.neg) > 0){
    L.neg <- sapply(j.neg, function(j){
      sapply(J.neg[[j]], function(J){
        fmatch(m,J)})
    })
    }
    
    if (length(j.pos) > 0 && length(j.neg) > 0){
      L <- cbind(L.pos, L.neg)
    } else if (length(j.pos) > 0){
      L <- L.pos
    } else{
      L <- L.neg
      }
    
    
    # indekser hvor alle identity functions ((xj > t) osv..) er oppfyllt
    Idx <- which(!is.na(rowSums(L)))
    J.idx <- c(j.pos, j.neg)
    
    for (j in J.idx){
      dirs.aux <- dirs[m,1:4]
      dirs.aux[j] <- 0

      k <- which(dirs.aux!=0)

      X <- sweep(as.matrix(x[,k],nrow=n), MARGIN=2, cuts[m,k], "-")
      diff <- sweep(X, MARGIN=2, dirs.aux[k], '*')
      diff <- cbind(diff, rep(1,n))
      derivatives[Idx,j] <- derivatives[Idx,j] + mars.model$coefficients[m]*dirs[m,j]*apply(diff[Idx,], 1, prod)
    }
    

  }
  
  return(derivatives[,1:3]) 
}



# Systematic errors ----
systematic.error.mars <- function(){
  library(latex2exp)
  
  system.time(mars.R0.bf.pred <- predict(object= mars.R0.bf, newdata = df.test)) 
  system.time(mars.G.bf.pred <- predict(object= mars.G.bf, newdata = df.test)) 
  mars.R0.bf.pred <- mars.R0.pred
  mars.G.bf.pred <- mars.G.pred
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xg_R0.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.gas, df.test$R0 - mars.R0.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_g$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xo_R0.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.oil, df.test$R0 - mars.R0.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_o$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xc_R0.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.clay, df.test$R0 - mars.R0.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_c$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xd_R0.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.depth, df.test$R0 - mars.R0.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_d$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xg_G.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.gas, df.test$G - mars.G.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_g$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xo_G.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.oil, df.test$G - mars.G.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_o$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xc_G.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.clay, df.test$G - mars.G.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_c$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
  pdf("/Users/karenauestad/Dokumenter/master/Alvheim_kode/h_hat_models/plots/error_xd_G.pdf", width=6, height=6)
  par(mar = c(5, 5, 4, 2))
  plot(df.test$x.depth, df.test$G - mars.G.bf.pred, col="brown2", pch=20, 
       xlab=TeX(r"($x_d$)"), ylab=TeX(r"($h - \widehat{h}$)"))
  abline(h=0, lty=3)
  dev.off()
  
}



