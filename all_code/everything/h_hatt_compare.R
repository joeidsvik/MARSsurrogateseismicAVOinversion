

source("h_hat_MARS.R")  
source("h_hat_NP.R")
source("pseudoBayes_aux.R")



library(latex2exp)

convrt <- 1.055
# Bechmark: 1.6 sec (Karens computer)
Average <- 10
sfm.time.vec <- rep(NaN, Average)
for (i in 1:Average){
     sfm.time.vec[i] <- system.time(sfm <- sim_forward_model(df.test[,1], df.test[,2], df.test[,3], df.test[,4]/convrt))[3]
}
sfm.time <- mean(sfm.time.vec)

sfm.R0 <- sfm[rep(c(TRUE, FALSE), length(sfm)/2)]
sfm.G <- sfm[rep(c(FALSE, TRUE), length(sfm)/2)]



# R0 ----
sfm.R0 <- df.test$R0

## MARS bf : R0 ----
#summary(mars.R0.bf)

time.mars.R0.bf.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.mars.R0.bf.vec[i] <- system.time(mars.R0.bf.pred <- predict(object= mars.R0.bf, newdata = df.test))[3]
}
time.mars.R0.bf <- mean(time.mars.R0.bf.vec)

corr.mars.bf.R0 <- corr(sfm.R0, mars.R0.bf.pred)
mse.mars.bf.R0 <- mean((sfm.R0 - mars.R0.bf.pred)^2)

#pdf("/Users/karenauestad/Dokumenter/master/plots_master/h_hat/sys_err_mars_bf_R0.pdf", width=6, height=6)
par(mar=c(5, 5, 4, 2))
plot(mars.R0.bf.pred, sfm.R0, pch = 21, col="brown2", xlim=c(-0.31, 0.15), ylim=c(-0.31, 0.15),
    xlab=TeX(r"($\widehat{h}_{MARS}$)"),  ylab=TeX(r'($h$)'), cex.lab=1.4, cex.axis =1.3)
arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)
#dev.off()


## MARS of : R0 ----
#summary(mars.R0.of)

time.mars.R0.of.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.mars.R0.of.vec[i] <- system.time(mars.R0.of.pred <- predict(object= mars.R0.of, newdata = df.test))[3]
}
time.mars.R0.of <- mean(time.mars.R0.of.vec)

corr.mars.of.R0 <- corr(sfm.R0, mars.R0.of.pred)
mse.mars.of.R0 <- mean((sfm.R0 - mars.R0.of.pred)^2)

#pdf("/Users/karenauestad/Dokumenter/master/plots_master/h_hat/sys_err_mars_of_R0.pdf", width=6, height=6)
par(mar=c(5, 5, 4, 2))
plot(mars.R0.of.pred, sfm.R0, pch = 21, col="brown2", xlim=c(-0.31, 0.15), ylim=c(-0.31, 0.15),
       xlab=TeX(r"($\widehat{h}_{MARS_{OF}}$)"),  ylab=TeX(r'($h$)'), cex.lab=1.4, cex.axis =1.3)
arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)
#dev.off() 

#pdf("/Users/karenauestad/Dokumenter/master/plots_master/h_hat/mars_of_R0_vs_mars_R0.pdf", width=6, height=6)
par(mar=c(5, 5, 4, 2))
plot(mars.R0.of.pred, mars.R0.bf.pred, pch = 21, col="dodgerblue", xlim=c(-0.31, 0.15), ylim=c(-0.31, 0.15),
     xlab=TeX(r"($\widehat{h}_{MARS_{OF}}$)"),  ylab=TeX(r'($\widehat{h}_{MARS}$)'), cex.lab=1.4, cex.axis =1.3)
arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)
#dev.off()

## NPR 1000 : R0 ----
#summary(npr.R0.1000)

time.npr.R0.1000.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.npr.R0.1000.vec[i] <- system.time(npr.R0.1000.pred <- predict(object= npr.R0.1000, newdata = df.test))[3]
}
time.npr.R0.1000 <- mean(time.npr.R0.1000.vec)

corr.npr.1000.R0 <- corr(sfm.R0, npr.R0.1000.pred)
mse.npr.1000.R0 <- mean((sfm.R0 - npr.R0.1000.pred)^2)

# par(mai=c(1, 1, 1, 1))
# plot(npr.R0.1000.pred, sfm.R0, pch = 21, col="brown2", xlim=c(-0.3, 0.3), ylim=c(-0.3, 0.3),
#       xlab=TeX(r"($\widehat{h}_{R_0}(x)$)"),  ylab=TeX(r'($h_{R_0}(x)$)'), cex.lab=2, cex.axis =2)
# arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)


## NPR 4000 : R0 ----
#summary(npr.R0.4000)

time.npr.R0.4000.vec <-  rep(NaN, Average)
for (i in 1:Average){
     time.npr.R0.4000.vec[i] <- system.time(npr.R0.4000.pred <- predict(object= npr.R0.4000, newdata = df.test))[3]
}
time.npr.R0.4000 <- mean(time.npr.R0.4000.vec)

corr.npr.4000.R0  <- corr(sfm.R0, npr.R0.4000.pred)
mse.npr.4000.R0 <- mean((sfm.R0 - npr.R0.4000.pred)^2)

# par(mai=c(1, 1, 1, 1))
# plot(npr.R0.4000.pred, sfm.R0, pch = 21, col="brown2", xlim=c(-0.3, 0.3), ylim=c(-0.3, 0.3),
#       xlab=TeX(r"($\widehat{h}_{R_0}(x)$)"),  ylab=TeX(r'($h_{R_0}(x)$)'), cex.lab=2, cex.axis =2)
# arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)


## Table ----
mse.vec <- c(mse.mars.bf.R0, mse.mars.of.R0, mse.npr.1000.R0, mse.npr.4000.R0)
corr.vec <- c(corr.mars.bf.R0, corr.mars.of.R0, corr.npr.1000.R0, corr.npr.4000.R0)
time.vec <- c(time.mars.R0.bf, time.mars.R0.of, time.npr.R0.1000, time.npr.R0.4000)
Table.R0 <- as.data.frame(matrix(cbind(corr.vec, mse.vec, time.vec), ncol=3), row.names =  c("mars bf", "mars of", "np 1000", "np 4000") )
colnames(Table.R0) <- c("corr", "mse", "time [sec]")
round(Table.R0, 6)




# G ----
sfm.G <- df.test$G


## MARS bf : G ----
#summary(mars.G.bf)

time.mars.G.bf.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.mars.G.bf.vec[i] <- time.mars.G.bf.vec[1] <- system.time(mars.G.bf.pred <- predict(object= mars.G.bf, newdata = df.test))[3]
}
time.mars.G.bf <- mean(time.mars.G.bf.vec)

corr.mars.bf.G <- corr(sfm.G, mars.G.bf.pred)
mse.mars.bf.G <- mean((sfm.G - mars.G.bf.pred)^2)

#pdf("/Users/karenauestad/Dokumenter/master/plots_master/h_hat/sys_err_mars_bf_G.pdf", width=6, height=6)
par(mar=c(5, 5, 4, 2))
plot(mars.G.bf.pred, sfm.G, pch = 21, col="brown2", xlim=c(-0.34, 0.15), ylim=c(-0.34, 0.15),
     xlab=TeX(r"($\widehat{h}_{MARS}$)"),  ylab=TeX(r'($h$)'), cex.lab=1.4, cex.axis=1.3)
arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)
#dev.off()


## MARS of : G ----
#summary(mars.G.of)

time.mars.G.of.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.mars.G.of.vec[i] <- system.time(mars.G.of.pred <- predict(object= mars.G.of, newdata = df.test))[3]
}
time.mars.G.of <- mean(time.mars.G.of.vec)

corr.mars.of.G <- corr(sfm.G, mars.G.of.pred)
mse.mars.of.G <- mean((sfm.G - mars.G.of.pred)^2)

#pdf("/Users/karenauestad/Dokumenter/master/plots_master/h_hat/sys_err_mars_of_G.pdf", width=6, height=6)
par(mar=c(5, 5, 4, 2))
plot(mars.G.of.pred, sfm.G, pch = 21, col="brown2", xlim=c(-0.34, 0.15), ylim=c(-0.34, 0.15),
      xlab=TeX(r"($\widehat{h}_{MARS_{OF}}$)"),  ylab=TeX(r'($h$)'), cex.lab=1.4, cex.axis =1.3)
arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)
#dev.off()


## NPR 1000 : G ----
#summary(npr.G.1000)

time.npr.G.1000.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.npr.G.1000.vec[i] <- system.time(npr.G.1000.pred <- predict(object= npr.G.1000, newdata = df.test))[3]
}
time.npr.G.1000 <- mean(time.npr.G.1000.vec)


corr.npr.1000.G <- corr(sfm.G, npr.G.1000.pred)
mse.npr.1000.G <- mean((sfm.G - npr.G.1000.pred)^2)

# par(mai=c(1, 1, 1, 1))
# plot(npr.G.1000.pred, sfm.G, pch = 21, col="brown2", xlim=c(-0.3, 0.3), ylim=c(-0.3, 0.3),
#       xlab=TeX(r"($\widehat{h}_{G}(x)$)"),  ylab=TeX(r'($h_{G}(x)$)'), cex.lab=2, cex.axis =2)
# arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)


## NPR 4000 : G ----
#summary(npr.G.4000)

time.npr.G.4000.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.npr.G.4000.vec[i] <- system.time(npr.G.4000.pred <- predict(object= npr.G.4000, newdata = df.test))[3]
}
time.npr.G.4000 <- mean(time.npr.G.4000.vec)

corr.npr.4000.G  <- corr(sfm.G, npr.G.4000.pred)
mse.npr.4000.G <- mean((sfm.G - npr.G.4000.pred)^2)

# par(mai=c(1, 1, 1, 1))
# plot(npr.G.4000.pred, sfm.G, pch = 21, col="brown2", xlim=c(-0.3, 0.3), ylim=c(-0.3, 0.3),
#      xlab=TeX(r"($\widehat{h}_{G}(x)$)"),  ylab=TeX(r'($h_{G}(x)$)'), cex.lab=2, cex.axis =2)
# arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)


## Table ----
mse.vec <- c(mse.mars.bf.G, mse.mars.of.G, mse.npr.1000.G, mse.npr.4000.G)
corr.vec <- c(corr.mars.bf.G, corr.mars.of.G, corr.npr.1000.G, corr.npr.4000.G)
time.vec <- c(time.mars.G.bf, time.mars.G.of, time.npr.G.1000, time.npr.G.4000)
Table.G <- as.data.frame(matrix(cbind(corr.vec, mse.vec, time.vec), ncol=3), row.names =  c("mars bf", "mars of", "np 1000", "np 4000") )
colnames(Table.G) <- c("corr", "mse", "time [sec]")
round(Table.G, 6)

# Gradient of R0 ----
# Pre-calculated numerical gradient (both R0 and G) using epsilon = 0.00001
load("h_hat_models/derivative_sfm_df_test.Rda")
X <- df.test[,1:4]

## MARS fitted to gradient : Gradient of R0 ----
load("h_hat_models/D_mars.Rda")

time.D.R0.mars.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.R0.mars.vec[i] <- system.time(D.R0 <- mars.derivatives.pred.R0(D.mars, df.test))[3]
}
time.D.R0.mars <- mean(time.D.R0.mars.vec)


corr.D.mars.R0.g <- corr(num.deriv.sfm$dR0dxg, predict(object = D.mars[[1]], newdata = df.test))
corr.D.mars.R0.o <- corr(num.deriv.sfm$dR0dxo, predict(object = D.mars[[2]], newdata = df.test))
corr.D.mars.R0.c <- corr(num.deriv.sfm$dR0dxc, predict(object = D.mars[[3]], newdata = df.test))
corr.D.mars.R0 <- mean(corr.D.mars.R0.g, corr.D.mars.R0.o, corr.D.mars.R0.c)

mse.D.mars.R0.g <- mean(((predict(object = D.mars[[1]], newdata = df.test)) - num.deriv.sfm$dR0dxg)^2) 
mse.D.mars.R0.o <- mean(((predict(object = D.mars[[2]], newdata = df.test)) - num.deriv.sfm$dR0dxo)^2) 
mse.D.mars.R0.c <- mean(((predict(object = D.mars[[3]], newdata = df.test)) - num.deriv.sfm$dR0dxc)^2)
mse.D.mars.R0 <- mean(mse.D.mars.R0.g, mse.D.mars.R0.o, mse.D.mars.R0.c)


# plot(df.test[,2], num.deriv.sfm[,2], col="brown2", ylim=c(-0.05, 0.02))
# points(df.test[,2], D.R0[,2], col="dodgerblue3")
# 
# par(mai=c(1, 1, 1, 1))
# plot(D.R0[,2], num.deriv.sfm[,2], pch = 21, col="brown2", xlim=c(-0.03, 0.03), ylim=c(-0.03, 0.03),
#      xlab="",  ylab="", cex.lab=2, cex.axis =2)
# arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)



## MARS bf binary : Gradient of R0 ----
time.D.mars.bf.R0.binary.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.mars.bf.R0.binary.vec[i] <- system.time(D.mars.bf.R0.binary <- mars.derivative.binary(x = X, mars.model = mars.R0.bf, cu.list=get.cu(mars.R0.bf)))[3]
}
time.D.mars.bf.R0.binary <- mean(time.D.mars.bf.R0.binary.vec)

corr.D.mars.bf.R0.binary.g <- corr(num.deriv.sfm[,1], D.mars.bf.R0.binary[,1])
corr.D.mars.bf.R0.binary.o <- corr(num.deriv.sfm[,2], D.mars.bf.R0.binary[,2])
corr.D.mars.bf.R0.binary.c <- corr(num.deriv.sfm[,3], D.mars.bf.R0.binary[,3])
corr.D.mars.bf.R0.binary <- mean(corr.D.mars.bf.R0.binary.g, corr.D.mars.bf.R0.binary.o, corr.D.mars.bf.R0.binary.c)

mse.D.mars.bf.R0.binary.g <- mean((num.deriv.sfm[,1]- D.mars.bf.R0.binary[,1])^2)
mse.D.mars.bf.R0.binary.o <- mean((num.deriv.sfm[,2]- D.mars.bf.R0.binary[,2])^2)
mse.D.mars.bf.R0.binary.c <- mean((num.deriv.sfm[,3]- D.mars.bf.R0.binary[,3])^2)
mse.D.mars.bf.R0.binary <- mean(mse.D.mars.bf.R0.binary.g, mse.D.mars.bf.R0.binary.o, mse.D.mars.bf.R0.binary.c) 


# If statements to show how much time is reduced by binary
time.D.mars.bf.R0.if.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.mars.bf.R0.if.vec[i] <- system.time(D.mars.bf.R0.if <- mars.derivative.X(X, mars.R0.bf, getCutsDirs(mars.R0.bf)$cuts, getCutsDirs(mars.R0.bf)$dirs))[3]
}
time.D.mars.bf.R0.if <- mean(time.D.mars.bf.R0.if.vec)

#(the same metrics)
corr.D.mars.bf.R0.if.g <- corr(num.deriv.sfm[,1], D.mars.bf.R0.if[1,])
corr.D.mars.bf.R0.if.o <- corr(num.deriv.sfm[,2], D.mars.bf.R0.if[2,])
corr.D.mars.bf.R0.if.c <- corr(num.deriv.sfm[,3], D.mars.bf.R0.if[3,])
corr.D.mars.bf.R0.if <- mean(corr.D.mars.bf.R0.if.g, corr.D.mars.bf.R0.if.o, corr.D.mars.bf.R0.if.c)

mse.D.mars.bf.R0.if.g <- mean((num.deriv.sfm[,1]- D.mars.bf.R0.if[1,])^2)
mse.D.mars.bf.R0.if.o <- mean((num.deriv.sfm[,2]- D.mars.bf.R0.if[2,])^2)
mse.D.mars.bf.R0.if.c <- mean((num.deriv.sfm[,3]- D.mars.bf.R0.if[3,])^2)
mse.D.mars.bf.R0.if <- mean(mse.D.mars.bf.R0.if.g, mse.D.mars.bf.R0.if.o, mse.D.mars.bf.R0.if.c) 



# worst (oil)
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.mars.bf.R0.binary[,2], col="dodgerblue3")


## MARS of binary : Gradient of R0 ----
time.D.mars.of.R0.binary.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.mars.of.R0.binary.vec[i] <- system.time(D.mars.of.R0.binary <- mars.derivative.binary(x = X, mars.model = mars.R0.of, cu.list=get.cu(mars.R0.of)))[3]
}
time.D.mars.of.R0.binary <- mean(time.D.mars.of.R0.binary.vec)


corr.D.mars.of.R0.binary.g <- corr(num.deriv.sfm[,1], D.mars.of.R0.binary[,1])
corr.D.mars.of.R0.binary.o <- corr(num.deriv.sfm[,2], D.mars.of.R0.binary[,2])
corr.D.mars.of.R0.binary.c <- corr(num.deriv.sfm[,3], D.mars.of.R0.binary[,3])
corr.D.mars.of.R0.binary <- mean(corr.D.mars.of.R0.binary.g, corr.D.mars.of.R0.binary.o, corr.D.mars.of.R0.binary.c)

mse.D.mars.of.R0.binary.g <- mean((num.deriv.sfm[,1]- D.mars.of.R0.binary[,1])^2)
mse.D.mars.of.R0.binary.o <- mean((num.deriv.sfm[,2]- D.mars.of.R0.binary[,2])^2)
mse.D.mars.of.R0.binary.c <- mean((num.deriv.sfm[,3]- D.mars.of.R0.binary[,3])^2)
mse.D.mars.of.R0.binary <- mean(mse.D.mars.of.R0.binary.g, mse.D.mars.of.R0.binary.o, mse.D.mars.of.R0.binary.c) 

# If statements to show how much time is reduced by binary
time.D.mars.of.R0.if.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.mars.of.R0.if.vec[i] <- system.time(D.mars.of.R0.if <- mars.derivative.X(X, mars.R0.of, getCutsDirs(mars.R0.of)$cuts, getCutsDirs(mars.R0.of)$dirs))[3]
}
time.D.mars.of.R0.if <- mean(time.D.mars.of.R0.if.vec)

#(the same metrics)
corr.D.mars.of.R0.if.g <- corr(num.deriv.sfm[,1], D.mars.of.R0.if[1,])
corr.D.mars.of.R0.if.o <- corr(num.deriv.sfm[,2], D.mars.of.R0.if[2,])
corr.D.mars.of.R0.if.c <- corr(num.deriv.sfm[,3], D.mars.of.R0.if[3,])
corr.D.mars.of.R0.if <- mean(corr.D.mars.of.R0.if.g, corr.D.mars.of.R0.if.o, corr.D.mars.of.R0.if.c)

mse.D.mars.of.R0.if.g <- mean((num.deriv.sfm[,1]- D.mars.of.R0.if[1,])^2)
mse.D.mars.of.R0.if.o <- mean((num.deriv.sfm[,2]- D.mars.of.R0.if[2,])^2)
mse.D.mars.of.R0.if.c <- mean((num.deriv.sfm[,3]- D.mars.of.R0.if[3,])^2)
mse.D.mars.of.R0.if <- mean(mse.D.mars.of.R0.if.g, mse.D.mars.of.R0.if.o, mse.D.mars.of.R0.if.c) 


# worst (oil)
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.mars.of.R0.binary[,2], col="dodgerblue3")



## NPR 1000 analytical : Gradient of R0 ----
time.D.npr.R0.1000.analytical.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.R0.1000.analytical.vec[i] <- system.time(D.npr.R0.1000.analytical <- derivative.npreg.anal(X, npr.R0.1000))[3]
}
time.D.npr.R0.1000.analytical <- mean(time.D.npr.R0.1000.analytical.vec)

corr.D.npr.R0.1000.g <- corr(num.deriv.sfm[,1], D.npr.R0.1000.analytical[1,])
corr.D.npr.R0.1000.o <- corr(num.deriv.sfm[,2], D.npr.R0.1000.analytical[2,])
corr.D.npr.R0.1000.c <- corr(num.deriv.sfm[,3], D.npr.R0.1000.analytical[3,])
corr.D.npr.R0.1000 <- mean(corr.D.npr.R0.1000.g, corr.D.npr.R0.1000.o, corr.D.npr.R0.1000.c)

mse.D.npr.R0.1000.g <- mean((num.deriv.sfm[,1]- D.npr.R0.1000.analytical[1,])^2)
mse.D.npr.R0.1000.o <- mean((num.deriv.sfm[,2]- D.npr.R0.1000.analytical[2,])^2)
mse.D.npr.R0.1000.c <- mean((num.deriv.sfm[,3]- D.npr.R0.1000.analytical[3,])^2)
mse.D.npr.R0.1000 <- mean(mse.D.npr.R0.1000.g, mse.D.npr.R0.1000.o, mse.D.npr.R0.1000.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.R0.1000.analytical[2,], col="dodgerblue3")



## NPR 1000 K means : Gradient of R0 ----
km.classif.R0.1000 <- kmeans(x=npr.R0.1000$eval, centers=10)
X.cl.R0.1000 <- clusterX(npr.R0.1000, km.classif.R0.1000)

time.D.npr.R0.1000.kmeans.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.R0.1000.kmeans.vec[i] <- system.time(D.npr.R0.1000.kmeans <- derivative.npreg.mult.data.points(X, X.cl.R0.1000, s = 5, km.classif.R0.1000, npr.R0.1000))[3]
}
time.D.npr.R0.1000.kmeans <- mean(time.D.npr.R0.1000.kmeans.vec)

corr.D.npr.R0.1000.kmeans.g <- corr(num.deriv.sfm[,1], D.npr.R0.1000.kmeans[1,])
corr.D.npr.R0.1000.kmeans.o <- corr(num.deriv.sfm[,2], D.npr.R0.1000.kmeans[2,])
corr.D.npr.R0.1000.kmeans.c <- corr(num.deriv.sfm[,3], D.npr.R0.1000.kmeans[3,])
corr.D.npr.R0.1000.kmeans <- mean(corr.D.npr.R0.1000.kmeans.g, corr.D.npr.R0.1000.kmeans.o, corr.D.npr.R0.1000.kmeans.c)

mse.D.npr.R0.1000.kmeans.g <- mean((num.deriv.sfm[,1]- D.npr.R0.1000.kmeans[1,])^2)
mse.D.npr.R0.1000.kmeans.o <- mean((num.deriv.sfm[,2]- D.npr.R0.1000.kmeans[2,])^2)
mse.D.npr.R0.1000.kmeans.c <- mean((num.deriv.sfm[,3]- D.npr.R0.1000.kmeans[3,])^2)
mse.D.npr.R0.1000.kmeans <- mean(mse.D.npr.R0.1000.kmeans.g, mse.D.npr.R0.1000.kmeans.o, mse.D.npr.R0.1000.kmeans.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.R0.1000.kmeans[2,], col="dodgerblue3")





## NPR 4000 analytical : Gradient of R0 ----
time.D.npr.R0.4000.analytical.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.R0.4000.analytical.vec[i] <- system.time(D.npr.R0.4000.analytical <- derivative.npreg.anal(X, npr.R0.4000))[3]
}
time.D.npr.R0.4000.analytical <- mean(time.D.npr.R0.4000.analytical.vec)

corr.D.npr.R0.4000.g <- corr(num.deriv.sfm[,1], D.npr.R0.4000.analytical[1,])
corr.D.npr.R0.4000.o <- corr(num.deriv.sfm[,2], D.npr.R0.4000.analytical[2,])
corr.D.npr.R0.4000.c <- corr(num.deriv.sfm[,3], D.npr.R0.4000.analytical[3,])
corr.D.npr.R0.4000 <- mean(corr.D.npr.R0.4000.g, corr.D.npr.R0.4000.o, corr.D.npr.R0.4000.c)

mse.D.npr.R0.4000.g <- mean((num.deriv.sfm[,1]- D.npr.R0.4000.analytical[1,])^2)
mse.D.npr.R0.4000.o <- mean((num.deriv.sfm[,2]- D.npr.R0.4000.analytical[2,])^2)
mse.D.npr.R0.4000.c <- mean((num.deriv.sfm[,3]- D.npr.R0.4000.analytical[3,])^2)
mse.D.npr.R0.4000 <- mean(mse.D.npr.R0.4000.g, mse.D.npr.R0.4000.o, mse.D.npr.R0.4000.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.R0.4000.analytical[2,], col="dodgerblue3")


## NPR 4000 K means : Gradient of R0 ----
km.classif.R0.4000 <- kmeans(x=npr.R0.4000$eval, centers=10)
X.cl.R0.4000 <- clusterX(npr.R0.4000, km.classif.R0.4000)

time.D.npr.R0.4000.kmeans.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.R0.4000.kmeans.vec[i] <- system.time(D.npr.R0.4000.kmeans <- derivative.npreg.mult.data.points(X, X.cl.R0.4000, s = 5, km.classif.R0.4000, npr.R0.4000))[3]
}
time.D.npr.R0.4000.kmeans <- mean(time.D.npr.R0.4000.kmeans.vec)

corr.D.npr.R0.4000.kmeans.g <- corr(num.deriv.sfm[,1], D.npr.R0.4000.kmeans[1,])
corr.D.npr.R0.4000.kmeans.o <- corr(num.deriv.sfm[,2], D.npr.R0.4000.kmeans[2,])
corr.D.npr.R0.4000.kmeans.c <- corr(num.deriv.sfm[,3], D.npr.R0.4000.kmeans[3,])
corr.D.npr.R0.4000.kmeans <- mean(corr.D.npr.R0.4000.kmeans.g, corr.D.npr.R0.4000.kmeans.o, corr.D.npr.R0.4000.kmeans.c)

mse.D.npr.R0.4000.kmeans.g <- mean((num.deriv.sfm[,1]- D.npr.R0.4000.kmeans[1,])^2)
mse.D.npr.R0.4000.kmeans.o <- mean((num.deriv.sfm[,2]- D.npr.R0.4000.kmeans[2,])^2)
mse.D.npr.R0.4000.kmeans.c <- mean((num.deriv.sfm[,3]- D.npr.R0.4000.kmeans[3,])^2)
mse.D.npr.R0.4000.kmeans <- mean(mse.D.npr.R0.4000.kmeans.g, mse.D.npr.R0.4000.kmeans.o, mse.D.npr.R0.4000.kmeans.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.R0.4000.kmeans[2,], col="dodgerblue3")


## Table ----
mse.vec <- c(mse.D.mars.R0, mse.D.mars.bf.R0.binary, mse.D.mars.of.R0.binary, 
             mse.D.npr.R0.1000, mse.D.npr.R0.1000.kmeans, 
             mse.D.npr.R0.4000, mse.D.npr.R0.4000.kmeans)
corr.vec <- c(corr.D.mars.R0, corr.D.mars.bf.R0.binary, corr.D.mars.of.R0.binary, 
              corr.D.npr.R0.1000, corr.D.npr.R0.1000.kmeans, 
              corr.D.npr.R0.4000, corr.D.npr.R0.4000.kmeans)
time.vec <- c(time.D.R0.mars, time.D.mars.bf.R0.binary, time.D.mars.of.R0.binary, 
              time.D.npr.R0.1000.analytical, time.D.npr.R0.1000.kmeans, 
              time.D.npr.R0.4000.analytical, time.D.npr.R0.4000.kmeans)
Table.D.R0 <- as.data.frame(matrix(cbind(corr.vec, mse.vec, time.vec), ncol=3), row.names =  c("mars", "mars bf",  "mars of", "np 1000", "np 1000 kmeans", "np 4000", "np 4000 kmeans") )
colnames(Table.D.R0) <- c("corr", "mse", "time [sec]")
round(Table.D.R0, 6)



# Gradient of G ----

# Pre-calculated numerical gradient (both R0 and G) using epsilon = 0.00001
load("h_hat_models/derivative_sfm_df_test.Rda")
X <- df.test[,1:4]

## MARS fitted to gradient : Gradient of G ----
# D.mars <- fit.mars.derivative()
# save(D.mars, file="D_mars.Rda")
load("h_hat_models/D_mars.Rda")


time.D.G.mars.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.G.mars.vec[i] <- system.time(D.G <- mars.derivatives.pred.G(D.mars, df.test))[3]
}
time.D.G.mars <- mean(time.D.G.mars.vec)

corr.D.mars.G.g <- corr(num.deriv.sfm$dGdxg, predict(object = D.mars[[4]], newdata = df.test))
corr.D.mars.G.o <- corr(num.deriv.sfm$dGdxo, predict(object = D.mars[[5]], newdata = df.test))
corr.D.mars.G.c <- corr(num.deriv.sfm$dGdxc, predict(object = D.mars[[6]], newdata = df.test))
corr.D.mars.G <- mean(corr.D.mars.G.g, corr.D.mars.G.o, corr.D.mars.G.c)

mse.D.mars.G.g <- mean(((predict(object = D.mars[[4]], newdata = df.test)) - num.deriv.sfm$dGdxg)^2) 
mse.D.mars.G.o <- mean(((predict(object = D.mars[[5]], newdata = df.test)) - num.deriv.sfm$dGdxo)^2) 
mse.D.mars.G.c <- mean(((predict(object = D.mars[[6]], newdata = df.test)) - num.deriv.sfm$dGdxc)^2)
mse.D.mars.G <- mean(mse.D.mars.G.g, mse.D.mars.G.o, mse.D.mars.G.c)


# plot(df.test[,2], num.deriv.sfm[,2], col="brown2", ylim=c(-0.05, 0.02))
# points(df.test[,2], D.G[,2], col="dodgerblue3")
# 
# par(mai=c(1, 1, 1, 1))
# plot(D.G[,2], num.deriv.sfm[,2], pch = 21, col="brown2", xlim=c(-0.03, 0.03), ylim=c(-0.03, 0.03),
#      xlab="",  ylab="", cex.lab=2, cex.axis =2)
# arrows(x0=-1, y0=-1, x1=1, y1=1, code=0, col="black", lwd=2)



## MARS bf binary : Gradient of G ----

time.D.mars.bf.G.binary.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.mars.bf.G.binary.vec[i] <- system.time(D.mars.bf.G.binary <- mars.derivative.binary(x = X, mars.model = mars.G.bf, cu.list=get.cu(mars.G.bf)))[3]
}
time.D.mars.bf.G.binary <- mean(time.D.mars.bf.G.binary.vec)

corr.D.mars.bf.G.binary.g <- corr(num.deriv.sfm$dGdxg, D.mars.bf.G.binary[,1])
corr.D.mars.bf.G.binary.o <- corr(num.deriv.sfm$dGdxo, D.mars.bf.G.binary[,2])
corr.D.mars.bf.G.binary.c <- corr(num.deriv.sfm$dGdxc, D.mars.bf.G.binary[,3])
corr.D.mars.bf.G.binary <- mean(corr.D.mars.bf.G.binary.g, corr.D.mars.bf.G.binary.o, corr.D.mars.bf.G.binary.c)

mse.D.mars.bf.G.binary.g <- mean((num.deriv.sfm$dGdxg- D.mars.bf.G.binary[,1])^2)
mse.D.mars.bf.G.binary.o <- mean((num.deriv.sfm$dGdxo- D.mars.bf.G.binary[,2])^2)
mse.D.mars.bf.G.binary.c <- mean((num.deriv.sfm$dGdxc- D.mars.bf.G.binary[,3])^2)
mse.D.mars.bf.G.binary <- mean(mse.D.mars.bf.G.binary.g, mse.D.mars.bf.G.binary.o, mse.D.mars.bf.G.binary.c) 

# worst (oil)
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.mars.bf.G.binary[,2], col="dodgerblue3")


## MARS of binary : Gradient of G ----
time.D.mars.of.G.binary.vec <-  rep(NaN, Average)
for (i in 1:Average){
     time.D.mars.of.G.binary.vec[i] <- system.time(D.mars.of.G.binary <- mars.derivative.binary(x = X, mars.model = mars.G.of, cu.list=get.cu(mars.G.of)))[3]
}
time.D.mars.of.G.binary <- mean(time.D.mars.of.G.binary.vec)


corr.D.mars.of.G.binary.g <- corr(num.deriv.sfm$dGdxg, D.mars.of.G.binary[,1])
corr.D.mars.of.G.binary.o <- corr(num.deriv.sfm$dGdxo, D.mars.of.G.binary[,2])
corr.D.mars.of.G.binary.c <- corr(num.deriv.sfm$dGdxc, D.mars.of.G.binary[,3])
corr.D.mars.of.G.binary <- mean(corr.D.mars.of.G.binary.g, corr.D.mars.of.G.binary.o, corr.D.mars.of.G.binary.c)

mse.D.mars.of.G.binary.g <- mean((num.deriv.sfm$dGdxg - D.mars.of.G.binary[,1])^2)
mse.D.mars.of.G.binary.o <- mean((num.deriv.sfm$dGdxo - D.mars.of.G.binary[,2])^2)
mse.D.mars.of.G.binary.c <- mean((num.deriv.sfm$dGdxc- D.mars.of.G.binary[,3])^2)
mse.D.mars.of.G.binary <- mean(mse.D.mars.of.G.binary.g, mse.D.mars.of.G.binary.o, mse.D.mars.of.G.binary.c) 

# worst (oil)
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.mars.of.G.binary[,2], col="dodgerblue3")



## NPR 1000 analytical : Gradient of G ----
time.D.npr.G.1000.analytical.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.G.1000.analytical.vec[i] <- system.time(D.npr.G.1000.analytical <- derivative.npreg.anal(X, npr.G.1000))[3]
}
time.D.npr.G.1000.analytical <- mean(time.D.npr.G.1000.analytical.vec)


corr.D.npr.G.1000.g <- corr(num.deriv.sfm$dGdxg, D.npr.G.1000.analytical[1,])
corr.D.npr.G.1000.o <- corr(num.deriv.sfm$dGdxo, D.npr.G.1000.analytical[2,])
corr.D.npr.G.1000.c <- corr(num.deriv.sfm$dGdxc, D.npr.G.1000.analytical[3,])
corr.D.npr.G.1000 <- mean(corr.D.npr.G.1000.g, corr.D.npr.G.1000.o, corr.D.npr.G.1000.c)

mse.D.npr.G.1000.g <- mean((num.deriv.sfm$dGdxg - D.npr.G.1000.analytical[1,])^2)
mse.D.npr.G.1000.o <- mean((num.deriv.sfm$dGdxo - D.npr.G.1000.analytical[2,])^2)
mse.D.npr.G.1000.c <- mean((num.deriv.sfm$dGdxc - D.npr.G.1000.analytical[3,])^2)
mse.D.npr.G.1000 <- mean(mse.D.npr.G.1000.g, mse.D.npr.G.1000.o, mse.D.npr.G.1000.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.G.1000.analytical[2,], col="dodgerblue3")



## NPR 1000 K means : Gradient of G ----
km.classif.G.1000 <- kmeans(x=npr.G.1000$eval, centers=10)
X.cl.G.1000 <- clusterX(npr.G.1000, km.classif.G.1000)

time.D.npr.G.1000.kmeans.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.G.1000.kmeans.vec[i] <- system.time(D.npr.G.1000.kmeans <- derivative.npreg.mult.data.points(X, X.cl.G.1000, s = 5, km.classif.G.1000, npr.G.1000))[3]
}
time.D.npr.G.1000.kmeans <- mean(time.D.npr.G.1000.kmeans.vec)


corr.D.npr.G.1000.kmeans.g <- corr(num.deriv.sfm$dGdxg, D.npr.G.1000.kmeans[1,])
corr.D.npr.G.1000.kmeans.o <- corr(num.deriv.sfm$dGdxo, D.npr.G.1000.kmeans[2,])
corr.D.npr.G.1000.kmeans.c <- corr(num.deriv.sfm$dGdxc, D.npr.G.1000.kmeans[3,])
corr.D.npr.G.1000.kmeans <- mean(corr.D.npr.G.1000.kmeans.g, corr.D.npr.G.1000.kmeans.o, corr.D.npr.G.1000.kmeans.c)

mse.D.npr.G.1000.kmeans.g <- mean((num.deriv.sfm$dGdxg - D.npr.G.1000.kmeans[1,])^2)
mse.D.npr.G.1000.kmeans.o <- mean((num.deriv.sfm$dGdxo - D.npr.G.1000.kmeans[2,])^2)
mse.D.npr.G.1000.kmeans.c <- mean((num.deriv.sfm$dGdxc - D.npr.G.1000.kmeans[3,])^2)
mse.D.npr.G.1000.kmeans <- mean(mse.D.npr.G.1000.kmeans.g, mse.D.npr.G.1000.kmeans.o, mse.D.npr.G.1000.kmeans.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.G.1000.kmeans[2,], col="dodgerblue3")





## NPR 4000 analytical : Gradient of G ----
time.D.npr.G.4000.analytical.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.G.4000.analytical.vec[i] <- system.time(D.npr.G.4000.analytical <- derivative.npreg.anal(X, npr.G.4000))[3]
}
time.D.npr.G.4000.analytical <- mean(time.D.npr.G.4000.analytical.vec)


corr.D.npr.G.4000.g <- corr(num.deriv.sfm$dGdxg, D.npr.G.4000.analytical[1,])
corr.D.npr.G.4000.o <- corr(num.deriv.sfm$dGdxo, D.npr.G.4000.analytical[2,])
corr.D.npr.G.4000.c <- corr(num.deriv.sfm$dGdxc, D.npr.G.4000.analytical[3,])
corr.D.npr.G.4000 <- mean(corr.D.npr.G.4000.g, corr.D.npr.G.4000.o, corr.D.npr.G.4000.c)

mse.D.npr.G.4000.g <- mean((num.deriv.sfm$dGdxg - D.npr.G.4000.analytical[1,])^2)
mse.D.npr.G.4000.o <- mean((num.deriv.sfm$dGdxo - D.npr.G.4000.analytical[2,])^2)
mse.D.npr.G.4000.c <- mean((num.deriv.sfm$dGdxc - D.npr.G.4000.analytical[3,])^2)
mse.D.npr.G.4000 <- mean(mse.D.npr.G.4000.g, mse.D.npr.G.4000.o, mse.D.npr.G.4000.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.G.4000.analytical[2,], col="dodgerblue3")



## NPR 4000 K means : Gradient of G ----
km.classif.G.4000 <- kmeans(x=npr.G.4000$eval, centers=10)
X.cl.G.4000 <- clusterX(npr.G.4000, km.classif.G.4000)

time.D.npr.G.4000.kmeans.vec <- rep(NaN, Average)
for (i in 1:Average){
     time.D.npr.G.4000.kmeans.vec[i] <- system.time(D.npr.G.4000.kmeans <- derivative.npreg.mult.data.points(X, X.cl.G.4000, s = 5, km.classif.G.4000, npr.G.4000))[3]
}
time.D.npr.G.4000.kmeans <- mean(time.D.npr.G.4000.kmeans.vec)


corr.D.npr.G.4000.kmeans.g <- corr(num.deriv.sfm$dGdxg, D.npr.G.4000.kmeans[1,])
corr.D.npr.G.4000.kmeans.o <- corr(num.deriv.sfm$dGdxo, D.npr.G.4000.kmeans[2,])
corr.D.npr.G.4000.kmeans.c <- corr(num.deriv.sfm$dGdxc, D.npr.G.4000.kmeans[3,])
corr.D.npr.G.4000.kmeans <- mean(corr.D.npr.G.4000.kmeans.g, corr.D.npr.G.4000.kmeans.o, corr.D.npr.G.4000.kmeans.c)

mse.D.npr.G.4000.kmeans.g <- mean((num.deriv.sfm$dGdxg - D.npr.G.4000.kmeans[1,])^2)
mse.D.npr.G.4000.kmeans.o <- mean((num.deriv.sfm$dGdxo - D.npr.G.4000.kmeans[2,])^2)
mse.D.npr.G.4000.kmeans.c <- mean((num.deriv.sfm$dGdxc - D.npr.G.4000.kmeans[3,])^2)
mse.D.npr.G.4000.kmeans <- mean(mse.D.npr.G.4000.kmeans.g, mse.D.npr.G.4000.kmeans.o, mse.D.npr.G.4000.kmeans.c) 

# oil
# plot(X[,2], num.deriv.sfm[,2], col="brown2")
# points(X[,2], D.npr.G.4000.kmeans[2,], col="dodgerblue3")


## Table ----
mse.vec <- c(mse.D.mars.G, mse.D.mars.bf.G.binary, mse.D.mars.of.G.binary, 
             mse.D.npr.G.1000, mse.D.npr.G.1000.kmeans, 
             mse.D.npr.G.4000, mse.D.npr.G.4000.kmeans)
corr.vec <- c(corr.D.mars.G, corr.D.mars.bf.G.binary, corr.D.mars.of.G.binary, 
              corr.D.npr.G.1000, corr.D.npr.G.1000.kmeans, 
              corr.D.npr.G.4000, corr.D.npr.G.4000.kmeans)
time.vec <- c(time.D.G.mars, time.D.mars.bf.G.binary, time.D.mars.of.G.binary, 
              time.D.npr.G.1000.analytical, time.D.npr.G.1000.kmeans, 
              time.D.npr.G.4000.analytical, time.D.npr.G.4000.kmeans)
Table.D.G <- as.data.frame(matrix(cbind(corr.vec, mse.vec, time.vec), ncol=3), row.names =  c("mars", "mars bf",  "mars of", "np 1000", "np 1000 kmeans", "np 1000", "np 4000 kmeans") )
colnames(Table.D.G) <- c("corr", "mse", "time [sec]")
round(Table.D.G, 5)


# Summary ----
round(Table.R0, 6)
round(Table.G, 6)
round(Table.D.R0, 6)
round(Table.D.G, 6)

Table.h <- cbind(round((Table.R0[,1]+Table.G[,1])/2, 4), round((Table.R0[,2]+Table.G[,2])/2, 7), round(Table.R0[,3]+Table.G[,3], 3))
rownames(Table.h) <- c("mars.bf", "mars.of", "npr.1000", "npr.4000")
colnames(Table.h) <- c("corr", "mse", "time [sec]")
Table.h


Table.grad.h <- cbind((Table.D.R0[,1:2] + Table.D.G[,1:2])/2, (Table.D.R0[,3] + Table.D.G[,3]))
colnames(Table.grad.h) <- c("corr", "mse", "time [sec]")
Table.grad.h

print(paste("Computation time using sim_forward model : ", sfm.time))

print(paste("Computation time using sim_forward model / copmutation time using MARS : ", sfm.time/(Table.R0[1,3]+Table.G[1,3])))



# Omega tile 

N <- dim(df.test)[1]
sample.var.R0 <- sum((mars.R0.bf.pred - sfm.R0)^2)/N
sample.var.G <- sum((mars.G.bf.pred - sfm.G)^2)/N
sample.cov.R0 <- sum((mars.R0.bf.pred-sfm.R0)*(mars.G.bf.pred-sfm.G))/N



# Select more spread points?  
npr.mod.R0 <- fit.NP.R0(100)
npr.mod.G <- fit.NP.G(100)
NPR.100.time.vec <- rep(NaN, Average)
for (i in 1:Average){
     NPR.100.time.vec[i] <- system.time(npr.R0.pred <- predict(object= npr.mod.R0, newdata = df.test))[3] + system.time(npr.G.pred <- predict(object= npr.mod.G, newdata = df.test))[3]
}
NPR.100.time <- mean(NPR.100.time.vec)

print(paste("Computation time of predicting with NPR model 'trained' on 100 data points: ", NPR.100.time))









