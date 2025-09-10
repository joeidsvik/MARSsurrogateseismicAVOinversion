

# This file contains functions to visualize the MCMC results and data for the Alvheim field. 
# Have not edited this much since specizalization project
# will comment on every function eventually 

#source("main_fm.R")
#library(LaplacesDemon)

library(fields)
library(latex2exp)
library(base)

source("Alvheim_transformations.R")
source("Alvheim_data.R")



result <- function(d.list, burn.in=1000){
  k <- 1
  M <- length(d.list)
  res <- data.frame(matrix(nrow = M, ncol = 2))
  d.names <- rep(NaN, M)
  for (d.name in sort(names(d.list))){
    d <- d.list[[d.name]]
    m <- length(d$gas$X1)
    a <- sum(d$accepted[burn.in:m, ])
    ess.gas <- ESS.average(d$gas, burn.in = burn.in)
    ess.oil <- ESS.average(d$oil, burn.in = burn.in)
    ess.clay <- ESS.average(d$clay, burn.in = burn.in)
    ess. <- mean(ess.gas, ess.oil, ess.clay)
    res[k,] <- c(round(100*a/(m-burn.in), 2), round(ess., 2))
    d.names[k] <- d.name
    k <- k+1
  }
  colnames(res) <- c("Acceptance rate %", "Effective sample size")
  rownames(res) <- d.names
  return(res)
}


acceptance.rate <- function(d){
  return(sum(d))
}

ESS.average <- function(d, burn.in=1000){
  m <- length(d$X1)
  n <- length(d[1,])
  #ess.d <- 0
  ess.vec <- rep(0,n)
  for (i in  1:n){
    #ess.d <- ess.d +  ESS(d[burn.in:m,i])
    ess.vec[i] <- ESS(d[burn.in:m,i])
  }
  return(mean(ess.vec))
}


# Plot travel times for selected xline and inline ----

  
plot.traveltimes <- function(inline.idx, xline.idx){
  gas <- "brown1"
  oil <- "darkorchid1"
  #well <- c(1026, 4924)
  convrt <- 1.055
  tt.area <- traveltimes[xline.idx, inline.idx]
  image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
             tt.area*convrt, 
             xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'),legend.args=list( text="Depth[m]"), 
             col = rev(viridis(256)), asp = 1, axes = F, cex.lab=1.3, cex.axis=1.2, 
             legend.width = 1.5, legend.shrink=1, axis.args=list(cex.axis=1.2))
  axis(1)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs))
  size <- 1
  rect(xleft=716-size, ybottom=-4914-size, xright=716+size, ytop=-4914+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 1
  rect(xleft=484-size, ybottom=-4798-size, xright=484+size, ytop=-4798+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 2
  rect(xleft=1300-size, ybottom=-5150-size, xright=1300+size, ytop=-5150+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 3
  rect(xleft=1026-size, ybottom=-4924-size, xright=1026+size, ytop=-4924+ size, 
       lwd=10, density=NA, border=TRUE, col=oil) # well 4
  legend("topright", legend = c("gas", "oil"), col=c(gas, oil), pch = c(19) )
}


plot.depth.with.area <- function(inline.idx, xline.idx){
  gas <- "brown1"
  oil <- "darkorchid1"
  xl <- inline[min(xline.idx), min(inline.idx)]
  xr <- inline[min(xline.idx), max(inline.idx)]
  yb <- -xline[min(xline.idx), max(inline.idx)]
  yt <- -xline[max(xline.idx), max(inline.idx)]
  
  inline.idx <- c(1:248)
  xline.idx <- c(1:178)
  
  #well <- c(1026, 4924)
  convrt <- 1.055
  tt.area <- traveltimes[xline.idx, inline.idx]
  depth <- tt.area*convrt
  zlim <- round(seq(min(depth), max(depth), 
              (max(depth)-min(depth))/5), 0)
  image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
             depth, 
             xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'),legend.args=list( text="Depth[m]"), 
             col = rev(viridis(256)), asp = 1, axes = F, cex.lab=1.3, cex.axis=1.2, 
             legend.width = 1.5, legend.shrink=1, 
             axis.args=list(at=zlim, labels=zlim, cex.axis=1.2))
  axis(1)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs))
  rect(xleft=xl, 
       ybottom=yb,
       xright=xr, 
       ytop=yt) # area
  size <- 2
  rect(xleft=716-size, ybottom=-4914-size, xright=716+size, ytop=-4914+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 1
  rect(xleft=484-size, ybottom=-4798-size, xright=484+size, ytop=-4798+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 2
  rect(xleft=1300-size, ybottom=-5150-size, xright=1300+size, ytop=-5150+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 3
  rect(xleft=1026-size, ybottom=-4924-size, xright=1026+size, ytop=-4924+ size, 
       lwd=10, density=NA, border=TRUE, col=oil) # well 4
  legend("topright", legend = c("gas", "oil"), col=c(gas, oil), pch = c(19,19) )
}



# Plot data
plot.R0 <- function(inline.idx, xline.idx){
  gas <- "brown1"
  oil <- "darkorchid1"
  xl <- inline[min(xline.idx), min(inline.idx)]
  xr <- inline[min(xline.idx), max(inline.idx)]
  yb <- -xline[min(xline.idx), max(inline.idx)]
  yt <- -xline[max(xline.idx), max(inline.idx)]
  
  inline.idx <- c(1:248)
  xline.idx <- c(1:178)
  R0 <- r0_data[xline.idx, inline.idx]
  image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
             R0, 
             xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
             col = viridis(128), asp = 1, axes = F, cex.lab=1.3, cex.axis=1.2, 
             legend.width = 1.5, legend.shrink=1, axis.args=list(cex.axis=1.2))
  axis(1)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs))
  #rect(xleft=xl, 
  #     ybottom=yb,
  #     xright=xr, 
  #     ytop=yt) # area
  size <- 1
  rect(xleft=716-size, ybottom=-4914-size, xright=716+size, ytop=-4914+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 1
  rect(xleft=484-size, ybottom=-4798-size, xright=484+size, ytop=-4798+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 2
  rect(xleft=1300-size, ybottom=-5150-size, xright=1300+size, ytop=-5150+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 3
  rect(xleft=1026-size, ybottom=-4924-size, xright=1026+size, ytop=-4924+ size, 
       lwd=10, density=NA, border=TRUE, col=oil) # well 4
  legend("topright", legend = c("gas", "oil"), col=c(gas, oil), pch = c(19,19) )
}


plot.G <- function(inline.idx, xline.idx){
  gas <- "brown1"
  oil <- "darkorchid1"
  xl <- inline[min(xline.idx), min(inline.idx)]
  xr <- inline[min(xline.idx), max(inline.idx)]
  yb <- -xline[min(xline.idx), max(inline.idx)]
  yt <- -xline[max(xline.idx), max(inline.idx)]
  
  #inline.idx <- c(1:248)
  #xline.idx <- c(1:178)
  
  G <- g_data[xline.idx, inline.idx]

  #par(mar=c(4.5,4.5,3,3))
  image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
             G, 
             xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
             col = viridis(128), asp = 1, axes = F, cex.lab=1.3, cex.axis=1.2, 
             legend.width = 1.5, legend.shrink=1, axis.args=list(cex.axis=1.2))
  axis(1)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs))
  #rect(xleft=xl, 
  #     ybottom=yb,
  #     xright=xr, 
  #     ytop=yt) # area
  size <- 1
  rect(xleft=716-size, ybottom=-4914-size, xright=716+size, ytop=-4914+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 1
  rect(xleft=484-size, ybottom=-4798-size, xright=484+size, ytop=-4798+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 2
  rect(xleft=1300-size, ybottom=-5150-size, xright=1300+size, ytop=-5150+size, 
       lwd=10, density=NA, border=TRUE, col=gas) # well 3
  rect(xleft=1026-size, ybottom=-4924-size, xright=1026+size, ytop=-4924+ size, 
       lwd=10, density=NA, border=TRUE, col=oil) # well 4
  #legend("topright", legend = c("oil"), col=c(oil), pch = c(19) )
  legend("topright", legend = c("gas", "oil"), col=c(gas, oil), pch = c(19,19) )
  
}




# Trace plots ----

trace.Gas <- function(d, filename=""){
  if (filename != ""){
    pdf(filename, width = 6, height = 6)
    par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    plot(d$gas[,3], type="l", main="", 
         xlab="iteration", ylab=TeX(r'($X_{3}$)'), cex.lab=1.4, cex.axis=1.3)
    dev.off()
  }
  else{
    #par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    plot(d$gas[,3], type="l", main="", 
         xlab="iteration", ylab=TeX(r'($X_{3}$)'), cex.lab=1.7, cex.axis=1.6)
  }
}

trace.Oil <- function(d, filename=""){
  if (filename != ""){
    pdf(filename, width = 6, height = 6)
    par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    plot(d$oil[,4], type="l", main="", 
         xlab="iteration", ylab=TeX(r'($X_{4}$)'), cex.lab=1.4, cex.axis=1.3)
    dev.off()
  }
  else{
    #par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    plot(d$oil[,4], type="l", main="", 
         xlab="iteration", ylab=TeX(r'($X_{4}$)'), cex.lab=1.7, cex.axis=1.6)
  }
}

trace.Clay <- function(d, filename=""){
  if (filename != ""){
    pdf(filename, width = 6, height = 6)
    #par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    plot(d$clay[,4], type="l", main="", 
         xlab="iteration", ylab=TeX(r'($X_{4}$)'), cex.lab=1.7, cex.axis=1.6)
    dev.off()
  }
  else{
    #par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    plot(d$clay[,4], type="l", main="", 
         xlab="iteration", ylab=TeX(r'($X_{4}$)'), cex.lab=1.7, cex.axis=1.6)
  }
}



# ACF plots ----

acf.Gas <- function(d, burn.in=1, max.lag=1000, filename=""){
  l <- length(d$gas[,1])
  m <- length(d$gas[burn.in:l, 100])
  acf100 <- acf(d$gas[burn.in:l, 100], lag.max = max.lag, plot = FALSE)$acf
  acf200 <- acf(d$gas[burn.in:l, 200], lag.max = max.lag, plot = FALSE)$acf
  acf300 <- acf(d$gas[burn.in:l, 300], lag.max = max.lag, plot = FALSE)$acf
  acf400 <- acf(d$gas[burn.in:l, 400], lag.max = max.lag, plot = FALSE)$acf
  acf500 <- acf(d$gas[burn.in:l, 500], lag.max = max.lag, plot = FALSE)$acf
  acf600 <- acf(d$gas[burn.in:l, 600], lag.max = max.lag, plot = FALSE)$acf
  acf700 <- acf(d$gas[burn.in:l, 700], lag.max = max.lag, plot = FALSE)$acf
  acf800 <- acf(d$gas[burn.in:l, 800], lag.max = max.lag, plot = FALSE)$acf
  if (filename != ""){
    pdf(filename, width = 6, height = 6)
    par(mar=c(5,5,2,1.5))
    plot(acf100, type = "l", col = "red3", ylim=c(-0.2, 1), 
         xlab="lag", ylab="ACF", 
         cex.lab=1.4, cex.axis=1.3)
    lines(acf200, col = "royalblue")
    lines(acf300, col = "hotpink")
    lines(acf400, col = "darkorange1")
    lines(acf500, col = "blueviolet")
    lines(acf600, col = "gold2")
    lines(acf700, col = "dodgerblue")
    lines(acf800, col = "springgreen4")
    abline(h=0, col="black", cex=1.2)
    abline(h=1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    abline(h=-1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    dev.off()
  }
  else {
    #par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    par(mar=c(5,5,2,1.5))
    plot(acf100, type = "l", col = "red3", ylim=c(-0.2, 1), 
         xlab="lag", ylab="ACF", 
         cex.lab=1.4, cex.axis=1.3)
    lines(acf200, col = "royalblue")
    lines(acf300, col = "hotpink")
    lines(acf400, col = "darkorange1")
    lines(acf500, col = "blueviolet")
    lines(acf600, col = "gold2")
    lines(acf700, col = "dodgerblue")
    lines(acf800, col = "springgreen4")
    abline(h=0, col="black", cex=1.2)
    abline(h=1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    abline(h=-1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
  }
}


acf.Oil <- function(d, burn.in=1, max.lag=1000, filename=""){
  l <- dim(d$gas)[1]
  m <- length(d$oil[burn.in:l, 100])
  acf100 <- acf(d$oil[burn.in:l, 100], lag.max = max.lag, plot = FALSE)$acf
  acf200 <- acf(d$oil[burn.in:l, 200], lag.max = max.lag, plot = FALSE)$acf
  acf300 <- acf(d$oil[burn.in:l, 300], lag.max = max.lag, plot = FALSE)$acf
  acf400 <- acf(d$oil[burn.in:l, 400], lag.max = max.lag, plot = FALSE)$acf
  acf500 <- acf(d$oil[burn.in:l, 500], lag.max = max.lag, plot = FALSE)$acf
  acf600 <- acf(d$oil[burn.in:l, 600], lag.max = max.lag, plot = FALSE)$acf
  acf700 <- acf(d$oil[burn.in:l, 700], lag.max = max.lag, plot = FALSE)$acf
  acf800 <- acf(d$oil[burn.in:l, 800], lag.max = max.lag, plot = FALSE)$acf
  if (filename != ""){
    pdf(filename, width = 6, height = 6)
    par(mar=c(5,5,2,1.5))
    plot(acf100, type = "l", col = "red3", ylim=c(-0.1, 1), 
         xlab="lag", ylab="ACF", 
         cex.lab=1.4, cex.axis=1.3)
    lines(acf200, col = "royalblue")
    lines(acf300, col = "hotpink")
    lines(acf400, col = "darkorange1")
    lines(acf500, col = "blueviolet")
    lines(acf600, col = "gold2")
    lines(acf700, col = "dodgerblue")
    lines(acf800, col = "springgreen4")
    abline(h=0, col="black", cex=1.2)
    abline(h=1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    abline(h=-1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    dev.off()
  }
  else {
    #par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    par(mar=c(5,5,2,1.5))
    plot(acf100, type = "l", col = "red3", ylim=c(-0.1, 1), 
         xlab="lag", ylab="ACF", 
         cex.lab=1.4, cex.axis=1.3)
    lines(acf200, col = "royalblue")
    lines(acf300, col = "hotpink")
    lines(acf400, col = "darkorange1")
    lines(acf500, col = "blueviolet")
    lines(acf600, col = "gold2")
    lines(acf700, col = "dodgerblue")
    lines(acf800, col = "springgreen4")
    abline(h=0, col="black", cex=1.2)
    abline(h=1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    abline(h=-1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
  }
}


acf.Clay <- function(d, burn.in=1, max.lag=1000, filename=""){
  l <- dim(d$gas)[1]
  m <- length(d$clay[burn.in:l, 100])
  acf100 <- acf(d$clay[burn.in:l, 100], lag.max = max.lag, plot = FALSE)$acf
  acf200 <- acf(d$clay[burn.in:l, 200], lag.max = max.lag, plot = FALSE)$acf
  acf300 <- acf(d$clay[burn.in:l, 300], lag.max = max.lag, plot = FALSE)$acf
  acf400 <- acf(d$clay[burn.in:l, 400], lag.max = max.lag, plot = FALSE)$acf
  acf500 <- acf(d$clay[burn.in:l, 500], lag.max = max.lag, plot = FALSE)$acf
  acf600 <- acf(d$clay[burn.in:l, 600], lag.max = max.lag, plot = FALSE)$acf
  acf700 <- acf(d$clay[burn.in:l, 700], lag.max = max.lag, plot = FALSE)$acf
  acf800 <- acf(d$clay[burn.in:l, 800], lag.max = max.lag, plot = FALSE)$acf
  if (filename != ""){
    pdf(filename, width = 6, height = 6)
    par(mar=c(5,5,2,1.5))
    plot(acf100, type = "l", col = "red3", ylim=c(-0.1, 1), 
         xlab="lag", ylab="ACF", 
         cex.lab=1.4, cex.axis=1.3)
    lines(acf200, col = "royalblue")
    lines(acf300, col = "hotpink")
    lines(acf400, col = "darkorange1")
    lines(acf500, col = "blueviolet")
    lines(acf600, col = "gold2")
    lines(acf700, col = "dodgerblue")
    lines(acf800, col = "springgreen4")
    abline(h=0, col="black", cex=1.2)
    abline(h=1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    abline(h=-1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    dev.off()
  }
  else {
    #par(mar=c(4,5,2,1.5), oma =c(0,0,0,1))
    par(mar=c(5,5,2,1.5))
    plot(acf100, type = "l", col = "red3", ylim=c(-0.1, 1), 
         xlab="lag", ylab="ACF", 
         cex.lab=1.4, cex.axis=1.3)
    lines(acf200, col = "royalblue")
    lines(acf300, col = "hotpink")
    lines(acf400, col = "darkorange1")
    lines(acf500, col = "blueviolet")
    lines(acf600, col = "gold2")
    lines(acf700, col = "dodgerblue")
    lines(acf800, col = "springgreen4")
    abline(h=0, col="black", cex=1.2)
    abline(h=1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
    abline(h=-1.96/sqrt(m), col="steelblue4", cex=1.2, lty=2)
  }
}



# Transformation of mean samples ----
# Transformation of mean samples ----
plot.Sg <- function(d, burn.in=1, filename =""){
  l <- dim(d$gas)[1]
  r <- length(xline.idx)
  k <- length(inline.idx)
  d.mean.gas <- apply(X=as.matrix(d$gas[burn.in:l,]), MAR=2, FUN=mean)
  d.mean.oil <- apply(d$oil[burn.in:l,], 2, mean)
  d.mean.clay <- apply(d$clay[burn.in:l,], 2, mean)
  
  x.mean <- cbind(d.mean.gas, d.mean.oil, d.mean.clay)
  S.mean.g <- S.g(x.mean)
  
  if (filename!=""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    map.mean.S.g <- matrix(as.matrix(S.mean.g), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.S.g, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               zlim = c(0,1),
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
  }
  else{
    map.mean.S.g <- matrix(as.matrix(S.mean.g), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.S.g, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               zlim = c(0,1),
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6)) 
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6) 
  }
}


plot.So <- function(d, burn.in=1, filename =""){
  l <- dim(d$oil)[1]
  r <- length(xline.idx)
  k <- length(inline.idx)
  d.mean.gas <- apply(d$gas[burn.in:l,], 2, mean)
  d.mean.oil <- apply(d$oil[burn.in:l,], 2, mean)
  d.mean.clay <- apply(d$clay[burn.in:l,], 2, mean)
  
  x.mean <- cbind(d.mean.gas, d.mean.oil, d.mean.clay)
  S.mean.o <- S.o(x.mean)
  
  if (filename!=""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    map.mean.S.o <- matrix(as.matrix(S.mean.o), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.S.o, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1,
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6)) 
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
  }
  else{
    map.mean.S.o <- matrix(as.matrix(S.mean.o), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.S.o, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6)) 
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
  }
}


plot.Sb <- function(d, burn.in=1, filename =""){
  l <- dim(d$gas)[1]
  r <- length(xline.idx)
  k <- length(inline.idx)
  d.mean.gas <- apply(d$gas[burn.in:l,], 2, mean)
  d.mean.oil <- apply(d$oil[burn.in:l,], 2, mean)
  d.mean.clay <- apply(d$clay[burn.in:l,], 2, mean)
  
  x.mean <- cbind(d.mean.gas, d.mean.oil, d.mean.clay)
  S.mean.b <- S.b(x.mean)
  
  if (filename!=""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    map.mean.S.b <- matrix(as.matrix(S.mean.b), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.S.b, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
  }
  else{
    map.mean.S.b <- matrix(as.matrix(S.mean.b), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.S.b, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
  }
}


plot.Vcl <- function(d, burn.in=1, filename =""){
  l <- dim(d$gas)[1]
  r <- length(xline.idx)
  k <- length(inline.idx)
  d.mean.gas <- apply(d$gas[burn.in:l,], 2, mean)
  d.mean.oil <- apply(d$oil[burn.in:l,], 2, mean)
  d.mean.clay <- apply(d$clay[burn.in:l,], 2, mean)
  
  x.mean <- cbind(d.mean.gas, d.mean.oil, d.mean.clay)
  V.mean.cl <- V.cl(x.mean)
  
  if (filename!=""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    map.mean.V.cl <- matrix(as.matrix(V.mean.cl), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.V.cl, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.4, cex.axis=1.3, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.3))
    axis(1, cex.axis=1.3)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.3)
    dev.off()
  }
  else{
    map.mean.V.cl <- matrix(as.matrix(V.mean.cl), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               map.mean.V.cl, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
  }
}

# Plot prior mu ----

# Plot mu prior gas
#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/muPriorGas.pdf"
plot.mu.prior.gas <- function(inline.idx, xline.idx, filename = "", transformation = 0){
  tt_area <- traveltimes[xline.idx, inline.idx]
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  mu.prior.gas <- mu.prior[,1]
  if (transformation){
    mu.prior.gas <- S.g(mu.prior)
  }
  
  if ( filename != ""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    prior.gas <- matrix(as.matrix(mu.prior.gas), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.gas, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main="",  col = viridis(128), axes=F, 
               zlim = c(0,1),
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1,
               axis.args=list(cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
  }
  else{
    par(mar=c(4,4.5,2,1.5), oma =c(0,0,0,1))
    prior.gas <-  matrix(as.matrix(mu.prior.gas), nrow = r, ncol = k)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.gas, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'),
               main="",  col = viridis(128), axes=F, 
               zlim = c(0,1),
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               axis.args=list(cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
  }
}


# Plot mu prior oil
#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/muPriorOil.pdf"
plot.mu.prior.oil <- function(inline.idx, xline.idx, filename = "", transformation = 0){
  tt_area <- traveltimes[xline.idx, inline.idx]
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  mu.prior.oil <- mu.prior[,2]
  if (transformation){
    mu.prior.oil <- S.o(mu.prior)
  }
  
  if (filename != ""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    prior.oil <-  matrix(as.matrix(mu.prior.oil), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.oil, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main="",  col = viridis(128), axes=F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
  }
  else{
    par(mar=c(4,4.5,2,1.5), oma =c(0,0,0,1))
    prior.oil <-  matrix(as.matrix(mu.prior.oil), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.oil, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main="",  col = viridis(128), axes=F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1,
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
  }
}



plot.mu.prior.brine <- function(inline.idx, xline.idx, filename = ""){
  tt_area <- traveltimes[xline.idx, inline.idx]
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  mu.prior.brine <- S.b(mu.prior)
  
  if (filename != ""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    prior.brine <-  matrix(as.matrix(mu.prior.brine), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.brine, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main="",  col = viridis(128), axes=F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1,
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
  }
  else{
    par(mar=c(4,4.5,2,1.5), oma =c(0,0,0,1))
    prior.brine <-  matrix(as.matrix(mu.prior.brine), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.brine, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main="",  col = viridis(128), axes=F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1,
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
  }
}



# Plot mu prior clay
#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/muPriorClay.pdf"
plot.mu.prior.clay <- function(inline.idx, xline.idx, filename = "", transformation = 0){
  tt_area <- traveltimes[xline.idx, inline.idx]
  mu.prior <- match_traveltimes_mean_trend(tt_area)
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  mu.prior.clay <- mu.prior[,3]
  if (transformation){
    mu.prior.clay <- V.cl(mu.prior)
  }
  
  if (filename != ""){
    pdf(filename)
    par(mar=c(5, 5, 4, 2))
    prior.clay <-  matrix(as.matrix(mu.prior.clay), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.clay, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main="",  col = viridis(128), axes=F,
               zlim = c(0,1),
               cex.lab=1.4, cex.axis=1.3, legend.width = 1.5, legend.shrink=1,
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.3))
    axis(1, cex.axis=1.2)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.3)
    dev.off()
  }
  else{
    par(mar=c(4,4.5,2,1.5), oma =c(0,0,0,1))
    prior.clay <-  matrix(as.matrix(mu.prior.clay), nrow = r, ncol = k)
    breaks <- seq(0, 1,  length.out=6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               prior.clay, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main="",  col = viridis(128), axes=F, 
               cex.lab=1.4, cex.axis=1.3, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.3))
    axis(1, cex.axis=1.2)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.3)
  }
}

# Percentile plots ----

# Minas results
realis.b <- function(realis.g, realis.o){
  S.b <- matrix(data=NA, nrow=100, ncol=900)
  for (i in 1:100){
    S.b[i,] <- 1 / (1 + realis.g[[i]] + realis.o[[i]])
  }
  return(S.b)
}


# MCMC results
Sb.MCMC <- function(x.g, x.o, inline.idx=inline.idx, xline.idx=xline.idx){
  Sb <- 1 / (1 + exp(x.g) + exp(x.o))
  return(Sb)
}



percentile.brine <- function(S.b, inline.idx, xline.idx, p=80, m=100){
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  S.b.sorted <- apply(S.b, 2, sort)
  
  n = floor((p/100) * m)
  lb <- floor((m - n)/2)
  ub <- floor(m - (m-n)/2)
  
  S.b.p <- S.b.sorted[ub,] - S.b.sorted[lb,]
  S.b.p.mat <- matrix(as.matrix(S.b.p), nrow = r, ncol = k)
  return(S.b.p.mat)
}


percentile <- function(d, inline.idx, xline.idx, bi=500, p=80, filename=""){
  m <- length(d$gas[,1])
  r <- length(xline.idx)
  k <- length(inline.idx)
  
  x.g <- d$gas[bi:m, ]
  x.o <- d$oil[bi:m, ]
  x.c <- d$clay[bi:m, ]
  Sb <- Sb.MCMC(x.g, x.o)
  
  m <- length(x.g[,1])
  
  x.g.sorted <- apply(x.g, 2, sort)
  x.o.sorted <- apply(x.o, 2, sort)
  x.c.sorted <- apply(x.c, 2, sort)
  Sb.sorted <- apply(Sb, 2, sort)
  
  n = floor((p/100) * m)
  lb <- floor((m - n)/2)
  ub <- floor(m - (m-n)/2)
  
  x.sorted.ub <- cbind(x.g.sorted[ub,], x.o.sorted[ub,], x.c.sorted[ub,])
  x.sorted.lb <- cbind(x.g.sorted[lb,], x.o.sorted[lb,], x.c.sorted[lb,])
  
  s.g.p <- abs(S.g(x.sorted.ub) - S.g(x.sorted.lb))
  s.o.p <- abs(S.o(x.sorted.ub) - S.o(x.sorted.lb))
  s.b.p <- abs(Sb.sorted[ub,] - Sb.sorted[lb,])
  s.c.p <- abs(V.cl(x.sorted.ub) - V.cl(x.sorted.lb))
  
  S.g.matrix <- matrix(as.matrix(s.g.p), nrow = r, ncol = k)
  S.o.matrix <- matrix(as.matrix(s.o.p), nrow = r, ncol = k)
  S.b.matrix <- matrix(as.matrix(s.b.p), nrow = r, ncol = k)
  V.c.matrix <- matrix(as.matrix(s.c.p), nrow = r, ncol = k)
  
  breaks <- seq(0, 1,  length.out=6)
  if (filename!=""){
    
    pdf(paste(filename, "_gas.pdf"), width = 6, height = 6)
    par(mar=c(5, 5, 4, 2))
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx],
               S.g.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
    
    pdf(paste(filename, "_oil.pdf"), width = 6, height = 6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               S.o.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
    
    pdf(paste(filename, "_brine.pdf"), width = 6, height = 6)
    par(mar=c(5, 5, 4, 2))
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               S.b.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    dev.off()
    
    pdf(paste(filename, "_clay.pdf"), width = 6, height = 6)
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               V.c.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "", col = viridis(128), axes = F, 
               cex.lab=1.4, cex.axis=1.3, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.3)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.3)
    dev.off()
  }
  else{
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               S.g.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "Gas", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               S.o.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "Oil", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               S.b.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "Brine", col = viridis(128), axes = F, 
               cex.lab=1.7, cex.axis=1.6, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.6)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.6)
    
    image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
               V.c.matrix, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
               main = "Clay", col = viridis(128), axes = F, 
               cex.lab=1.4, cex.axis=1.3, legend.width = 1.5, legend.shrink=1, 
               zlim = c(0,1),
               axis.args=list(at=breaks, labels=breaks, cex.axis=1.6))
    axis(1, cex.axis=1.3)
    box()
    ylabs <- axis(2, labels = FALSE, tick = FALSE)
    axis(2, at = ylabs, labels=-(ylabs), cex.axis=1.3)
  }
}

mat <- matrix(c(1,2,3,3,2,1,2,3,1), nrow = 3)
mat.sorted <- apply(mat, 2, sort)
mat
mat.sorted


# Plot Minas results ----

plot.func <- function(thing, inline.idx, xline.idx){
  set.cex.axis=1.3
  set.cex.lab=1.4
  breaks <- seq(0, 1,  length.out=6)
  image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
             thing, xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'), 
             main = "", col = viridis(128), axes = F,
             cex.lab=set.cex.lab, cex.axis=set.cex.axis, 
             legend.width = 1.5, legend.shrink=1, 
             zlim = c(0,1),
             axis.args=list(at=breaks, labels= breaks, cex.axis=set.cex.axis))
  axis(1, cex.axis=set.cex.axis)
  box()
  ylabs <- axis(2, labels = FALSE, tick = FALSE)
  axis(2, at = ylabs, labels=-(ylabs), cex.axis=set.cex.axis)
}

