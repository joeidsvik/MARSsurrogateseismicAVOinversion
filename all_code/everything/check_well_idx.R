
source("Alvheim_data.R")

# small area
#inline.idx = seq(156, 215, 2)
#xline.idx = seq(71, 130, 2)

# total area
inline.idx = seq(1, 248, 1)
xline.idx = seq(1, 178, 1)




# Wells 
r <- length(xline.idx)
k <- length(inline.idx)

# well 1 (gas) (484,4798) = (inline, xline)
j.well.1 <- which(inline[1,] == 484)
i.well.1 <- c(which(xline[,1] == 4796), which(xline[,1] == 4800))
j.1 <- which(inline.idx == j.well.1)
i.1 <- c(which(xline.idx == i.well.1[1]), which(xline.idx == i.well.1[2]))
idx.well.1 <- c(((j.1-1)*r+i.1),((j.1)*r+i.1))

# well 2 (gas) (716,4914)
j.well.2 <- which(inline[1,] == 716)
i.well.2 <- c(which(xline[,1] == 4912), which(xline[,1] == 4916))
j.2 <- which(inline.idx == j.well.2)
i.2 <- c(which(xline.idx == i.well.2[1]), which(xline.idx == i.well.2[2]))
idx.well.2 <- c(((j.2-1)*r+i.2),((j.2)*r+i.2))

# well 3 (gas) (1300,5150)
j.well.3 <- which(inline[1,] == 1296) # highest inline
i.well.3 <- c(which(xline[,1] == 5148), which(xline[,1] == 5152))
j.3 <- which(inline.idx == j.well.3)
i.3 <- c(which(xline.idx == i.well.3[1]), which(xline.idx == i.well.3[2]))
idx.well.3 <- c(((j.3-1)*r+i.3),((j.3)*r+i.3))[1:2] 

# well 4 (oil) (1026,4924)
j.well.4 <- c(which(inline[1,] == 1024), which(inline[1,] == 1028))
i.well.4 <- which(xline[,1] == 4924)
j.4 <- c(which(inline.idx == j.well.4[1]), which(inline.idx == j.well.4[2]))
i.4 <- which(xline.idx == i.well.4)
idx.well.4 <- c(((j.4-1)*r+i.4),((j.4)*r+i.4))


gas <- "brown1"
oil <- "darkorchid1"

tt.area <- traveltimes[xline.idx, inline.idx]
depth <- match_traveltimes_mean_trend(tt.area)[,4]

depth[idx.well.1] <- 1
depth[idx.well.2] <- 1
depth[idx.well.3] <- 1
depth[idx.well.4] <- 1

depth.mat <- matrix(depth, ncol=k, nrow=r)

#par(mar=c(4.5,4.5,3,3))
image.plot(x=inline[xline.idx,inline.idx], y = -xline[xline.idx,inline.idx], 
           depth.mat, 
           xlab=TeX(r'($inline$)'), ylab=TeX(r'($xline$)'),legend.args=list( text="Depth[m]"), 
           col = rev(viridis(256)), asp = 1, axes = F, cex.lab=1.3, cex.axis=1.2, 
           legend.width = 1.5, legend.shrink=1, axis.args=list(cex.axis=1.2))
axis(1)
box()
ylabs <- axis(2, labels = FALSE, tick = FALSE)
axis(2, at = ylabs, labels=-(ylabs))
size <- 1
rect(xleft=484-size, ybottom=-4798-size, xright=484+size, ytop=-4798+size, 
     lwd=10, density=NA, border=TRUE, col=gas) # well 2
rect(xleft=716-size, ybottom=-4914-size, xright=716+size, ytop=-4914+size, 
     lwd=10, density=NA, border=TRUE, col=gas) # well 1
rect(xleft=1300-size, ybottom=-5150-size, xright=1300+size, ytop=-5150+size, 
     lwd=10, density=NA, border=TRUE, col=gas) # well 3
rect(xleft=1026-size, ybottom=-4924-size, xright=1026+size, ytop=-4924+ size, 
     lwd=10, density=NA, border=TRUE, col=oil) # well 4
legend("topright", legend = c("gas", "oil"), col=c(gas, oil), pch = c(19) )














