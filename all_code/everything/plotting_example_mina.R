# Example figure plotting with well-locations
library(viridis)
fig_color <- viridis(256) 


#### Depth plot with well locations ####
png("minor_revision_figures/depth.png",units="px", width = 2086, height= 1500, pointsize = 40)

par(mfrow=c(1,1), mar=c(5,4,4,6))
image.plot(x=x_grid_full, y = -y_grid_full, depth_full, xlab="inline", ylab="crossline", main="",legend.args=list( text="Depth (m)"), col = rev(viridis(256)),xaxt = "n", yaxt = "n", yaxs = "i")
axis(1)
box()
ylabs <- axis(2, labels = FALSE, tick = FALSE)
axis(2, at = ylabs, labels=-(ylabs), las = 2)

# Kameleon - 24/6-2
points(x = 716, y = -4914, bg ="red", col = "black", pch = 22, cex = 1.7)
# BOA 24/6-4
points(x = 484, y = -4798, bg = "red", col = "black", pch = 22, cex = 1.7)
# Kneler 25/4-7
points(x = 1026, y = -4924, bg = "green", col = "black", pch = 24, cex = 1.7)
# Gekko (part of) 25/4-8
points(x = 1300, y = -5150, bg = "red", col = "black", pch = 22, cex = 1.7)  
# 25/4-J1-AH Inline 1115 and Xline 4610
#points(x = 1115, y = 4610, col = "green", pch = 8 )  

legend("topright", legend=c("Gas","Oil"),
       pt.bg=c("red","green"), col = c("black", "black"),pch=c(22,24), cex = 1.7)
dev.off()