# Plot Alvheim on total area

source("Alvheim_data.R")
source("Alvheim_plot_aux.R")

library(LaplacesDemon)


inmin <- 1
inmax <- 248
xmin <- 1
xmax <- 178
by <- 1

inline.idx = seq(inmin, inmax, by)
xline.idx = seq(xmin, xmax, by)

N <- length(inline.idx) * length(xline.idx)

m <- 10000
save.every <- 10 


# Results ----
d1$accepted/m

d2$accepted/m

d3$accepted/m

d4$accepted/m

l <- 210
ESS.d3 <- rep(NaN, N)
for(l in 1:N){
  ESS.d3[l] <- mean(c(ESS(d3$gas[500:1000,l]), ESS(d3$oil[500:1000,l]), ESS(d3$clay[500:1000,l])))
}
mean(ESS.d3)
max(ESS.d3)


ESS(d1$gas[500:1000,l])
ESS(d2$gas[5000:10000,l])
ESS(d3$gas[500:1000,l])

ESS(d4$gas[5000:10000,l])

# Visualization of MCMC runs ----

par(mfrow = c(2, 2))
plot(d1$oil[,200], type = "l")
plot(d1$oil[,400], type = "l")
plot(d1$oil[,600], type = "l")
plot(d1$oil[,800], type = "l")

plot(d2$gas[,200], type = "l")
plot(d2$gas[,400], type = "l")
plot(d2$gas[,600], type = "l")
plot(d2$gas[,800], type = "l")

plot(d3$gas[,200], type = "l")
plot(d3$gas[,400], type = "l")
plot(d3$gas[,600], type = "l")
plot(d3$gas[,800], type = "l")

plot(d4$gas[,200], type = "l")
plot(d4$gas[,400], type = "l")
plot(d4$gas[,600], type = "l")
plot(d4$gas[,800], type = "l")




b <- 500

plot.Sg(d1, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q1_gas.pdf")  
plot.So(d1, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q1_oil.pdf")  
plot.Sb(d1, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q1_brine.pdf")  
plot.Vcl(d1, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q1_clay.pdf")  

plot.Sg(d2, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q2_gas.pdf")  
plot.So(d2, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q2_oil.pdf")   
plot.Sb(d2, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q2_brine.pdf") 
plot.Vcl(d2, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q2_clay.pdf") 

plot.Sg(d3, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q3_230_gas.pdf") 
plot.So(d3, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q3_230_oil.pdf")  
plot.Sb(d3, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q3_230_brine.pdf") 
plot.Vcl(d3, burn.in = b) #, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q3_230_clay.pdf") 

plot.Sg(d4, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q4_gas.pdf") 
plot.So(d4, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q4_oil.pdf") 
plot.Sb(d4, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q4_brine.pdf") 
plot.Vcl(d4, burn.in = b) # , filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/q4_clay.pdf") 





par(mfrow = c(3, 4))
filename <- "/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/total_area/uncertainty_q3_230"
percentile(d=d3, inline.idx, xline.idx, bi=500, p=80, filename=filename)














