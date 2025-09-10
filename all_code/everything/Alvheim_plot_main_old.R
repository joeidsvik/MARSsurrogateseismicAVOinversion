

# In this file one can create plots for MCMC runs and other plots related to Alvheim. 
# Will cleen this up eventually :) 


source("Alvheim_transformations.R")
source("Alvheim_plot_aux.R")
source("Alvheim_data.R")

# not sure if I need these packages
library(Matrix)
library(nlme)
library(mcmcse)



#pdf("/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/depth.pdf", width=6, height=6)
plot.traveltimes(inline.idx = c(1:248), xline.idx=c(1:178))
#dev.off()

#pdf("/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/depth_w_area.pdf", width=6, height=6)
plot.depth.with.area(inline.idx, xline.idx)
#dev.off()


#pdf("/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/R0.pdf", width=6, height=6)
plot.R0(inline.idx = c(1:248), xline.idx=c(1:178))
#dev.off()


#pdf("/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/G.pdf", width=6, height=6)
plot.G(inline.idx = c(1:248), xline.idx=c(1:178))
#dev.off()


# load MCMC runs ----



# Convergence analysis ----

# Trace plots ----

par(mfrow=c(1,2), mar=c(4.5,4.5,3,3))

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/traceGas.pdf"
trace.Gas(d1, filename = "")
trace.Gas(d2, filename = "")
trace.Gas(d2_mars, filename = "")
trace.Gas(d2_mars_tilde, filename = "")
trace.Gas(d3, filename = "")

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/traceOil.pdf"
trace.Oil(d1, filename = "")
trace.Oil(d2, filename = "")
trace.Oil(d3, filename = "")

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/traceClay.pdf"
trace.Clay(d1, filename = "")
trace.Clay(d2, filename = "")
trace.Clay(d3, filename = "")


# Remove burn in
bi <- 5000

l <- length(d2$accepted[,1])
(sum(d1$accepted[bi:l,])/(l-bi))*100
(sum(d2$accepted[bi:l,])/(l-bi))*100
(sum(d3$accepted[bi:l,])/(l-bi))*100



# ACF plots ----
source("Alvheim_visualize.R")
#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/acfGas.pdf"
acf.Gas(d1, burn.in=bi, max.lag=1000, filename= "")
acf.Gas(d2, burn.in=bi, max.lag=200, filename= "")
acf.Gas(d2_mars, burn.in=bi, max.lag=200, filename= "")
acf.Gas(d2_mars_tilde, burn.in=bi, max.lag=200, filename= "")
acf.Gas(d3, burn.in=bi, max.lag=1000, filename= "")

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/acfOil.pdf"
acf.Oil(d1, burn.in=bi, max.lag=1000, filename= "")
acf.Oil(d2, burn.in=bi, max.lag=600, filename= "")
acf.Oil(d3, burn.in=bi, max.lag=1000, filename= "")

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/acfClay.pdf"
acf.Clay(d1, burn.in=bi, max.lag=1000, filename="")
acf.Clay(d2, burn.in=bi, max.lag=1000, filename="")
acf.Clay(d3, burn.in=bi, max.lag=1000, filename="")


# Transformation of mean plot ----
par(mfrow = c(1, 3))
filename <- "/Users/karenauestad/Dokumenter/master/masterpresentasjon.pdf"
#plot.Sg(d1, burn.in=bi, filename ="")
plot.Sg(d2, burn.in=bi, filename ="")
plot.Sg(d2_mars, burn.in=bi, filename ="")
plot.Sg(d2_mars_tilde, burn.in=bi, filename ="")
#plot.Sg(d3, burn.in=bi, filename ="")


filename <- "/Users/karenauestad/Dokumenter/master/masterpresentasjon.pdf"
#plot.So(d1, burn.in=bi, filename ="")
plot.So(d2, burn.in=bi, filename ="")
plot.So(d2_mars, burn.in=bi, filename ="")
plot.So(d2_mars_tilde, burn.in=bi, filename ="")
#plot.So(d3, burn.in=bi, filename ="")


#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/SmeanBrine.pdf"
#plot.Sb(d1, burn.in=bi, filename = "")
plot.Sb(d2, burn.in=bi, filename = "")
plot.Sb(d2_mars, burn.in=bi, filename = "")
plot.Sb(d2_mars_tilde, burn.in=bi, filename = "")
#plot.Sb(d3, burn.in=bi, filename = "")
#

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/VmeanClay.pdf"
#plot.Vcl(d1, burn.in=bi, filename = "")
plot.Vcl(d2, burn.in=bi, filename = "")
plot.Vcl(d2_mars, burn.in=bi, filename = "")
plot.Vcl(d2_mars_tilde, burn.in=bi, filename = "")
#plot.Vcl(d3, burn.in=bi, filename = "")



# Plot mu_prior ----
filename <- "/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/muPriorGas.pdf"
plot.mu.prior.gas(inline.idx, xline.idx, filename = "", transformation = 1)


filename <- "/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/muPriorOil.pdf"
plot.mu.prior.oil(inline.idx, xline.idx, filename = "", transformation = 1)


filename <- "/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/muPriorBrine.pdf"
plot.mu.prior.brine(inline.idx, xline.idx, filename = "")


filename <- "//Users/karenauestad/Dokumenter/master/plots_master/Alvheim/muPriorClay.pdf"
plot.mu.prior.clay(inline.idx, xline.idx, filename = "", transformation = 1)


par(mfrow = c(3, 4))
#filename <- "/Users/karenauestad/Dokumenter/master/masterpresentasjon"
percentile(d=d2, inline.idx, xline.idx, bi=5000, p=80, filename="")
percentile(d=d2_mars, inline.idx, xline.idx, bi=5000, p=80, filename="")
percentile(d=d2_mars_tilde, inline.idx, xline.idx, bi=5000, p=80, filename="")


# Results from Mina ----

realis.g <- readRDS("results_karen/realis_sat_g_ver1.rds") 
realis.o <- readRDS("results_karen/realis_sat_o_ver1.rds")
realis.b <- realis.b(realis.g, realis.o)
realist.c <- readRDS("results_karen/realis_clay_ver1.rds") 


percentile.sat.g <- readRDS("results_karen/uncert_sat_g_ver1.rds")
percentile.sat.o <- readRDS("results_karen/uncert_sat_o_ver1.rds")
percentile.sat.b <- percentile.brine(realis.b, inline.idx, xline.idx)
percentile.sat.c <- readRDS("results_karen/uncert_vclay_ver1.rds")



pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanGasUncertainty.pdf", width=6, height=6)
plot.func(percentile.sat.g, inline.idx, xline.idx)
dev.off()

pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanOilUncertainty.pdf", width=6, height=6)
plot.func(percentile.sat.o, inline.idx, xline.idx)
dev.off()

pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanBrineUncertainty.pdf", width=6, height=6)
plot.func(percentile.sat.b, inline.idx, xline.idx)
dev.off()

pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanClayUncertainty.pdf", width=6, height=6)
plot.func(percentile.sat.c, inline.idx, xline.idx)
dev.off()




sat.g <- readRDS("results_karen/mean_sat_g_ver1.rds") 
sat.o <- readRDS("results_karen/mean_sat_o_ver1.rds") 
sat.b <- 1 - sat.g - sat.o 
sat.c <- readRDS("results_karen/mean_clay_ver1.rds") 


#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanGas.pdf", width=6, height=6)
plot.func(sat.g, inline.idx, xline.idx)
#dev.off()

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanOil.pdf", width=6, height=6)
plot.func(sat.o, inline.idx, xline.idx)
#dev.off()

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanBrine.pdf", width=6, height=6)
plot.func(sat.b, inline.idx, xline.idx)
#dev.off()

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanClay.pdf", width=6, height=6)
plot.func(sat.c, inline.idx, xline.idx)
#dev.off()















