
source("Alvheim_prior.R")
source("Alvheim_likelihood.R")
source("Alvheim_proposal.R")
source("Alvheim_transformations.R")
source("Alvheim_visualize.R")

library(Matrix)
library(nlme)
library(mcmcse)


inmin <- 171 - 15
inmax <- 200 + 15
xmin <- 86 - 15
xmax <- 115 + 15
by <- 2

inline.idx = seq(inmin, inmax, 2)
xline.idx = seq(xmin, xmax, 2)

inline.idx = seq(156, 215, 2)
xline.idx = seq(71, 130, 2)


#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/depthArea.pdf", width=6, height=6)
plot.traveltimes(inline.idx, xline.idx)
#dev.off()

#pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/depth.pdf", width=6, height=6)
plot.depth.with.area(inline.idx, xline.idx)
#dev.off()


#load("Dq1.rda")
#result(Dq1, burn.in = 500)

#load("Dq2.rda")
#result(Dq2, burn.in = 500)

#load("Dq3.rda")
#result(Dq3, burn.in = 500)

#d1 <- Dq1$d005
#d2 <- Dq2$d035
#d3 <- Dq3$d02
#par(mfrow = c(1, 2))

load("~/Desktop/var2023/prosjektoppgave/Alvheim/forwardmodel_AF/dq2_035.rda")
Dq2.035 <- list("dq2.035" = dq2.035)
result(Dq2.035, burn.in = 500)

d2 <- dq2.035
rm(Dq2.035)
rm(dq2.035)


# Convergence analysis ----

# Trace plots ----

par(mfrow=c(1,2), mar=c(4.5,4.5,3,3))

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/traceGas.pdf"
trace.Gas(d1, filename = "")
trace.Gas(d2, filename = "")
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
bi <- 500

l <- length(d2$accepted[,1])
(sum(d1$accepted[bi:l,])/(l-bi))*100
(sum(d2$accepted[bi:l,])/(l-bi))*100
(sum(d3$accepted[bi:l,])/(l-bi))*100



# ACF plots ----
source("Alvheim_visualize.R")
#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/acfGas.pdf"
acf.Gas(d1, burn.in=bi, max.lag=1000, filename= "")
acf.Gas(d2, burn.in=bi, max.lag=200, filename= "")
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

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/SmeanGas.pdf"
#plot.Sg(d1, burn.in=bi, filename ="")
#plot.Sg(d2, burn.in=bi, filename =filename)
#plot.Sg(d3, burn.in=bi, filename ="")


#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/SmeanOil.pdf"
#plot.So(d1, burn.in=bi, filename ="")
#plot.So(d2, burn.in=bi, filename =filename)
#plot.So(d3, burn.in=bi, filename ="")


#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/SmeanBrine.pdf"
#plot.Sb(d1, burn.in=bi, filename = "")
#plot.Sb(d2, burn.in=bi, filename = filename)
#plot.Sb(d3, burn.in=bi, filename = "")
#

#filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/VmeanClay.pdf"
#plot.Vcl(d1, burn.in=bi, filename = "")
#plot.Vcl(d2, burn.in=bi, filename = "")
#plot.Vcl(d3, burn.in=bi, filename = "")



# Plot mu_prior ----
filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/muPriorGas.pdf"
plot.mu.prior.gas(inline.idx, xline.idx, filename = filename, transformation = 1)


filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/muPriorOil.pdf"
plot.mu.prior.oil(inline.idx, xline.idx, filename = filename, transformation = 1)


filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/muPriorBrine.pdf"
plot.mu.prior.brine(inline.idx, xline.idx, filename =filename)


filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/muPriorClay.pdf"
plot.mu.prior.clay(inline.idx, xline.idx, filename = filename, transformation = 1)



source("Alvheim_visualize.R")
filename <- "/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/Uncertainty"
percentile(d=d2, inline.idx, xline.idx, bi=500, p=80, filename=filename)


# Results from Mina ----

realis.g <- readRDS("results_karen/realis_sat_g_ver1.rds") 
realis.o <- readRDS("results_karen/realis_sat_o_ver1.rds")
realis.b <- realis.b(realis.g, realis.o)
realist.c <- readRDS("results_karen/realis_clay_ver1.rds") 


percentile.sat.g <- readRDS("results_karen/uncert_sat_g_ver1.rds")
percentile.sat.o <- readRDS("results_karen/uncert_sat_o_ver1.rds")
percentile.sat.b <- percentile.brine(realis.b, inline.idx, xline.idx)
percentile.sat.c <- readRDS("results_karen/uncert_vclay_ver1.rds")


source("Alvheim_visualize.R")


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


pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanGas.pdf", width=6, height=6)
plot.func(sat.g, inline.idx, xline.idx)
dev.off()

pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanOil.pdf", width=6, height=6)
plot.func(sat.o, inline.idx, xline.idx)
dev.off()

pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanBrine.pdf", width=6, height=6)
plot.func(sat.b, inline.idx, xline.idx)
dev.off()

pdf("/Users/karenauestad/Desktop/var2023/prosjektoppgave/Alvheim/plots/KalmanClay.pdf", width=6, height=6)
plot.func(sat.c, inline.idx, xline.idx)
dev.off()















