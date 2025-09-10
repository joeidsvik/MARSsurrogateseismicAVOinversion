


# Plot Alvheim on small area

source("Alvheim_data.R")
source("Alvheim_plot_aux.R")

library(LaplacesDemon)

# 30x30 area of interest ----

inline.idx = seq(156, 215, 2)
xline.idx = seq(71, 130, 2)

N <- length(inline.idx) * length(xline.idx)

m <- 1000000
save.every <- 10 


# Results ----
# s = 0.027
b <- m/(2*save.every) 

## MARS ----
load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/MCMC_runs/small_area/RW2_MARS_027.Rda")
d2$accepted/m 
plot.Sg(d2, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_233_gas.pdf")  
plot.So(d2, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_233_oil.pdf")   
plot.Sb(d2, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_233_brine.pdf") 
plot.Vcl(d2, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_233_clay.pdf") 
percentile(d=d2, inline.idx, xline.idx, bi=b, p=80, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_233_uncertainty")
rm(d2)

## MARS w Omega tilde ----
load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/MCMC_runs/small_area/RW2_MARS_tilde_027.Rda")
d2.tilde$accepted/m
plot.Sg(d2.tilde, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_tilde_234_gas.pdf")  
plot.So(d2.tilde, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_tilde_234_oil.pdf")   
plot.Sb(d2.tilde, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_tilde_234_brine.pdf") 
plot.Vcl(d2.tilde, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_tilde_234_clay.pdf") 
percentile(d=d2.tilde, inline.idx, xline.idx, bi=b, p=80, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_tilde_234_uncertainty")
rm(d2.tilde)

## SFM ----
load("/Users/karenauestad/Dokumenter/master/Karen_codes_3/MCMC_runs/small_area/RW2_sfm_027.Rda")
d2.sfm$accepted/m
plot.Sg(d2.sfm, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_sfm_235_gas.pdf")  
plot.So(d2.sfm, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_sfm_235_oil.pdf")   
plot.Sb(d2.sfm, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_sfm_235_brine.pdf") 
plot.Vcl(d2.sfm, burn.in = b, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_sfm_235_clay.pdf") 
percentile(d=d2.sfm, inline.idx, xline.idx, bi=b, p=80, filename="/Users/karenauestad/Dokumenter/master/plots_master/Alvheim/small_area/q2_sfm_235_uncertainty")
rm(d2.sfm)


