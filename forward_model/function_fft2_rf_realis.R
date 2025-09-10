
fft2_rf_realis <- function(nt,nx,sig,ksi,phi){
  
  #set.seed(26)
  #nt <- 181
  #nx <- 251
  nv <- c(nx,nt)
  
  n <- nx*nt
  
  # Spacing 
  
  delta1 <- 1
  # delta_2 <- 0.1
  delta2 <- 1
  
  # mean increasing with depth
  incdepth <- 0.015*(1:nx)
  
  # dimensions seem incorrect
  #mu <- 0 
  #mu <- -1*matrix(1,nt,nx)+matrix(1,nt,1)%*%matrix(incdepth,1,nx)
  #mu <- t_cov
  
  # variance parameters
  #sig <- 1
  # sig for porosity
  #sig <- 0.008
  # sigma otherwise
  #sig <- 0.5
  #ksi <- 0.01
  #phi <- 20
  
  
  # prior variance
  sig2 <- sig*sig
  
  # hyper prior variance
  ksi2 <- ksi*ksi
  # base of Toeplitz circular matrix 
  n1 <- nv[2]
  n2 <- nv[1]
  cm <- matrix(0,n1,n2)
  cmf <- matrix(0,n1,n2)
  
  for( i1 in (0:(n1-1))){
    tau1 <- delta1*(i1)
    if((i1)>((n1-1)/2)){
      tau1 <- delta1*(n1-i1)
    }
    for(j1 in (0:(n2-1))){
      tau2 <- delta2*(j1)
      if((j1) > ((n2-1)/2)){
        tau2 <- delta2*(n2-j1)
      }
      
      if(phi!=0){
        cm[i1+1,j1+1] <- sig2*exp(-0.5*(tau1^2+tau2^2)/phi^2)
      }else{
       print("Error: corr=0")
        break
      }
      
      
    }
  }
  
  #cmf <- Re(fft(cm))
  cmf <- Re(fft(cm))
  #cm[1,1] <- cm[1,1]+ksi2
  
  uv <- rnorm(n)
  u <- matrix(0,n1,n2)
  u <- matrix(uv, n1,n2)
  
  Ltoep <- sqrt(matrix(as.complex(cmf),n1,n2))
  # Here maybe I need to divide mftt(u) with its length
  zf <- Re(fft(Ltoep*fft(u,inverse=TRUE)/sqrt(n))/sqrt(n))
  
  return(zf)
}

