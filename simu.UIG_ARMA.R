setwd("~/estatistica/semestre 7/modelo taxas proporcoes")
source("UIG_functions.R")


simu.UIG_ARMA <- function(n, sigma, alpha, phi=NULL, theta=NULL, freq=12, link="logit"){
  
  ifelse((length(phi) == 1), ar <- 1, stop("phi deve receber um valor"))
  
  ifelse((length(theta) == 1), ma <- 1, stop("theta deve receber um valor"))
  
    
  link <- make.link(link)
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rUIG(1,mu[i], sigma)
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]   
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
}

yhat <- simu.UIG_ARMA(100, 2, 0.4, -0.5, alpha = 0.5)


sd(yhat)
