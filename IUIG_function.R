library(gamlss)
# Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023
#------------------------------------------------------------------------------------------ 
# initial values

source("UIG_functions.R")

#------------------------------------------------------------------------------------------
mu = 0.58; sigma = 1; nu = 0.1; x=.22
#------------------------------------------------------------------------------------------
# Kumaraswamy distribution - basic functions
#------------------------------------------------------------------------------------------
# Zero-inflated Gaussiana Inversa distribution - basic functions
#------------------------------------------------------------------------------------------
# density function for x in [0,1)

dIUIG <- function(x, mu = 0.5, sigma = 2, nu = 0.1, log = FALSE){
  
  if (any(x < 0 | x > 1)) stop("x must be in (0, 1)")
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(sigma <= 0)) stop("lambda must be positive")
  
  # if (any(x == 0)) {
  #   if (lambda / (2 * mu^2) < 1) return(Inf)
  #   else return(0)
  # }
  
  if (any(x == 1)) return(0)
  

  log.d = (1/2) * (log(sigma) - log(2*pi)) - (log(x *(-log(x))^(3/2))) + 
    ((sigma * (log(x) + mu) ^ 2) / (2 * mu ^ 2 * log(x)))
  
  
  log.lik <- ifelse(x == 0, log(nu), log(1 - nu) + log.d)
  
  if (log == FALSE) 
    fy <- exp(log.lik)
  else fy <- log.lik
    fy <- ifelse(x < 0 | x >= 1, 0, fy)
  
  return(fy)
}
  
nu = 0.1
dIUIG(0, nu = nu) # = nu

dIUIG(x, mu, sigma)
(1-nu)*dUIG(x, mu, sigma)
dIUIG(1) # 0 because it is not inflated in one
#------------------------------------------------------------------------------------------ 
# cumulative distribution function


pIUIG<-function (q,  mu = 0.5, sigma = 2, nu = 0.1, 
                 lower.tail = TRUE, log.p = FALSE) 
{
  if (any(q < 0 | q >= 1)) stop("x must be in (0, 1)")
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(sigma <= 0)) stop("sigma must be positive")
  
  if(q > 0 & q < 1) {
  term1 <- sqrt(sigma / -log(q)) * ((-log(q) / mu) - 1)
  term2 <- sqrt(sigma / -log(q)) * ((-log(q) / mu) + 1)
  
  cdf <-  nu + (1 - nu) * 
          (1 - (pnorm(term1) + exp(2 * sigma / mu) * pnorm(-term2)))
  
  } else 0

  cdf <- ifelse((q == 0), nu, cdf)
  cdf <- ifelse((q >= 1), 1, cdf)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf <- ifelse(q < 0, 0, cdf)
  cdf <- ifelse(q >= 1, 1, cdf)
  cdf
}
# checking pZIKUM
pIUIG(0)
x=0.4
pIUIG(x, mu = mu, sigma = sigma)
nu+(1-nu)*pUIG(x, mu, lambda = sigma)
#------------------------------------------------------------------------------------------
# quantile function
qIUIG<-function (p, mu = 0.7, sigma = 2.1, nu = 0.1, lower.tail = TRUE, 
                  log.p = FALSE) {
  
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be positive"))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive"))
  if (any(nu <= 0) | any(nu >= 1)) 
    stop(paste("nu must be beetwen 0 and 1 "))
  
  if (log.p == TRUE) p <- exp(p)
  else p <- p
  
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  
    u <- (p - nu)/(1 - nu)
  
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  suppressWarnings(

    if(nu >= p) 0
    else {     
    
      ext <- uniroot(function(x) pUIG(x, mu, lambda) - u, 
                   interval = c(1e-10, 1 - 1e-10),
                   tol = 1e-300, extendInt = "yes")$root
    
      q <- abs(ext)
    
    }
  )
  return(q)
}


# checking qIUIG
x
u = pIUIG(x, mu = mu, sigma)

qIUIG(u, mu, sigma)

qUIG((u - nu)/(1 - nu),mu = mu, lambda = sigma)

#------------------------------------------------------------------------------------------
# inversion method for random generation

rIUIG<-function (n, mu = 0.7, sigma = 2.1, nu = 0.1){
  
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0) | any(nu >= 1)) 
    stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))

  if(n != length(mu)) mu = rep(mu, n)
  if(n != length(lambda)) lambda = rep(lambda, n)
  if(n != length(nu)) nu = rep(nu, n)
  
  n <- ceiling(n)
  p <- runif(n)
  

  r <- rep(NA, n)
  for (i in 1:n) {
    r[i] <- qIUIG(p[i], mu[i], lambda[i], nu[i])
  }
  r
  

  
  return(r)

}


rIUIG(1)
