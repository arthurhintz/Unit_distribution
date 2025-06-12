# Distribuição Gaussiana Unitária Inflacionada em 1

source("UIG_functions.R")

# função densidade inflacionada em 1
dOIUIG <- function(x, mu = 0.5, sigma = 2, nu = 0.1, log = FALSE){
  
  if (any(x < 0 | x > 1)) stop("x must be in (0, 1)")
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(sigma <= 0)) stop("lambda must be positive")
  
  log.d = (1/2) * (log(sigma) - log(2*pi)) - (log(x *(-log(x))^(3/2))) + 
    ((sigma * (log(x) + mu) ^ 2) / (2 * mu ^ 2 * log(x)))
  
  log.lik <- ifelse(x == 1, log(nu), log(1 - nu) + log.d)
  
  ifelse((log == FALSE), fx <- exp(log.lik), fx <- log.lik)
  
  fx <- ifelse(x <= 0 | x > 1, 0, fx)
  
  return(fx)
}

# nu = 0.2
# dOIUIG(1, nu = nu)
# 
# mu = 0.58; sigma = 1; x=.22
# dOIUIG(x, mu, sigma, nu)
# (1-nu)*dUIG(x, mu, sigma)


#==========/==========/==========/==========/==========/==========/==========/==========/


# Função Acumulada

pOIUIG<-function (q,  mu = 0.5, sigma = 2, nu = 0.1, 
                  lower.tail = TRUE, log.p = FALSE) 
{
  if (any(q < 0 | q > 1)) stop("x must be in (0, 1)")
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(sigma <= 0)) stop("sigma must be positive")
  
  if(any(q > 0 & q < 1)) {
    term1 <- sqrt(sigma / -log(q)) * ((-log(q) / mu) - 1)
    term2 <- sqrt(sigma / -log(q)) * ((-log(q) / mu) + 1)
    
    cdf <-  (1 - nu) * 
      (1 - (pnorm(term1) + exp(2 * sigma / mu) * pnorm(-term2)))
    
  } else {cdf <- 1}
  
  ifelse(lower.tail == TRUE, cdf, cdf <- 1 - cdf)
  ifelse(log.p == FALSE, cdf, cdf <- log(cdf))
  
  return(cdf)
}
# checking pOIUIG
# nu
# pOIUIG(0.99,mu,sigma,nu)
# 
# x=0.4
# pOIUIG(x, mu = mu, sigma = sigma, nu)
# (1-nu)*pUIG(x, mu, lambda = sigma)

#==========/==========/==========/==========/==========/==========/==========/==========/
# quantile function

qOIUIG<-function (p, mu = 0.7, sigma = 2.1, nu = 0.1, lower.tail = TRUE, 
                 log.p = FALSE) {
  
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be positive"))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive"))
  if (any(nu <= 0) | any(nu >= 1)) 
    stop(paste("nu must be beetwen 0 and 1 "))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  
  ifelse((log.p == FALSE), p, p <- exp(p)) 
  ifelse((lower.tail == TRUE), p, p <- 1 - p)
  
  u <- (p) / (1 - nu) 
  
  if(p > (1 - nu)) {q <- 1
  
  }else{
    ext <- uniroot(function(x) pUIG(x, mu, sigma) - u, 
                   interval = c(1e-10, 1 - 1e-10),
                   tol = 1e-300, extendInt = "yes")$root
    
    q <- abs(ext)
  }
            
  return(q)
}


# checking qOIUIG
# x = 0.2
# nu
# u = pOIUIG(x, mu = mu, sigma,nu)
# 
# qOIUIG(u, mu, sigma, nu)
# 
# qUIG((u)/(1 - nu),mu = mu, lambda = sigma)


#==========/==========/==========/==========/==========/==========/==========/==========/


rOIUIG<-function (n, mu = 0.7, sigma = 2.1, nu = 0.1){
  
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0) | any(nu >= 1)) 
    stop(paste("nu must be beetwen 0 and 1 ", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  
  if(n != length(mu)) mu = rep(mu, n)
  if(n != length(sigma)) sigma = rep(sigma, n)
  if(n != length(nu)) nu = rep(nu, n)
  
  n <- ceiling(n)
  p <- runif(n)
  
  
  r <- rep(NA, n)
  for (i in 1:n) {
    r[i] <- qOIUIG(p[i], mu[i], sigma[i], nu[i])
  }
  r
  
  return(r)
  
}

# sampling <- rOIUIG(1000, mu, sigma,nu)
# hist(sampling, breaks = 10, xlim = c(0,1))