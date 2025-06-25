# Modelo ARMA(1,1) da Distribuição Gauissiana Unitária

UIG_ARMA.fit <- function(y, link = "logit") {
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  if (!is.ts(y)) stop("Data must be a time-series object")
  
  n <- length(y)
  m <- 1  
  
  stats <- make.link(link)
  linkfun <- stats$linkfun
  linkinv <- stats$linkinv
  
  ynew <- linkfun(y)
  error <- rep(0, n)
  eta <- rep(0, n)
  
  # Ajuste inicial via Mínimos Quadrados
  Ynew <- ynew[(m + 1):n]
  x <- matrix(ynew[m:(n - 1)], n - m, 1) 
  x <- cbind(rep(1,(n-m)), x)
  ajuste <- lm.fit(x, Ynew)
  a_phi <- ajuste$coefficients
  reg <- c(a_phi, 0, log(1))  #alpha, phi, theta, log_sigma
  
  loglik <- function(z) {
    alpha <- z[1]
    phi <- z[2]
    theta <- z[3]
    sigma <- exp(z[4])
    
    for (i in (m + 1):n) {
      eta[i] <- alpha + as.numeric(phi%*%ynew[i-1]) + as.numeric(theta%*%error[i-1])
      error[i] <- ynew[i] - eta[i]
    }
    
    mu <- linkinv(eta[(m + 1):n])
    yt <- y[(m + 1):n]
    
    ll <- (1/2) * (log(sigma) - log(2 * pi)) -
      log(yt * (-log(yt))^(3/2)) +
      ((sigma * (log(yt) + mu)^2) / (2 * mu^2 * log(yt)))
    
    sum(ll)
  }
  
  opt <- optim(reg, loglik, method = "BFGS", hessian = TRUE,
               control = list(fnscale = -1, maxit = 50, reltol = 1e-10))
  
  z <- list()
  
  if (opt$conv != 0){warning("FUNCTION DID NOT CONVERGE!")}
    
  z$conv <- opt$convergence
  
  par_opt <- opt$par
  
  coef <- c(par_opt[1:3], sigma = exp(par_opt[4]))
  names(coef) <- c("alpha", "phi", "theta", "sigma")
  
  z$coef <- coef
  J_inv <- try(solve(opt$hessian), silent = TRUE)
  
  z$stderror <- sqrt(diag(J_inv))
  z$zstat <- coef / z$stderror
  z$pvalues <- 2 * (1 - pnorm(abs(z$zstat)))
  
  
  z$loglik <- opt$value
  k <- length(coef)
  z$aic <- -2 * z$loglik + 2 * k
  z$bic <- -2 * z$loglik + log(n) * k
  z$hq <- -2 * z$loglik + log(log(n)) * k
  
  z$model <- cbind(
    Estimate = round(z$coef, 4),
    Std.Error = round(z$stderror, 4),
    z.value = round(z$zstat, 4),
    `Pr(>|z|)` = round(z$pvalues, 4)
  )
  
  return(z)
}

# R <- 1000
# coeff<-matrix(NA,R,4)
# for(i in 1:R){
#   y <- simu.UIG_ARMA(500, sigma = 2.5, alpha = .2, phi = .3, theta = 0.2)
#   fit <- UIG_ARMA.fit(y)
#   coeff[i,]<-fit$model[,1]
# }
# 
# apply(coeff,2,mean)

