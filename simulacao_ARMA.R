ECONV<-function(estimator)
{	
  if(is.vector(estimator)==TRUE) {
    re <- length(estimator)
    m <- rep(0, times = re)
    for (j in 1:re) {
      m[j] <- mean(estimator[1:j])
    }
    par(mfrow = c(1, 1))
    return(plot(m))
  }
  else{
    re <- length(estimator[, 1])
    ne <- length(estimator[1, ])
    m <- matrix(nrow = re, ncol = ne)
    for (j in 1:re) {
      for (k in 1:ne) {
        m[j, k] <- mean(estimator[1:j, k])
      }
    }
  }
  return(m)
}


#----------------------


set.seed(8)
R <- 500
alpha = .2; phi = .3; theta = 0.2; sigma = 2.5

vp <- c(alpha, phi, theta, sigma)
mf <- matrix(NA, ncol = length(vp), nrow = R)

nvalores = c(800)

results <- array(NA, dim = c(4, length(vp), length(nvalores)))
dimnames(results) <- list(c("parâmetros", "par_estimados", "viés_relativo", "EQM"), 
                          c("alpha", "phi", "theta", "sigma"), as.character(nvalores))

for (k in nvalores) {

  for (i in 1:R) {
    y<-simu.UIG_ARMA(k, sigma, alpha, phi, theta)
    fit <- UIG_ARMA.fit(y)
    mf[i, 1:length(vp)] <- fit$coef
  
  }
  
  par_est <- apply(mf, 2, mean)
  vviesm <- (par_est - vp)
  viest <- (par_est - vp) / vp
  vvarm <- apply(mf, 2, var)
  veqmm <- vviesm^2 + vvarm
  
  results[, , as.character(k)] <- rbind(vp, round(par_est, digits = 4), round(viest, digits = 5), round(veqmm, digits = 3))
}

results


ma<-ECONV(mf)
#----------------------
plot.ts(ma[,1],ylab="Convergência",xlab=expression(b1))
abline(h=alpha,col="red")

plot.ts(ma[,2],ylab="Convergência",xlab=expression(b2))
abline(h=phi,col="red")

plot.ts(ma[,3],ylab="Convergência",xlab=expression(b2))
abline(h=theta,col="red")

plot.ts(ma[,4],ylab="Convergência",xlab=expression(b2))
abline(h=sigma,col="red")

