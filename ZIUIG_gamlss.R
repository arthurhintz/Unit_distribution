library(gamlss)
#------------------------------------------------------------------------------------------
# ZIUIG in gamlss.family
#------------------------------------------------------------------------------------------


source("IUIG_function.R")
# derivatives of log.kum
UIG <- expression(
  (1/2) * (log(sigma) - log(2*pi)) - (log(y *(-log(y))^(3/2))) + 
    ((sigma * (log(y) + mu) ^ 2) / (2 * mu ^ 2 * log(y)))
)

m1<-D(UIG,"mu")
s1<-D(UIG,"sigma")
ms2<-D(m1,"sigma")
#

ZIUIG<-function (mu.link = "log", sigma.link = "log", nu.link = "logit") 
{
  mstats <- checklink("mu.link", "ZIUIG", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "log", "own"))
  dstats <- checklink("sigma.link", "ZIUIG", substitute(sigma.link), 
                      c("inverse", "log", "identity"))
  vstats <- checklink("nu.link", "ZIUIG", substitute(nu.link), 
                      c("logit", "probit", "cloglog", "log", "own"))
  structure(list(family = c("ZIUIG", "Zero Inflated Gaussian Inverse"), 
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), 
                 nopar = 3, type = "Mixed", 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)), 
                 mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, nu.linkfun = vstats$linkfun, 
                 mu.linkinv = mstats$linkinv, sigma.linkinv = dstats$linkinv, nu.linkinv = vstats$linkinv, 
                 mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta, nu.dr = vstats$mu.eta, 
                 dldm = function(y, mu, sigma) {
                   tau=.5   
                   dldm <- ifelse((y == 0), 0, eval(m1))
                   dldm
                 }, d2ldm2 = function(y, mu, sigma) { 
                   tau=.5   
                   dldm <- eval(m1)  
                   d2ldm2 <- ifelse((y == 0), 0, -dldm * dldm)      
                   d2ldm2
                 }, dldd = function(y, mu, sigma) {
                   tau=.5
                   dldd <- ifelse((y == 0), 0, eval(s1))
                   dldd
                 }, d2ldd2 = function(y, mu, sigma) { 
                   tau=.5                  
                   dldd <- ifelse(is.nan(eval(s1)), 0, eval(s1))  
                   d2ldd2 <- ifelse((y == 0), 0, -(dldd * dldd))
                   d2ldd2
                 }, dldv = function(y, nu) {
                   dldv <- ifelse(y == 0, 1/nu, -1/(1 - nu))
                   dldv
                 }, d2ldv2 = function(nu) {
                   d2ldv2 <- -1/(nu * (1 - nu))
                   d2ldv2
                 }, d2ldmdd = function(y, mu, sigma) {
                   tau=.5   
                   dldm <- eval(m1)
                   dldd <- ifelse(is.nan(eval(s1)), 0, eval(s1))  
                   d2ldmdd <- ifelse((y == 0), 0, -(dldm * dldd))
                   d2ldmdd
                 }, d2ldmdv = function(y) {
                   d2ldmdv <- rep(0, length(y))
                   d2ldmdv
                 }, d2ldddv = function(y) {
                   d2ldddv <- rep(0, length(y))
                   d2ldddv
                 }, G.dev.incr = function(y, mu, sigma, nu, ...) {
                   -2 * dIUIG(y, mu, sigma, nu, log = TRUE)
                 }, rqres = expression({
                   uval <- ifelse(y == 0, nu * runif(length(y), 0, 1), 
                                  (1 - nu) * pIUIG(y, mu, sigma, nu))
                   rqres <- qnorm(uval)
                 }), 
                 mu.initial = expression(mu <- rep(median(y),length(y))),
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial = expression(nu <- rep(0.3, length(y))), 
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu) all(nu > 0 & nu < 1), 
                 y.valid = function(y) all(y >= 0 & y < 1)), 
            class = c("gamlss.family", "family")
  )
}

# Checking the results
library(gamlss)
set.seed(10)
n<-50
# Case 1: without regressors
mu_true<- 0.9
sigma_true<-2.6
nu_true<-.4
mu_result<-sigma_result<-nu_result<-c()
for (i in 1:100) {
  y<-rIUIG(n,mu_true,sigma_true,nu_true)
  fit1<-gamlss(y~1,family="ZIKUM",trace=F)
  mu_result[i]<-fit1$mu.fv[1]
  sigma_result[i]<-fit1$sigma.fv[1]
  nu_result[i]<-fit1$nu.fv[1]
}

result1<- matrix(c(mu_true, mean(mu_result),
                   sigma_true, mean(sigma_result),
                   nu_true, mean(nu_result)),2,3)
colnames(result1)<-c("mu","sigma","nu")
rownames(result1)<-c("true value","mean")
print(round(result1,2))