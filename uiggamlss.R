#setwd("~/estatistica/semestre 7/modelo taxas proporcoes")
#source("UIG_functions.R")
#library(gamlss)

#log-likelihood of the unit gaussian inverse
UIG <- expression(
 
  log(sigma/(2*pi)) * (1 / 2) - log(y *(-log(y))^(3/2)) + 
    (sigma * (log(y) + mu) ^ 2) / (2 * mu ^ 2 * log(y))
)

# teste
# sigma = 2
# mu = 0.5
# y = 0.4
# log(dUIG(y, mu, sigma))

#==========/==========/==========/==========/==========/==========/==========/==========/

d1m <- D(UIG,"mu")
d1l <- D(UIG,"sigma")
d2ml <- D(d1m,"sigma")


UIG <- function (mu.link = "logit", sigma.link = "identity"){
  
  mstats <- checklink("mu.link", "UIG", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "UIG", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("UIG", "Unit-Inverse-Gaussiana"),
                 parameters = list(mu = TRUE, sigma = TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 dldm = function(y, mu, sigma) {
                   dldm <- eval(d1m)
                   dldm
                 },
                 d2ldm2 = function(y,mu, sigma) {
                   dldm <- eval(d1m)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
                   d2ldm2
                 },
                 dldd = function(y, mu, sigma) {
                   dldd <- eval(d1l)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   dldd <- eval(d1l)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   dldm <- eval(d1m)
                   dldd <- eval(d1l)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd
                 },
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * 
                   log(dUIG(x = y, mu = mu, lambda = sigma)),
                 rqres = expression(
                   rqres(pfun = "pUIG", type = "Continuous", y = y, mu = mu, lambda = sigma)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),
                 sigma.initial = expression(sigma<- rep(4, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}