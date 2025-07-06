library(gamlss)
source("OIUIG.R")

# OIUIG in gamlss.family----------------------------------------------------

# derivatives of log.UIG
UIG <- expression(
  (1/2) * (log(sigma) - log(2*pi)) - (log(y *(-log(y))^(3/2))) + 
    ((sigma * (log(y) + mu) ^ 2) / (2 * mu ^ 2 * log(y)))
)
m1 <- D(UIG, "mu")
s1 <- D(UIG, "sigma")
ms2 <- D(m1, "sigma")

# Gamlss family function
OIUIG <- function(mu.link = "logit", sigma.link = "log", nu.link = "logit") {
  mstats <- checklink(
    "mu.link", "OIUIG", substitute(mu.link),
    c("logit", "probit", "cloglog", "log", "own")
  )
  dstats <- checklink(
    "sigma.link", "OIUIG", substitute(sigma.link),
    c("inverse", "log", "identity")
  )
  vstats <- checklink(
    "nu.link", "OIUIG", substitute(nu.link),
    c("logit", "probit", "cloglog", "log", "own", "identity")
  )
  structure(
    list(
      family = c("OIUIG", "One Inflated Unit Gaussian Inversa"),
      parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
      nopar = 3, type = "Mixed",
      mu.link = as.character(substitute(mu.link)),
      sigma.link = as.character(substitute(sigma.link)),
      nu.link = as.character(substitute(nu.link)),
      mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, nu.linkfun = vstats$linkfun,
      mu.linkinv = mstats$linkinv, sigma.linkinv = dstats$linkinv, nu.linkinv = vstats$linkinv,
      mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta, nu.dr = vstats$mu.eta,
      dldm = function(y, mu, sigma) {
        dldm <- ifelse((y == 1), 0, eval(m1))
        dldm
      }, d2ldm2 = function(y, mu, sigma) {
        dldm <- eval(m1)
        d2ldm2 <- ifelse((y == 1), 0, -dldm * dldm)
        d2ldm2
      }, dldd = function(y, mu, sigma) {
        dldd <- ifelse((y == 1), 0, eval(s1))
        dldd
      }, d2ldd2 = function(y, mu, sigma) {
        dldd <- ifelse(is.nan(eval(s1)), 0, eval(s1))
        d2ldd2 <- ifelse((y == 1), 0, -(dldd * dldd))
        d2ldd2
      }, dldv = function(y, nu) {
        dldv <- ifelse(y == 1, 1 / nu, -1 / (1 - nu))
        dldv
      }, d2ldv2 = function(nu) {
        d2ldv2 <- -1 / (nu * (1 - nu))
        d2ldv2
      }, d2ldmdd = function(y, mu, sigma) {
        dldm <- eval(m1)
        dldd <- ifelse(is.nan(eval(s1)), 0, eval(s1))
        d2ldmdd <- ifelse((y == 1), 0, -(dldm * dldd))
        d2ldmdd
      }, d2ldmdv = function(y) {
        d2ldmdv <- rep(0, length(y))
        d2ldmdv
      }, d2ldddv = function(y) {
        d2ldddv <- rep(0, length(y))
        d2ldddv
      }, G.dev.incr = function(y, mu, sigma, nu, ...) {
        -2 * dOIUIG(y, mu, sigma, nu, log = TRUE)
      }, rqres = expression({
        uval <- ifelse(y == 1, nu * runif(length(y), 0, 1),
                       (1 - nu) * pOIUIG(y, mu, sigma, nu)
        )
        rqres <- qnorm(uval)
      }),
      mu.initial = expression(mu <- rep(median(y), length(y))),
      sigma.initial = expression(sigma <- rep(1, length(y))),
      nu.initial = expression(nu <- rep(0.3, length(y))),
      mu.valid = function(mu) all(mu > 0 & mu < 1),
      sigma.valid = function(sigma) all(sigma > 0),
      nu.valid = function(nu) all(nu > 0 & nu < 1),
      y.valid = function(y) all(y > 0 & y <= 1)
    ),
    class = c("gamlss.family", "family")
  )
}


