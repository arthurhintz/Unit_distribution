# Função Densidade de Probabilidade
dUIG <- function(x, mu = 0.5, lambda = 2){
  
  if (any(x < 0 | x > 1)) stop("x must be in (0, 1)")
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(lambda <= 0)) stop("lambda must be positive")
  
  if (any(x == 0)) {
    if (lambda / (2 * mu^2) < 1) return(Inf)
    else return(0)
  }
  
  if (any(x == 1)) return(0)
  
  d <-  sqrt(lambda / (2 * pi)) * (1 / (x *(-log(x))^(3/2))) * 
        exp((lambda * (log(x) + mu) ^ 2) / (2 * mu ^ 2 * log(x)))
  
  return(d)

}

# testes:
# dUIG(2)
# dUIG(0.2)
# integrate(dUIG, 0, 1)


# Função Acumulada
pUIG <- function(q, mu = 0.5, lambda = 2) {
  
  if (any(q <= 0 | q >= 1)) stop("x must be in (0, 1)")
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(lambda <= 0)) stop("lambda must be positive")
  
  term1 <- sqrt(lambda / -log(q)) * ((-log(q) / mu) - 1)
  term2 <- sqrt(lambda / -log(q)) * ((-log(q) / mu) + 1)
  
  cdf <-  1 - (pnorm(term1) + exp(2 * lambda / mu) * pnorm(-term2))
  
  return(cdf)
}

#testes
# x = .4
# integrate(dUIG,0,x, mu = 0.5, lambda = 3)
# pUIG(x, mu = 0.5, lambda = 3)


# Função Quantílica
qUIG <- function(u, mu = 0.5, lambda = 2){
  
  if (any(u <= 0 | u >= 1)) stop("u must be in (0, 1)")
  if (any(mu <= 0)) stop("mu must be positive")
  if (any(lambda <= 0)) stop("lambda must be positive")
  
    ext <- uniroot(function(x) pUIG(x, mu, lambda) - u, 
                   interval = c(1e-10, 1 - 1e-10),
                   tol = 1e-300, extendInt = "yes")$root
    return(abs(ext))
  }

# testes
# u = pUIG(0.5)  
# qUIG(u)

# inversion method for random generation
rUIG <- function(n, mu, lambda) {
  if(n != length(mu)) mu = rep(mu, n)
  if(n != length(lambda)) lambda = rep(lambda, n)
  
  u <- runif(n)
  r <- rep(NA, n)
  for (i in 1:n) {
    r[i] <- qUIG(u[i], mu[i], lambda[i])
  }
  r
}

# Verossimilhança
log_vero <- function(pars, x) {
  mu <- pars[1]
  lambda <- pars[2]
  n <- length(x)
  
  lv <- (n/2)*log(lambda / (2 * pi)) -
    sum(log(x * (-log(x))^(3/2))) +
    (lambda / (2 * mu^2)) * sum((log(x) + mu)^2 / log(x))
  
  return(-lv)  
}


# # Teste
# lambda <- 2
# mu <- 0.5
# x <- rUIG(1000, mu, lambda)
# sum(log(dUIG(x, mu, lambda)))
# - log_vero(c(mu, lambda), x)


#==========/==========/==========/==========/==========/==========/==========/==========/
# Estimação Máxima Verossimilhança
estim <- function(x) {
  
  result <- optim(par = c(1,1),
                  fn = log_vero,
                  x = x,
                  method = "SANN")
  
  return(result$par)
}

#estim(x)

#==========/==========/==========/==========/==========/==========/==========/==========/
# Simulação

# set.seed(1248)
# nrep <- 500
# n <- 100
# mu_true <- 0.8
# lambda_true <- 2.5
# 
# mu_hat <- numeric(nrep)
# lambda_hat <- numeric(nrep)
# 
# for(i in 1:nrep){
#   amostra <- rUIG(n, mu_true, lambda_true)
# 
#     est <- estim(amostra)
# 
#     mu_hat[i] <- est[1]
#     lambda_hat[i] <- est[2]
# }
# 
# media_mu <- mean(mu_hat)
# media_lambda <- mean(lambda_hat)
# 
# # viés
# bies_mu <- media_mu - mu_true
# bies_lambda <- media_lambda - lambda_true
# 
# # EQM
# eqm_mu <- mean((mu_hat - mu_true)^2)
# eqm_lambda <- mean((lambda_hat - lambda_true)^2)
# 
# # Resultados
# resultados <- data.frame(
#   Parametro = c("mu", "lambda"),
#   Verdadeiro = c(mu_true, lambda_true),
#   Estimado = round(c(media_mu, media_lambda),4),
#   Vies = round(c(bies_mu, bies_lambda), 4),
#   EQM = round(c(eqm_mu, eqm_lambda), 4)
# )
# 
# resultados
