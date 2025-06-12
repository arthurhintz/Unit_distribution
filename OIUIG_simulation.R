# Loading packages and functions -----------------------------------------------

library(gamlss)
source("OIUIG_gamlss.R")

# Simulation -------------------------------------------------------------------
set.seed(8)

# Case 1: without regressors
n <- 100
mu_true <- 0.8
sigma_true <- 1.2
nu_true <- 0.2
mu_result <- sigma_result <- nu_result <- c()
for (i in 1:100) {
  y <- rOIUIG(n, mu_true, sigma_true, nu_true)
  fit1 <- gamlss(y ~ 1, family = "OIUIG", trace = F)
  mu_result[i] <- fit1$mu.fv[1]
  sigma_result[i] <- fit1$sigma.fv[1]
  nu_result[i] <- fit1$nu.fv[1]
}

result1 <- matrix(c(
  mu_true, mean(mu_result),
  sigma_true, mean(sigma_result),
  nu_true, mean(nu_result)
), 2, 3)
colnames(result1) <- c("mu", "sigma", "nu")
rownames(result1) <- c("true value", "mean")
print(round(result1, 2))


# Case 2: with regressors
n <- 1000
X <- runif(n)
log_link <- make.link("log")
logit_link <- make.link("logit")
b1 <- 1.6
b2 <- -0.7
mu_true <- logit_link$linkinv(b1 + b2 * X)
g1 <- 0.6
g2 <- -1
sigma_true <- log_link$linkinv(g1 + g2 * X)
nu_true <- 0.15
R <- 300
mu_result <- matrix(NA, R, 2)
sigma_result <- matrix(NA, R, 2)
nu_result <- matrix(NA, R, 1)

for (i in 1:R) {
  y <- rOIUIG(n, mu_true, sigma_true, nu_true)
  fit_reg <- gamlss(
    y ~ X, sigma.formula = ~X, nu.formula = ~1,
    family = OIUIG(mu.link = "logit", sigma.link = "log", nu.link = "identity"),
    trace = F, method = RS()
  )
  mu_result[i, ] <- fit_reg$mu.coefficients
  sigma_result[i, ] <- fit_reg$sigma.coefficients
  nu_result[i, ] <- fit_reg$nu.coefficients
}

true_values <- c(b1, b2, g1, g2, nu_true)
mean_values <- c(
  apply(mu_result, 2, mean),
  apply(sigma_result, 2, mean),
  apply(nu_result, 2, mean)
)
b_values <- (true_values - mean_values) / true_values * 100
eqm_values <- c(
  apply(mu_result, 2, var),
  apply(sigma_result, 2, var),
  apply(nu_result, 2, var)
) + (true_values - mean_values)^2
result1 <- cbind(
  true_values,
  mean_values,
  b_values,
  eqm_values
)
colnames(result1) <- c("true value", "mean", "bias", "eqm")
rownames(result1) <- c("b1", "b2", "g1", "g2", "nu")
print(round(result1, 2))
