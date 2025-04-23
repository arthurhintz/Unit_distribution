# Simulação Sem Regressor
source("uiggamlss.R")

set.seed(1248)
n<-1000
R <- 100
# Case 1: without regressors

logit_link <- make.link("logit")
mu_true <- 0.5
sigma_true <- 2
mu_result <- sigma_result <- c()

for (i in 1:R) {
  y <- rUIG(n,mu_true, sigma_true)
  fit <- gamlss(y~1, family="UIG", trace = F)
  mu_result[i] <- logit_link$linkinv(fit$mu.coefficients)
  sigma_result[i] <- fit$sigma.coefficients
}
result1<- matrix(c(mu_true, mean(mu_result),
                   sigma_true, mean(sigma_result)),2,2)
colnames(result1)<-c("mu","sigma")
rownames(result1)<-c("true value","mean")
print(round(result1,4))


# Com Regressor

X <- runif(n)
logit_link <- make.link("logit")
log_link <- make.link("identity")
b1 <- .6
b2 <- 2.5
mu_true <- logit_link$linkinv(b1+b2*X)

g1 <- .5
g2 <- 1.5
sigma_true <- log_link$linkinv(g1+g2*X)
mu_result <- sigma_result <- matrix(NA,R,2)

for (i in 1:R) {
  y <- rUIG(n, mu_true, sigma_true)
  fit1 <- gamlss(y~X,sigma.formula = ~ X, family=UIG(), trace = F)
  mu_result[i,] <- fit1$mu.coefficients
  sigma_result[i,] <- fit1$sigma.coefficients
}

true_values<-c(b1,b2, g1,g2)
mean_values<-c(apply(mu_result,2,mean),
               apply(sigma_result,2,mean))
b_values<-(true_values-mean_values)/true_values*100
eqm_values<-c(apply(mu_result,2,var),
              apply(sigma_result,2,var))+(true_values-mean_values)^2
result1<- cbind(true_values,
                mean_values,
                b_values,
                eqm_values
)
colnames(result1)<-c("true value","mean","bias","eqm")
rownames(result1)<-c("b1","b2","g1","g2")
print(round(result1,2))
