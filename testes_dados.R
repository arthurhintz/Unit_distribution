library(gamlss)
library(car)
library(dplyr)


setwd("~/estatistica/semestre 7/modelo taxas proporcoes")
dados_base <- read.csv("tobacco_data.csv")


table(dados_base$YEAR)


anos <- unique(dados_base$YEAR)

str(dados_base)

summary(dados_base)

hist(dados_base$Resposta)

modelo_lm <- lm(Resposta ~ ., data = dados_base)
car::vif(modelo_lm)


# Data frame para armazenar os resultados
results <- data.frame(
  ANO = character(),
  Modelo = character(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

for(k in seq_along(anos)){
  
  c <- anos[k]
  
  banco <- dados_base |> 
    filter(YEAR == c) |> 
    filter(Resposta > 0 & Resposta < 100) |> 
    mutate(Resposta = Resposta/100) |> 
    mutate(Latitude = Latitude/100) |>
    mutate(Longitude = Longitude/100) |> 
    select(3:5, 7:11,13,14, 16) |> 
    na.omit()
  
  try({
    fit1 <- gamlss(Resposta ~ ., 
                   sigma.formula = ~., 
                   family = BE(), 
                   method = RS(),
                   data = banco, trace = FALSE)
    
    results <- rbind(results, data.frame(ANO = c, Modelo = "BE", AIC = AIC(fit1)))
  }, silent = TRUE)
  
  try({
    fit2 <- gamlss(Resposta ~ .,
                   sigma.formula = ~.,
                   family = UIG(),
                   method = RS(),
                   data = banco, trace = FALSE)
    
    results <- rbind(results, data.frame(ANO = c, Modelo = "UIG", AIC = AIC(fit2)))
  }, silent = TRUE)
  
  try({
    fit3 <- gamlss(Resposta ~ ., 
                   sigma.formula = ~., 
                   family = CUW(), 
                   method = RS(),
                   data = banco, trace = FALSE)
    
    results <- rbind(results, data.frame(ANO = c, Modelo = "CUW", AIC = AIC(fit3)))
  }, silent = TRUE)
  
  try({
    fit4 <- gamlss(Resposta ~ ., 
                   family = LB(), 
                   data = banco, trace = FALSE)
    
    results <- rbind(results, data.frame(ANO = c, Modelo = "LB", AIC = AIC(fit4)))
  }, silent = TRUE)
}

results

melhores_modelos <- results %>%
  group_by(ANO) %>%
  slice_min(order_by = AIC, with_ties = FALSE) %>%
  ungroup()



write.table(results, file = "aic.txt")

completo_beta<-gamlss(TX_EVASAO ~ .,
                      sigma.formula = TX_EVASAO ~ . , 
                      family=BE(),
                      data = banco, trace=T)

completo_beta$aic


completo_UIG <- gamlss(TX_EVASAO ~ .,
                       sigma.formula = ~. ,
                       family = UIG(mu.link = "log", sigma.link = "log"),
                       method = CG(),
                       control = gamlss.control(c.crit = 0.001, n.cyc = 200,
                                              mu.step = .5, sigma.step = .5,
                                              trace= T),
                       data=banco)

summary(completo_UIG)

y <- dados$TX_EVASAO
a <- dUIG(y)
summary(a)


completo_CUW <- gamlss(TX_EVASAO ~.,
                       sigma.formula = ~.,
                       family = CUW(sigma.link = "log"),
                       data = banco,trace=F)
completo_CUW$aic


completo_LB <- gamlss(TX_EVASAO~.,
                      family=LB(mu.link = "logit"),
                      data=banco, trace=F)


completo_LB

completo_simplex<-gamlss(TX_EVASAO~.,
                         sigma.formula = ~.,family=SIMPLEX(),
                         data=banco,trace=F)




#==========/==========/==========/==========/==========/==========/==========/==========/

setwd("~/estatistica/semestre 7/modelo taxas proporcoes")
dados_base <- read.csv("tobacco_data.csv")



table(dados_base$YEAR)


dados <- dados_base |> 
  filter(YEAR == "2000") |> 
  filter(Resposta > 0 & Resposta < 100) |> 
  mutate(Resposta = Resposta/100) |> 
  mutate(Latitude = Latitude/100) |>
  mutate(Longitude = Longitude/100) |> 
  select(3:5, 7:11,13,14, 16) |> 
  na.omit()


summary(dados)


completo_beta<-gamlss(Resposta ~ .,
                      sigma.formula = ~.,
                      family=BE(),
                      data = dados, trace=T)

completo_beta$aic

summary(completo_beta)

completo_UIG <- gamlss(Resposta ~ .,
                       sigma.formula = ~.,
                       family = UIG(),
                       method = RS(),
                       control = gamlss.control(c.crit = 0.001, n.cyc = 200,
                                                mu.step = .1, sigma.step = .1,
                                                trace= T),
                       data=dados)

summary(completo_UIG)

completo_UIG$aic


completo_CUW <- gamlss(Resposta ~.,
                       sigma.formula = ~.,
                       family = CUW(sigma.link = "log"),
                       data = dados,trace=F)
completo_CUW$aic


completo_LB <- gamlss(Resposta~.,
                      family=LB(mu.link = "logit"),
                      data=dados, trace=F)


completo_LB$aic

completo_simplex<-gamlss(Resposta~.,
                         sigma.formula = ~.,family=SIMPLEX(),
                         data=dados,trace=F)


summary(completo_simplex)
completo_simplex$aic


#==========/==========/==========/==========/==========/==========/==========/==========/

completo_UIG <- gamlss(Resposta ~ .,
                       sigma.formula = ~.,
                       family = UIG(),
                       method = RS(),
                       control = gamlss.control(c.crit = 0.001, n.cyc = 200,
                                                mu.step = .1, sigma.step = .1,
                                                trace= T),
                       data=dados)


final_UIG <- stepGAIC(completo_UIG)

summary(final_UIG)

final_UIG$call

final <- gamlss(formula = Resposta ~ Longitude + TopicDesc_Cigarette.Use..Youth. + 
         TopicDesc_Smokeless.Tobacco.Use..Youth. + MeasureDesc_Percent.of.Current.Smokers.Who.Want.to.Quit + 
         Education_High.School, 
         sigma.formula = ~ TopicDesc_Cigarette.Use..Youth. + Gender_Male + Education_High.School, family = UIG(), 
       data = dados, method = RS(),
       trace = FALSE)

summary(final)

plot(final)
wp(final,ylim.all = T)
shapiro.test(final$residuals)

hist(final$residuals,freq = F)
curve(dnorm,add=T)    
