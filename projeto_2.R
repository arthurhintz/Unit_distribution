library(gamlss)
library(car)
library(dplyr)
library(fastDummies)
library(tidyverse)
library(yardstick)


# Arthur UIG
source("https://raw.githubusercontent.com/arthurhintz/Unit_distribution/refs/heads/main/UIG_functions.R")
source("https://raw.githubusercontent.com/arthurhintz/Unit_distribution/main/uiggamlss.R")


# Joca
source("https://raw.githubusercontent.com/JoaquimReis/LB-gamlss/refs/heads/main/LB_gamlss.R")
source("https://raw.githubusercontent.com/JoaquimReis/LB-gamlss/refs/heads/main/Functions_LBgamlss.R")

# Renata Kumaraswamy
source("https://raw.githubusercontent.com/renata-rojasg/UnitDistForGAMLSS/refs/heads/main/KW_reg.R")

# Renata RUBXII
source("https://raw.githubusercontent.com/renata-rojasg/UnitDistForGAMLSS/refs/heads/main/RUBXII_reg.R")

# UG
source("https://raw.githubusercontent.com/renata-rojasg/UnitDistForGAMLSS/refs/heads/main/UG_reg.R")

# ULindley
source("https://raw.githubusercontent.com/renata-rojasg/UnitDistForGAMLSS/refs/heads/main/ULindley_reg.R")

# UW
source("https://raw.githubusercontent.com/renata-rojasg/UnitDistForGAMLSS/refs/heads/main/UW_reg.R")

#==========/==========/==========/==========/==========/==========/==========/==========/

dados_base <- read.csv("tobacco_data.csv")

plot(dados_base$Ano, dados_base$Resposta)

str(dados_base)

dados_base$classes_anos <- cut(
  dados_base$Ano,
  breaks = c(1998, 2002, 2007, 2012, 2018),
  labels = c("99-02", "03-07", "08-12", "13-17"),
  right = TRUE
)


dados_base_dummies <- dummy_cols(
  dados_base,
  select_columns = "classes_anos",
  remove_selected_columns = TRUE,  # remove a variável categórica original
  remove_first_dummy = TRUE        # evita multicolinearidade (opcional)
)


dados <- dados_base_dummies |> 
  filter(Resposta > 0 & Resposta < 100) |> 
  mutate(Resposta = Resposta/100) |> 
  mutate(Latitude = Latitude/100) |>
  mutate(Longitude = Longitude/100) |> 
  mutate(Tamanho_Amostra = log(Tamanho_Amostra)) |> 
  select(2:13) |> 
  na.omit()


summary(dados)
  
modelo_lm <- lm(Resposta ~ ., data = dados)
vif(modelo_lm)
  
# AJustando Todos Os Modelos Completos

# beta
completo_beta<-gamlss(Resposta ~ .,
                      sigma.formula = ~.,
                      family=BE(),
                      control = gamlss.control(c.crit = 0.001, n.cyc = 200,
                                               mu.step = .1, sigma.step = .1,
                                               trace= T),
                      data = dados)

completo_beta$aic
plot(completo_beta)
Rsq(completo_beta)

png("plot_beta.png", width = 8, height = 5, units = "in",res = 200)
plot(completo_beta)
dev.off()

summary(completo_beta)

# UIG
completo_UIG <- gamlss(Resposta ~ .,
                       sigma.formula = ~.,
                       family = UIG(mu.link = "log", sigma.link = "log"),
                       method = RS(),
                       control = gamlss.control(c.crit = 0.001, n.cyc = 200,
                                                mu.step = .1, sigma.step = .1,
                                                trace= T),
                       data=dados)

summary(completo_UIG)
plot(completo_UIG)
completo_UIG$aic

summary(fitted(completo_UIG))     # valores ajustados
summary(residuals(completo_UIG))  # resíduos

#LB
completo_LB <- gamlss(Resposta~.,
                      family=LB(mu.link = "logit"),
                      data=dados, trace=F)


completo_LB$aic
Rsq(completo_LB)

png("plot_logbilal.png", width = 8, height = 5, units = "in",res = 200)
plot(completo_LB)
dev.off()

# Simplex
completo_simplex<-gamlss(Resposta~.,sigma.formula = ~.,
                         family=SIMPLEX(),
                         control = gamlss.control(c.crit = 0.001, n.cyc = 200,
                                                  mu.step = .1, sigma.step = .1,
                                                  trace= T),
                         data = dados)



summary(completo_simplex)
Rsq(completo_simplex)
completo_simplex$aic
plot(completo_simplex)


# Kumaraswamy

kw <- gamlss(Resposta ~ ., sigma.formula = ~ .,
             family = KW(sigma.link = "log"),
             control = gamlss.control(c.crit = 0.001, n.cyc = 200,
                                       mu.step = .1, sigma.step = .1,
                                      trace= T),
             data=dados)



kw$aic  
Rsq(kw)

png("plot_kw.png", width = 8, height = 5, units = "in",res = 200)
plot(kw)
dev.off()


# BurXII
completo_bxii <-gamlss(Resposta ~ .,
                      sigma.formula = ~.,
                      family=RUBXII(), mu.start = 0.3, 
                      control = gamlss.control(n.cyc = 200, mu.step = .1,
                                               sigma.step = .1, trace= T),
                      data = dados)


completo_bxii$aic

png("plot_burr.png", width = 8, height = 5, units = "in",res = 200)
plot(completo_bxii)
dev.off()





plot(completo_bxii)
summary(completo_bxii$y)
Rsq(completo_bxii)



# UG

completo_ug <-gamlss(Resposta ~ .,
                       sigma.formula = ~.,
                       family=UG(), 
                       control = gamlss.control(n.cyc = 200, mu.step = .2,
                                                sigma.step = .2, trace= T),
                       data = dados)

completo_ug$aic

png("plot_ug.png", width = 8, height = 5, units = "in",res = 200)
plot(completo_ug)
dev.off()

Rsq(completo_ug)

rmse <- sqrt(mean(completo_ug$residuals^2))
rmse

# Lindley

completo_lindley <-gamlss(Resposta ~ .,
                     family=ULindley(), mu.start = 0.2,
                     control = gamlss.control(n.cyc = 200, mu.step = .2,
                                              trace= T),
                     data = dados)

completo_lindley$aic

png("plot_lindley.png", width = 8, height = 5, units = "in",res = 200)
plot(completo_lindley)
dev.off()

Rsq(completo_lindley)
completo_lindley$y


# UW

completo_uw <-gamlss(Resposta ~ ., sigma.formula = ~.,
                          family=UW(), 
                          control = gamlss.control(n.cyc = 200, mu.step = .2,
                                                   trace= T),
                          data = dados)

completo_uw$aic

png("plot_uw.png", width = 8, height = 5, units = "in",res = 200)
plot(completo_uw)
dev.off()

plot(completo_uw)
Rsq(completo_uw)

#==========/==========/==========/==========/==========/==========/==========/==========/


final_beta <- stepGAIC(completo_beta)

summary(final_beta)

final_beta$call

final <- gamlss(formula = Resposta ~ Tamanho_Amostra + Latitude + Tabaco_Mascavel + 
                   Ja_fumou + Fuma_regularmente + Masculino + Ensino_Fund + 
                   `classes_anos_03-07` + `classes_anos_08-12` + `classes_anos_13-17`, 
                 sigma.formula = Resposta ~ Latitude + Tabaco_Mascavel + 
                   Ja_fumou + Fuma_regularmente + Ensino_Fund + 
                   `classes_anos_03-07` + `classes_anos_08-12` + `classes_anos_13-17`, 
                 family = BE(), data = dados,
                 control = gamlss.control(c.crit = 0.001,
                                           n.cyc = 200, mu.step = 0.1, 
                                          sigma.step = 0.1, trace = T), 
                 trace = FALSE)

summary(final)


png("final_beta.png", width = 8, height = 5, units = "in",res = 200)
plot(final)
dev.off()


png("wp_fig.png", width = 8, height = 5, units = "in",res = 200)
wp(final,ylim.all = T)
dev.off()




hist(final$residuals,freq = F)
curve(dnorm,add=T)    


h <- hatvalues(final)
plot(h,pch="+") 

Rsq(final)

#==========/==========/==========/==========/==========/==========/==========/==========/
# Gráfico final dos coeficientes estimados

estim_mu <- final$mu.coefficients
estim_sigma <- final$sigma.coefficients


# Dados dos efeitos significativos
df_mu <- data.frame(
  Variable = names(estim_mu)[-1],
  Effect = estim_mu[-1],
  Type = "Mu (Logit)"
)

df_sigma <- data.frame(
  Variable = names(estim_sigma)[-1],
  Effect = estim_sigma[-1],
  Type = "Sigma (Log)"
)

# Combina os dois dataframes
df_combined <- bind_rows(df_mu, df_sigma)

# Gráfico
gg_estim <- ggplot(df_combined, aes(x = reorder(Variable, Effect), y = Effect, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  labs(title = "",
       x = "Variável",
       y = "Coeficiente Estimado",
       fill = "Parâmetro") +
  theme_minimal()

ggsave("plot_estim.png", plot = gg_estim, width = 6, height = 4, dpi = 200, bg = "white")

#==========/==========/==========/==========/==========/==========/==========/==========/
# Gerando a tabela para o Latex

library(xtable)

# Extraindo coeficientes para mu e sigma
mu_coefs <- summary(final)

sigma_coefs <- summary(final)

# Tabela LaTeX para o parâmetro mu
xtable(mu_coefs,
       digits = c(2, 4, 4, 3, 4),
       caption = "Coeficientes para o parâmetro \\(\\mu\\)",
       label = "tab:mu")

# Tabela LaTeX para o parâmetro sigma
xtable(sigma_coefs,
       digits = c(2, 4, 4, 3, 4),
       caption = "Coeficientes para o parâmetro \\(\\sigma\\)",
       label = "tab:sigma")

#==========/==========/==========/==========/==========/==========/==========/==========/
# Gráfico da dos grupos dos anos com a proporção

gg_anos <- ggplot(dados_base, aes(x = Ano, y = Resposta, color = classes_anos)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Ano",
    y = "Resposta",
    color = "Grupo de anos"
  )

ggsave("plot_anos.png", plot = gg_anos, width = 7, height = 4, dpi = 200, bg = "white")

#==========/==========/==========/==========/==========/==========/==========/==========/
# Gráfico da Freq anos

frequencia <- dados_base %>%
  count(Ano, classes_anos)

# Plotar gráfico de barras
ggplot(frequencia, aes(x = as.factor(Ano), y = n, fill = classes_anos)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Frequência de observações por ano",
    x = "Ano",
    y = "Frequência",
    fill = "Grupo de anos"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
