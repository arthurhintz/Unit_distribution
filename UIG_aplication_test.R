setwd("~/estatistica/semestre 7/modelo taxas proporcoes")

library(tidyverse)
library(gamlss)
data(cable)
attach(cable)

# fitting the beta regression
completo_beta<-gamlss(pen5~.,
                      sigma.formula =~.,family=BE(),
                      data=cable,trace=F)
final_beta<-stepGAIC(completo_beta)

final_beta <- gamlss(formula = pen5 ~ lin + ltv + agehe, 
                     sigma.formula = ~ lin, 
                     family = BE(), data = cable, trace = FALSE)
summary(final_beta)

#==========/==========/==========/==========/==========/==========/==========/==========/
# fitting the simplex regression
completo_simplex<-gamlss(pen5~.,
                         sigma.formula =~.,family=SIMPLEX(),
                         data=cable,trace=F)
final_simplex<-stepGAIC(completo_simplex)

final_simplex <- gamlss(formula = pen5 ~ lin + ltv + agehe, 
                        sigma.formula = ~ child + ltv, 
                       family = SIMPLEX(), data = cable, trace = FALSE)

summary(final_simplex)

# fitting the UG regression

source("uiggamlss.R")

completo_UIG <- gamlss(pen5~.,
                    sigma.formula =~.,
                    family=UIG(sigma.link = "log"),
                    data=cable,trace=F)


final_UIG<-stepGAIC(completo_UIG)
final_UIG <- gamlss(formula = pen5 ~ lin + ltv + agehe, 
                    sigma.formula = ~ child + ltv + agehe, 
                    family = UIG(sigma.link = "log"), data = cable, trace = FALSE)

summary(final_UIG)

#two equivalent ways to obtain predictions
cbind(final_beta$mu.fv,
      predict(final_beta,what="mu",newdata = cable,
              type = "response"))

# some goodness-of-fit measures
IC <-data.frame(
  AIC=c(final_beta$aic,final_simplex$aic,final_UIG$aic),
  BIC=c(final_beta$sbc,final_simplex$sbc,final_UIG$sbc),
  Rsq=c(Rsq(final_beta),Rsq(final_simplex),Rsq(final_beta))
)
results<-cbind(IC,
               rbind(
                 forecast::accuracy(final_beta$mu.fv,pen5),
                 forecast::accuracy(final_simplex$mu.fv,pen5),
                 forecast::accuracy(final_UIG$mu.fv,pen5)
               )
)
rownames(results)<-c("Beta","Simplex","UIG")
results

# residual analysis for the beta
plot(final_beta)
wp(final_beta)
shapiro.test(final_beta$residuals)

hist(final_beta$residuals,freq = F)
curve(dnorm,add=T)    

# residual analysis for the simplex
# plot(final_simplex)
wp(resid=na.omit(final_simplex$residuals))
shapiro.test(final_simplex$residuals)

hist(final_simplex$residuals,freq = F)
curve(dnorm,add=T)  

# residual analysis for the UG
plot(final_UIG)
wp(final_UIG,ylim.all = T)
shapiro.test(final_UIG$residuals)

hist(final_UIG$residuals,freq = F)
curve(dnorm,add=T)    

