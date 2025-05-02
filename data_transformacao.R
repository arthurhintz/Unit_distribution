# Pacotes
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(fastDummies)

# Importar base de tabaco
tobacco_data_raw <- read_csv("Tobacco_Survey_Data.csv")

tobacco_data <- na.omit(tobacco_data_raw)



# Seleção das variáveis
tobacco_data <- tobacco_data %>% 
  select(
    -LocationAbbr,
    -Data_Value_Std_Err,
    -Low_Confidence_Limit,
    -High_Confidence_Limit,
    -Sample_Size,
    -StratificationID1,
    -StratificationID4,
    -DisplayOrder,
    -SubMeasureID,
    -MeasureDesc
  )
# Transformação das cordenadas
tobacco_data <- tobacco_data %>%
  mutate(GeoLocation = str_remove_all(GeoLocation, "[()]")) %>%
  separate_wider_delim(
    cols = "GeoLocation",
    delim = ",",
    names = c("Latitude", "Longitude")
  ) %>%
  mutate(
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

colnames(tobacco_data)

# Transformação das variáveis dummies
tobacco_data <- dummy_cols(tobacco_data, select_columns = c("TopicDesc","Response","Gender","Education"),
                           remove_selected_columns = TRUE, remove_first_dummy = T)



write.csv(tobacco_data, "tobacco_data.csv", row.names = FALSE)

cor(tobacco_data)
#==========/==========/==========/==========/==========/==========/==========/==========/

modelo_lm <- lm(Resposta ~ ., data = tobacco_data)
car::vif(modelo_lm)


banco <- tobacco_data |> 
  filter(YEAR == 2000) |> 
  filter(Resposta > 0 & Resposta < 100) |> 
  mutate(Resposta = Resposta/100) |> 
  mutate(Latitude = Latitude/100) |>
  mutate(Longitude = Longitude/100)



be <- gamlss(Resposta ~ ., 
             sigma.formula = ~., 
             family = BE(), 
             method = RS(),
             data = banco, trace = FALSE)


be$aic



uig <- gamlss(Resposta ~ .,
              sigma.formula = ~.,
              family = UIG(),
              method = RS(),
              data = banco, trace = FALSE)


uig$aic

cuw <- gamlss(Resposta ~ ., 
              sigma.formula = ~., 
              family = CUW(), 
              method = RS(),
              data = banco, trace = FALSE)


cuw$aic  

lb <- gamlss(Resposta ~ ., 
             family = LB(), 
             data = banco, trace = FALSE)

lb$aic  


simplex <- gamlss(Resposta ~ ., 
                  family = SIMPLEX(), 
                  data = banco, trace = FALSE)


simplex$aic  


#==========/==========/==========/==========/==========/==========/==========/==========/




final_UIG <- stepGAIC(uig)


summary(final_UIG)

final_UIG$call

final <- gamlss(formula = Resposta ~ Latitude + `TopicDesc_Smokeless Tobacco Use (Youth)` + 
                  Response_Ever + Response_Frequent + Gender_Male + Gender_Overall + 
                  `Education_Middle School` + YEAR_2000 + YEAR_2001 + YEAR_2002 + 
                  YEAR_2003 + YEAR_2004 + YEAR_2005 + YEAR_2006 + YEAR_2007 + 
                  YEAR_2008 + YEAR_2009 + YEAR_2010 + YEAR_2011 + YEAR_2012 + 
                  YEAR_2013 + YEAR_2014 + YEAR_2015 + YEAR_2016 + YEAR_2017, 
                sigma.formula = Resposta ~ Latitude + `TopicDesc_Smokeless Tobacco Use (Youth)` + 
                  Response_Ever + Response_Frequent + Gender_Overall + 
                  `Education_Middle School` + YEAR_2000 + YEAR_2002 + 
                  YEAR_2009 + YEAR_2010 + YEAR_2011 + YEAR_2012 + YEAR_2014 + YEAR_2015, 
                family = UIG(), data = banco, method = RS(), 
                trace = FALSE)

summary(final)

plot(final)
wp(final,ylim.all = T)
shapiro.test(final$residuals)

hist(final$residuals,freq = F)
curve(dnorm,add=T)    

