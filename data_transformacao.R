# Pacotes
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(fastDummies)
library(ggplot2)


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
    #-Sample_Size,
    -StratificationID1,
    -StratificationID4,
    -DisplayOrder,
    -SubMeasureID,
    -MeasureDesc
  ) |> 
  filter(Gender != "Overall")


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

set_data <- tobacco_data

tobacco_data <- set_data

# Transformação das variáveis dummies
tobacco_data_dummy <- dummy_cols(tobacco_data, select_columns = c("TopicDesc","Response","Gender","Education"),
                           remove_selected_columns = TRUE, remove_first_dummy = T)



write.csv(tobacco_data_dummy, "tobacco_data.csv", row.names = FALSE)
