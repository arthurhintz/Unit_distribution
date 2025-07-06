
source("data_transformacao.R")

# Gráficos Análise Descritiva

grafico <- ggplot(tobacco_data, aes(x = Education, y = Resposta, fill = Gender)) +
  geom_boxplot() +
  labs(
    title = "",
    x = "Educação",
    y = "Consumo de cigarro (%)",
    fill = "Gênero"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Centraliza o texto no eixo x
    legend.title = element_text(size = 12),  # Tamanho do título da legenda
    legend.text = element_text(size = 11)  # Tamanho do texto da legenda
  ) +
  scale_fill_brewer(palette = "Paired")

ggsave("graph1.png", plot = grafico, width = 8, height = 5, dpi = 200, bg = "white")
#==========/==========/==========/==========/==========/==========/==========/==========/

ggplot(tobacco_data, aes(x = factor(YEAR), y = Resposta, fill = TopicDesc)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "",
    x = "Ano",
    y = "Consumo de cigarro (%)",
    fill = "TopicDesc"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(legend.position = c(0.8, 0.9)) +
  scale_x_discrete(drop = FALSE)


#==========/==========/==========/==========/==========/==========/==========/==========/

#corplot


colnames(tobacco_data_dummy) <- c("Ano", "Resposta", "Tamanho_Amostra", "Latitude", "Longitude", 
                                  "Tabaco_Mascavel", "Ja_fumou", "Fuma_regularmente",
                                  "Masculino","Ensino_Fund")


png("corrplot_tobacco.png", width = 8, height = 5, units = "in", res = 200)

corrplot::corrplot(cor(tobacco_data_dummy), type = "lower")

dev.off()
