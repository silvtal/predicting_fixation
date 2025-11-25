# made with cursor

library(dilgrowth)
library(ggplot2)
library(tidyverse)

# Función auxiliar para crear matriz de interacciones
create_interaction_matrix <- function(n, value) {
  mat <- matrix(0, n, n)
  mat[3,1] <- value  # Interacción de P1 sobre P3
  return(mat)
}

# Configuración inicial
set.seed(123)
initial_abundance <- c(100, 100, 100, 100)  # P1, P2, P3, P4 con abundancia inicial 100
names(initial_abundance) <- c("P1", "P2", "P3", "P4")

# Definir capacidades de carga para cada grupo funcional
# Grupo 1 (P1, P2) tiene mayor capacidad de carga que Grupo 2 (P3, P4)
carrying_capacities <- c(0.6, 0.6, 0.4, 0.4)  # 60% para Grupo 1, 40% para Grupo 2
names(carrying_capacities) <- c("G1", "G1", "G2", "G2")

dilfactor <- 0.8
numdil <- 10

# 1. Sin interacciones
sim_no_interaction <- simulate_timeseries(
  dilution = dilfactor,
  force_continue = TRUE,
  counts_data = initial_abundance,
  carrying_capacities = carrying_capacities,
  no_of_dil = numdil,
  keep_all_timesteps = TRUE
)

# 2. Con interacción negativa
sim_negative <- simulate_timeseries(
  dilution = dilfactor,
  counts_data = initial_abundance,
  carrying_capacities = carrying_capacities,
  interactions = create_interaction_matrix(4, -1),
  no_of_dil = numdil,
  keep_all_timesteps = TRUE
)

# 3. Con interacción positiva
sim_positive <- simulate_timeseries(
  dilution = dilfactor,
  counts_data = initial_abundance,
  carrying_capacities = carrying_capacities,
  interactions = create_interaction_matrix(4, 0.5),
  no_of_dil = numdil,
  keep_all_timesteps = TRUE
)

# Preparar datos para ggplot
prepare_data <- function(sim_data, title) {
  colnames(sim_data) <- names(initial_abundance)
  df <- data.frame(
    time = 0:numdil,
    P1 = sim_data$P1,
    P2 = sim_data$P2,
    P3 = sim_data$P3,
    P4 = sim_data$P4
  )
  
  # Reemplazar NA con el último valor no-NA manualmente
  for(col in c("P1", "P2", "P3", "P4")) {
    last_val <- NA
    for(i in 1:nrow(df)) {
      if(is.na(df[i,col])) {
        df[i,col] <- last_val
      } else {
        last_val <- df[i,col]
      }
    }
  }
  
  df %>%
    pivot_longer(
      cols = c(P1, P2, P3, P4),
      names_to = "Population",
      values_to = "Abundance"
    ) %>%
    mutate(
      Scenario = title,
      Group = ifelse(Population %in% c("P1", "P2"), "Group 1 (CC: 60%)", "Group 2 (CC: 40%)")
    )
}

# Combinar todos los datos
all_data <- rbind(
  prepare_data(sim_no_interaction, "𝗔                No interactions                 "),
  prepare_data(sim_negative, "𝗕   Negative interaction (P1→P3 = -1)  "),
  prepare_data(sim_positive, "𝗖   Positive interaction (P1→P3 = 0.5)  ")
)  %>% 
  mutate(Scenario = factor(Scenario, 
                          levels = c("𝗔                No interactions                 ", 
                                   "𝗕   Negative interaction (P1→P3 = -1)  ",
                                   "𝗖   Positive interaction (P1→P3 = 0.5)  ")))

# Crear las figuras
p <- ggplot(all_data, aes(x = time, y = Abundance, color = Population, linetype = Group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Scenario) +
  theme_minimal() +
  scale_x_continuous(breaks = 0:10) +  # Forzar números enteros en el eje x
  labs(
    x = "Number of cycles",
    y = "Abundance",
    # title = "Simulación de crecimiento con diferentes tipos de interacciones",
    # subtitle = "Grupo Funcional 1 (P1, P2): Capacidad de carga 60%\nGrupo Funcional 2 (P3, P4): Capacidad de carga 40%"
  ) +
  scale_color_manual(values = c(
    "P1" = "#E41A1C",  # Rojo
    "P2" = "#FB8072",  # Rojo claro
    "P3" = "#377EB8",  # Azul
    "P4" = "#80B1D3"   # Azul claro
  )) +
  theme(
    strip.text = element_text(size = 13),
    # plot.title = element_text(hjust = 0.5, size = 16),
    # plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Guardar la figura
ggsave("SupplFig1_interaction_simulations.png", p, width = 14, height = 5)