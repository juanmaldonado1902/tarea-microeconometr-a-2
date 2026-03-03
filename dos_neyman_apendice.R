# Neyman

library(tidyverse)
library(lmtest)
library(sandwich)

# Cargar 
df <- readr::read_csv("/Users/juanpablo/Desktop/india_base_final.csv",
                      show_col_types = FALSE)

# Limpieza 
df <- df %>%
  mutate(
    treat = as.numeric(treat),
    total_expenditure = as.numeric(total_expenditure)
  ) %>%
  filter(!is.na(treat), !is.na(total_expenditure)) %>%
  filter(treat %in% c(0, 1))


# Estimador de Neyman 

y1 <- df %>% filter(treat == 1) %>% pull(total_expenditure)
y0 <- df %>% filter(treat == 0) %>% pull(total_expenditure)

tau_hat_neyman <- mean(y1) - mean(y0)

cat("E[Y|T=1] =", round(mean(y1), 4), "\n")
cat("E[Y|T=0] =", round(mean(y0), 4), "\n")
cat("Tau_hat  =", round(tau_hat_neyman, 4), "\n")

# (Opcional útil) EE clásico de diferencia de medias (Neyman)
n1 <- length(y1); n0 <- length(y0)
s1 <- var(y1);    s0 <- var(y0)
se_neyman <- sqrt(s1/n1 + s0/n0)


# coeficiente de OLS 

m <- lm(total_expenditure ~ treat, data = df)

print(summary(m)$coefficients)

beta_hat <- coef(m)["treat"]



# EE robustos HC1 

rob <- coeftest(m, vcov. = vcovHC(m, type = "HC1"))

print(rob)


#  magnitud 

mean_control <- mean(y0)
effect_level <- tau_hat_neyman
effect_pct <- 100 * effect_level / mean_control

cat("\n==============================\n")
cat("INTERPRETACIÓN DE MAGNITUD\n")
cat("==============================\n")
cat("Media control (T=0):", round(mean_control, 4), "\n")
cat("Efecto en nivel (rupias/mes):", round(effect_level, 4), "\n")
cat("Efecto como % de la media control:", round(effect_pct, 2), "%\n")

#Neyman (con datos de varianza y las ns para generar la tabla)

# Cargar base
df <- read.csv("/Users/juanpablo/Desktop/india_base_final.csv")

df <- df[!is.na(df$total_expenditure) & !is.na(df$treat), ]

y <- df$total_expenditure

Tt <- df$treat

# Separar 
y1 <- y[Tt == 1]
y0 <- y[Tt == 0]

# Tamaños 
n_T <- length(y1)
n_C <- length(y0)

# Medias
mean_T <- mean(y1)
mean_C <- mean(y0)

# Varianzas muestrales
var_T <- var(y1)
var_C <- var(y0)

# Estimador de Neyman
tau_hat <- mean_T - mean_C

# Varianza 
var_tau_hat <- var_T/n_T + var_C/n_C

# Error estándar
se_tau_hat <- sqrt(var_tau_hat)

# Estadístico t
t_stat <- tau_hat / se_tau_hat

# tabla
tabla_neyman <- data.frame(
  n_tratamiento = n_T,
  n_control = n_C,
  var_tratamiento = round(var_T,4),
  var_control = round(var_C,4),
  tau_hat = round(tau_hat,4),
  var_tau_hat = round(var_tau_hat,6),
  se_tau_hat = round(se_tau_hat,4),
  t_stat = round(t_stat,4)
)

tabla_neyman