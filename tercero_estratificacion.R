# ESTIMADOR ESTRATIFICADO

library(tidyverse)
library(lmtest)
library(sandwich)

# Cargar
df <- read_csv("/Users/juanpablo/Desktop/india_base_final.csv",
               show_col_types = FALSE)
# Limpieza 
df <- df %>%
  mutate(
    treat = as.numeric(treat),
    total_expenditure = as.numeric(total_expenditure),
    age = as.numeric(age)
  ) %>%
  filter(!is.na(treat),
         !is.na(total_expenditure),
         !is.na(age))

# Definir estratos

df <- df %>%
  mutate(
    estrato = case_when(
      age <= 35 ~ "Joven",
      age >= 36 & age <= 50 ~ "Mediana",
      age >= 51 ~ "Mayor"
    )
  )

# Diferencias dentro de cada estrato


estrato_results <- df %>%
  group_by(estrato) %>%
  summarise(
    N = n(),
    media_T1 = mean(total_expenditure[treat == 1], na.rm = TRUE),
    media_T0 = mean(total_expenditure[treat == 0], na.rm = TRUE),
    diff = media_T1 - media_T0
  ) %>%
  ungroup()

# Proporciones 
estrato_results <- estrato_results %>%
  mutate(
    peso = N / sum(N)
  )

print(estrato_results)


# Estimador agregado estratificado


tau_estratificado <- sum(estrato_results$diff * estrato_results$peso)

cat("Tau_hat estratificado =", round(tau_estratificado, 4), "\n")


# Regresión dummies de estrato

model_strata <- lm(total_expenditure ~ treat + factor(estrato), data = df)

robust_strata <- coeftest(model_strata,
                          vcov = vcovHC(model_strata, type = "HC1"))

print(robust_strata)


# Neyman estratificado 2

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Cargar 
df <- read.csv("/Users/juanpablo/Desktop/india_base.csv")

# Variables
Yvar <- "total_expenditure"
Tvar <- "treat"
Avar <- "age"

# Limpiar NA
df <- df %>%
  filter(!is.na(.data[[Yvar]]),
         !is.na(.data[[Tvar]]),
         !is.na(.data[[Avar]]))

# Crear estratos
df <- df %>%
  mutate(
    estrato = case_when(
      .data[[Avar]] <= 35 ~ "Joven",
      .data[[Avar]] >= 36 & .data[[Avar]] <= 50 ~ "Mediana",
      .data[[Avar]] >= 51 ~ "Mayor",
      TRUE ~ NA_character_
    ),
    estrato = factor(estrato, levels = c("Joven", "Mediana", "Mayor"))
  ) %>%
  filter(!is.na(estrato))

# Resumen por estrato y grupo 
by_stratum <- df %>%
  group_by(estrato, treat = .data[[Tvar]]) %>%
  summarise(
    n = n(),
    mean_y = mean(.data[[Yvar]]),
    var_y  = var(.data[[Yvar]]),
    .groups = "drop"
  )

wide <- by_stratum %>%
  mutate(grupo = ifelse(treat == 1, "T", "C")) %>%
  select(-treat) %>%
  pivot_wider(
    names_from = grupo,
    values_from = c(n, mean_y, var_y),
    names_sep = "_"
  )

# Pesos por estrato
weights <- df %>%
  group_by(estrato) %>%
  summarise(n_h = n(), .groups = "drop") %>%
  mutate(w_h = n_h / sum(n_h))

wide <- wide %>%
  left_join(weights, by = "estrato") %>%
  mutate(
    tau_h = mean_y_T - mean_y_C,
    var_tau_h = (var_y_T / n_T) + (var_y_C / n_C),
    contrib_var = (w_h^2) * var_tau_h
  )

# Estimador Neyman estratificado 
tau_strat <- sum(wide$w_h * wide$tau_h, na.rm = TRUE)
var_strat <- sum(wide$contrib_var, na.rm = TRUE)
se_strat  <- sqrt(var_strat)
t_strat   <- tau_strat / se_strat

# Tabla 
tabla_estratos <- wide %>%
  transmute(
    estrato,
    n_tratamiento = n_T,
    n_control = n_C,
    var_tratamiento = var_y_T,
    var_control = var_y_C,
    tau_h = tau_h,
    peso_wh = w_h
  )

print(tabla_estratos)

library(dplyr)

# Cargar base
df <- read.csv("/Users/juanpablo/Desktop/india_base.csv")

# Limpiar NA
df <- df %>%
  filter(!is.na(total_expenditure),
         !is.na(treat),
         !is.na(age))

# Crear estratos
df <- df %>%
  mutate(
    estrato = case_when(
      age <= 35 ~ "Joven",
      age >= 36 & age <= 50 ~ "Mediana",
      age >= 51 ~ "Mayor"
    ),
    estrato = factor(estrato, levels = c("Joven","Mediana","Mayor"))
  )

# Calcular tamaño por estrato
n_total <- nrow(df)

pesos_estrato <- df %>%
  group_by(estrato) %>%
  summarise(n_h = n(), .groups="drop") %>%
  mutate(weight = n_total / n_h)

# Asignar peso individual
df <- df %>%
  left_join(pesos_estrato, by="estrato")


# Regresión ponderada (WLS)

modelo_wls <- lm(total_expenditure ~ treat + estrato,
                 data = df,
                 weights = weight)

summary(modelo_wls)

modelo_interacciones <- lm(total_expenditure ~ 0 + estrato + estrato:treat,
                           data = df)

summary(modelo_interacciones)

coef_h <- coef(modelo_interacciones)[grep("estrato.*:treat", names(coef(modelo_interacciones)))]

# pesos (generar la t estratificada pondereda)
pesos <- df %>%
  group_by(estrato) %>%
  summarise(n_h = n()) %>%
  mutate(w_h = n_h / sum(n_h))

tau_strat_reg <- sum(pesos$w_h * coef_h)
tau_strat_reg