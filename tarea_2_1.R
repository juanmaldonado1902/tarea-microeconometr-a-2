# tarea microeconometria (todo junto)
# Pregunta 1

library(tidyverse)
library(broom)
library(lmtest)
library(sandwich)
library(car)

# Cargar
df <- read_csv("/Users/juanpablo/Desktop/india_base_final.csv",
               show_col_types = FALSE)

covars <- c("gender", "age", "religion", "caste", "education", "homeBuilt")

# Convertir a numéricas 
df <- df %>%
  mutate(
    treat = as.numeric(treat),
    across(all_of(covars), ~ as.numeric(.x))
  ) %>%
  filter(!is.na(treat))

# Balance

balance_table <- map_dfr(covars, function(var) {
  
  d <- df %>% select(treat, all_of(var)) %>% drop_na()
  
  mean_all  <- mean(d[[var]])
  mean_t1   <- mean(d[[var]][d$treat == 1])
  mean_t0   <- mean(d[[var]][d$treat == 0])
  diff      <- mean_t1 - mean_t0
  
  pval <- t.test(d[[var]] ~ d$treat)$p.value
  
  tibble(
    Variable = var,
    Media_Total = mean_all,
    Media_T1 = mean_t1,
    Media_T0 = mean_t0,
    Diferencia = diff,
    p_value = pval
  )
})

balance_table <- balance_table %>%
  mutate(across(where(is.numeric), round, 4))

print(balance_table)


# regresion

formula_reg <- as.formula(
  paste("treat ~", paste(covars, collapse = " + "))
)

model <- lm(formula_reg, data = df)

# Errores robustos HC1
robust_se <- vcovHC(model, type = "HC1")
reg_results <- coeftest(model, vcov = robust_se)

print(reg_results)


# (F-test)
hypothesis <- paste0(covars, " = 0")
joint_test <- linearHypothesis(model, hypothesis, vcov = robust_se)

print(joint_test)


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

# Controles y Double Lasso

library(fixest)
library(dplyr)

# Cargar base 
india_base_final <- read.csv("/Users/juanpablo/Desktop/india_base_final.csv")

stopifnot("Pole" %in% names(india_base_final))

# Asegurar FE 
india_base_final$Pole <- as.factor(india_base_final$Pole)

# Controles 
controls <- c("gender","birthplace","age","religion","caste","education","homeBuilt")
controls <- controls[controls %in% names(india_base_final)]  # por si alguno no existe

# Fórmulas con FE = Pole
fml_adult <- as.formula(paste0(
  "adult_activity ~ treat + ", paste(controls, collapse=" + "), " | Pole"
))
fml_child <- as.formula(paste0(
  "child_activity ~ treat + ", paste(controls, collapse=" + "), " | Pole"
))

# Estimar OLS + FE 
ols_adult <- feols(fml_adult, data = india_base_final, cluster = ~Pole)
ols_child <- feols(fml_child, data = india_base_final, cluster = ~Pole)

# 7) Tabla 
etable(ols_adult, ols_child, keep = "treat", se = "cluster")

library(fixest)
library(dplyr)
library(glmnet)
library(sandwich)
library(lmtest)

# Cargar base
india_base_final <- read.csv("/Users/juanpablo/Desktop/india_base_final.csv")
india_base_final$Pole <- as.factor(india_base_final$Pole)

# Controles candidatos 
cand <- c("gender","birthplace","age","religion","caste","education","homeBuilt")
cand <- cand[cand %in% names(india_base_final)]

double_lasso_fe_PO <- function(yname, data=india_base_final, fe="Pole", dname="treat", xvars=cand){
  
  df <- data %>%
    select(all_of(c(yname, dname, fe, xvars))) %>%
    na.omit()
  
  fe_id <- df[[fe]]
  
  for(v in xvars){
    if (is.logical(df[[v]])) df[[v]] <- as.integer(df[[v]])  
  }
  if (is.logical(df[[dname]])) df[[dname]] <- as.integer(df[[dname]])
  
  demean <- function(v) v - ave(v, fe_id, FUN=function(x) mean(x, na.rm=TRUE))
  
  y <- demean(df[[yname]])
  d <- demean(df[[dname]])
  
  X <- model.matrix(as.formula(paste0("~ 0 + ", paste(xvars, collapse=" + "))), data=df)
  X <- apply(X, 2, demean)
  
  cv_y <- cv.glmnet(X, y, alpha=1)
  cv_d <- cv.glmnet(X, d, alpha=1)
  
  by <- as.vector(coef(cv_y, s="lambda.min"))[-1]
  bd <- as.vector(coef(cv_d, s="lambda.min"))[-1]
  
  sel <- sort(unique(c(which(by != 0), which(bd != 0))))
  X_sel <- if(length(sel)==0) NULL else X[, sel, drop=FALSE]
  
  if (is.null(X_sel)) {
    y_res <- y
    d_res <- d
  } else {
    y_res <- resid(lm(y ~ X_sel))
    d_res <- resid(lm(d ~ X_sel))
  }
  
  final <- lm(y_res ~ d_res)
  V <- vcovCL(final, cluster = fe_id, type="HC1")
  ct <- coeftest(final, V)
  
  list(
    final_lm = final,
    coeftest = ct,
    selected_cols = if(is.null(X_sel)) character(0) else colnames(X_sel),
    n = nrow(df)
  )
}

# Correr para ambos 
dl_adult <- double_lasso_fe_PO("adult_activity")
dl_child <- double_lasso_fe_PO("child_activity")

# Resultados 
dl_adult$coeftest
dl_child$coeftest

# columnas   LASSO
dl_adult$selected_cols
dl_child$selected_cols


# Atrición / Lee Bounds

library(tidyverse)

# Cargar
df <- readr::read_csv("/Users/juanpablo/Desktop/india_base_final.csv",
                      show_col_types = FALSE)

# Limpieza 
df <- df %>%
  mutate(
    treat = as.numeric(treat),
    total_expenditure_obs = as.numeric(total_expenditure_obs), 
    total_expenditure_isobs = ifelse(is.na(total_expenditure_obs), 0, 1)  
  ) %>%
  filter(!is.na(treat), treat %in% c(0, 1))


# Tasas de atrición 

attrition_tbl <- df %>%
  group_by(treat) %>%
  summarise(
    N = n(),
    tasa_observado = mean(total_expenditure_isobs),  
    tasa_atricion = 1 - tasa_observado              
  ) %>%
  ungroup()

print(attrition_tbl)

atr_t1 <- attrition_tbl %>% filter(treat == 1) %>% pull(tasa_atricion)
atr_t0 <- attrition_tbl %>% filter(treat == 0) %>% pull(tasa_atricion)

cat("\nGrupo con mayor atrición: ",
    ifelse(atr_t1 > atr_t0, "Tratamiento (treat=1)", "Control (treat=0)"),
    "\n", sep="")

# Lee Bounds

df_obs <- df %>% filter(total_expenditure_isobs == 1)

tau_uncorrected <- with(df_obs,
                        mean(total_expenditure_obs[treat == 1]) -
                          mean(total_expenditure_obs[treat == 0]))

cat("Tau sin corrección =", round(tau_uncorrected, 4), "\n")


# Lee trimming

p1 <- mean(df$total_expenditure_isobs[df$treat == 1])  # P(obs=1|T=1)
p0 <- mean(df$total_expenditure_isobs[df$treat == 0])  # P(obs=1|T=0)

cat("\nTasas de observación:\n")
cat("p1 = P(obs=1|T=1) =", round(p1, 4), "\n")
cat("p0 = P(obs=1|T=0) =", round(p0, 4), "\n")

trim_group <- function(vec, trim_share, drop = c("top","bottom")) {
  drop <- match.arg(drop)
  vec <- sort(vec)
  n <- length(vec)
  n_trim <- floor(trim_share * n)
  if(n_trim <= 0) return(vec)
  
  if(drop == "top") {
    vec[1:(n - n_trim)]       
  } else {
    vec[(n_trim + 1):n]        
  }
}

if (abs(p1 - p0) < 1e-12) {
  
  lee_lower <- tau_uncorrected
  lee_upper <- tau_uncorrected
  
} else if (p1 > p0) {
  
  trim_share <- (p1 - p0) / p1
  
  treated_y <- df %>%
    filter(treat == 1, total_expenditure_isobs == 1) %>%
    pull(total_expenditure_obs)
  
  control_y <- df %>%
    filter(treat == 0, total_expenditure_isobs == 1) %>%
    pull(total_expenditure_obs)
  
  treated_for_lower <- trim_group(treated_y, trim_share, drop = "top")
  lee_lower <- mean(treated_for_lower) - mean(control_y)
  
  treated_for_upper <- trim_group(treated_y, trim_share, drop = "bottom")
  lee_upper <- mean(treated_for_upper) - mean(control_y)
  
} else {
  
  trim_share <- (p0 - p1) / p0
  
  treated_y <- df %>%
    filter(treat == 1, total_expenditure_isobs == 1) %>%
    pull(total_expenditure_obs)
  
  control_y <- df %>%
    filter(treat == 0, total_expenditure_isobs == 1) %>%
    pull(total_expenditure_obs)
  
  control_for_lower <- trim_group(control_y, trim_share, drop = "bottom")
  lee_lower <- mean(treated_y) - mean(control_for_lower)
  
  control_for_upper <- trim_group(control_y, trim_share, drop = "top")
  lee_upper <- mean(treated_y) - mean(control_for_upper)
}

# Tabla final

lee_tbl <- tibble(
  Estimador = c("Sin corrección (obs)",
                "Lee - Límite inferior",
                "Lee - Límite superior"),
  Valor = c(tau_uncorrected, lee_lower, lee_upper)
) %>%
  mutate(Valor = round(Valor, 4))

print(lee_tbl)

# Poder estadístico 

library(tidyverse)

# Cargar 
df <- readr::read_csv("/Users/juanpablo/Desktop/india_base_final.csv",
                      show_col_types = FALSE) %>%
  mutate(
    treat = as.numeric(treat),
    total_expenditure = as.numeric(total_expenditure)
  ) %>%
  filter(!is.na(treat), !is.na(total_expenditure), treat %in% c(0,1))

# Estimación
y1 <- df %>% filter(treat == 1) %>% pull(total_expenditure)
y0 <- df %>% filter(treat == 0) %>% pull(total_expenditure)

tau_hat <- mean(y1) - mean(y0)

# Sigma para estandarizar
s1 <- sd(y1); s0 <- sd(y0)
n1 <- length(y1); n0 <- length(y0)

sigma_pooled <- sqrt(((n1 - 1)*s1^2 + (n0 - 1)*s0^2) / (n1 + n0 - 2))

# Efecto estandarizado 
delta <- tau_hat / sigma_pooled

cat("tau_hat =", round(tau_hat, 4), "\n")
cat("sigma_pooled =", round(sigma_pooled, 4), "\n")
cat("delta =", round(delta, 4), "\n")

# Función de poder 
power_diffmeans <- function(n, delta, alpha = 0.07, p = 0.5) {
  zcrit <- qnorm(1 - alpha/2)
  delta_n <- delta / sqrt(1/(p*n) + 1/((1-p)*n))
  1 - (pnorm(zcrit - delta_n) - pnorm(-zcrit - delta_n))
}

# Construir curva 
alpha <- 0.07
p_treat <- 0.5    
n_grid <- 20:5000  

psi <- tibble(
  n = n_grid,
  power = power_diffmeans(n, delta = delta, alpha = alpha, p = p_treat)
)

# Encontrar n mínimo 
target_power <- 0.83
n_star <- psi %>% filter(power >= target_power) %>% slice(1) %>% pull(n)

cat("n* =", n_star, "\n")

# Gráfica
plot(psi$n, psi$power, type = "l",
     xlab = "Tamaño muestral total (n)",
     ylab = "Poder estadístico ψ(n)",
     main = paste0("Curva de poder: α = ", alpha, ", delta = ", round(delta, 3), ", p = ", p_treat))

abline(h = target_power, lty = 2)
abline(v = n_star, lty = 2)
text(n_star, target_power,
     labels = paste0(" n* ≈ ", n_star, "\n power = ", target_power),
     pos = 4)

# Mostrar alrededor de n*
psi %>%
  filter(n %in% (n_star-5):(n_star+5)) %>%
  print(n = 11)

#Apéndice Double LASSO

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(sandwich)
  library(lmtest)
})

#  Helpers
mode_impute <- function(x) {
  ux <- x[!is.na(x)]
  if (length(ux) == 0) return(x)
  m <- names(sort(table(ux), decreasing = TRUE))[1]
  x[is.na(x)] <- m
  x
}
median_impute <- function(x) { m <- median(x, na.rm=TRUE); x[is.na(x)] <- m; x }
demean_vec <- function(v, fe_id) v - ave(v, fe_id, FUN=function(x) mean(x, na.rm=TRUE))
demean_mat <- function(X, fe_id) apply(X, 2, demean_vec, fe_id = fe_id)

find_by_regex <- function(nms, patterns) {
  hits <- character(0)
  for (pat in patterns) hits <- union(hits, grep(pat, nms, ignore.case=TRUE, value=TRUE))
  hits
}

find_perfect_dupes_numeric <- function(df, target_vec, vars, tol=1e-12) {
  hits <- character(0)
  y <- target_vec
  ok_y <- !is.na(y)
  for (v in vars) {
    x <- df[[v]]
    if (!is.numeric(x) && !is.integer(x)) next
    ok <- ok_y & !is.na(x)
    if (sum(ok) < 10) next
    if (sd(x[ok]) < tol) next
    corv <- suppressWarnings(cor(x[ok], y[ok]))
    if (is.finite(corv) && abs(1 - abs(corv)) < 1e-10) hits <- c(hits, v)
  }
  unique(hits)
}

prep_X_matrix <- function(dfX, fe_id,
                          impute_X=TRUE,
                          max_levels=200,
                          drop_high_cardinality=TRUE) {
  for (v in names(dfX)) {
    if (is.logical(dfX[[v]])) dfX[[v]] <- as.integer(dfX[[v]])
    if (is.character(dfX[[v]])) dfX[[v]] <- as.factor(dfX[[v]])
    if (inherits(dfX[[v]], "Date") || inherits(dfX[[v]], "POSIXct") || inherits(dfX[[v]], "POSIXt")) {
      dfX[[v]] <- as.numeric(dfX[[v]])
    }
  }
  
  if (drop_high_cardinality) {
    high_card <- names(dfX)[sapply(dfX, function(z) {
      (is.factor(z) || is.character(z)) && length(unique(z[!is.na(z)])) > max_levels
    })]
    if (length(high_card) > 0) {
      message("Dropping high-cardinality vars (> ", max_levels, " levels): ",
              paste(high_card, collapse=", "))
      dfX <- dfX[, setdiff(names(dfX), high_card), drop=FALSE]
    }
  }
  
  if (impute_X) {
    for (v in names(dfX)) {
      if (is.numeric(dfX[[v]]) || is.integer(dfX[[v]])) dfX[[v]] <- median_impute(dfX[[v]])
      if (is.factor(dfX[[v]]) || is.character(dfX[[v]])) dfX[[v]] <- as.factor(mode_impute(as.character(dfX[[v]])))
    }
  }
  
  fml <- as.formula(paste0("~ 0 + ", paste(names(dfX), collapse=" + ")))
  X <- model.matrix(fml, data=dfX)
  
  sds <- apply(X, 2, sd)
  keep <- which(is.finite(sds) & sds > 0)
  X <- X[, keep, drop=FALSE]
  
  demean_mat(X, fe_id)
}

# Double LASSO ( buscando generar después de un análisis de las variables y que sea robusto)

double_lasso_fe_robust <- function(data,
                                   yname,
                                   dname="treat",
                                   fe="Pole",
                                   exclude_vars=character(0),
                                   # patrones base a excluir
                                   id_patterns = c("(^id$|_id$|^ID$|ID$|Id$)", "uuid", "household", "^hh", "respondent"),
                                   instrument_patterns = c("^instrument$", "assign", "assignment", "encourag", "^z_", "^iv"),
                                   attrition_patterns = c("^respond$", "respond", "attrit", "nonresponse", "missing"),
                                   # (opcional) “post-treatment” generales, por si quieres mantenerlos fuera
                                   post_treat_patterns = character(0),
                                   # opciones técnicas
                                   impute_X=TRUE,
                                   drop_high_cardinality=TRUE,
                                   max_levels=200,
                                   lambda_rule=c("lambda.min","lambda.1se"),
                                   standardize=TRUE) {
  
  lambda_rule <- match.arg(lambda_rule)
  stopifnot(yname %in% names(data), dname %in% names(data), fe %in% names(data))
  
  df <- data %>% filter(!is.na(.data[[yname]]), !is.na(.data[[dname]]), !is.na(.data[[fe]]))
  fe_id <- as.factor(df[[fe]])
  
  if (is.logical(df[[dname]])) df[[dname]] <- as.integer(df[[dname]])
  if (is.factor(df[[dname]]) || is.character(df[[dname]])) df[[dname]] <- as.numeric(as.factor(df[[dname]])) - 1
  
  y <- demean_vec(df[[yname]], fe_id)
  d <- demean_vec(df[[dname]], fe_id)
  
  reserved <- unique(c(yname, dname, fe, exclude_vars))
  xvars <- setdiff(names(df), reserved)
  
  # excluir por patrones
  auto_excl <- unique(c(
    find_by_regex(xvars, id_patterns),
    find_by_regex(xvars, instrument_patterns),
    find_by_regex(xvars, attrition_patterns),
    find_by_regex(xvars, post_treat_patterns)
  ))
  xvars <- setdiff(xvars, auto_excl)
  
  # excluir duplicados perfectos de y o d (numéricos), ( tratar de eliminar variables que sean construcciones)
  dup_y <- find_perfect_dupes_numeric(df, df[[yname]], xvars)
  dup_d <- find_perfect_dupes_numeric(df, df[[dname]], xvars)
  xvars <- setdiff(xvars, unique(c(dup_y, dup_d)))
  
  # X matrix
  X <- prep_X_matrix(df[, xvars, drop=FALSE], fe_id,
                     impute_X=impute_X,
                     max_levels=max_levels,
                     drop_high_cardinality=drop_high_cardinality)
  
  # LASSO
  cv_y <- cv.glmnet(X, y, alpha=1, standardize=standardize)
  cv_d <- cv.glmnet(X, d, alpha=1, standardize=standardize)
  
  by <- as.vector(coef(cv_y, s=lambda_rule))[-1]
  bd <- as.vector(coef(cv_d, s=lambda_rule))[-1]
  
  sel <- sort(unique(c(which(by != 0), which(bd != 0))))
  selected_cols <- if (length(sel)==0) character(0) else colnames(X)[sel]
  
  # partial-out
  if (length(sel)==0) {
    y_res <- y; d_res <- d
  } else {
    X_sel <- X[, sel, drop=FALSE]
    y_res <- resid(lm(y ~ X_sel))
    d_res <- resid(lm(d ~ X_sel))
  }
  
  final <- lm(y_res ~ d_res)
  V <- vcovCL(final, cluster=fe_id, type="HC1")
  ct <- coeftest(final, V)
  
  list(
    y=yname, d=dname, fe=fe,
    n=nrow(df),
    selected_cols=selected_cols,
    excluded_auto=auto_excl,
    excluded_dup=unique(c(dup_y, dup_d)),
    coeftest=ct,
    final_lm=final
  )
}

# Wrapper

run_double_lasso_activity <- function(data, yname, dname="treat", fe="Pole",
                                      lambda_rule="lambda.min") {
  
  nms <- names(data)
  
  outcome_component_excl <- character(0)
  
  if (yname == "adult_activity") {
    outcome_component_excl <- union(outcome_component_excl,
                                    grep("^adult_home_", nms, value=TRUE))
    outcome_component_excl <- union(outcome_component_excl,
                                    grep("^adult_", nms, value=TRUE))
    outcome_component_excl <- setdiff(outcome_component_excl, yname)
  }
  
  if (yname == "child_activity") {
    outcome_component_excl <- union(outcome_component_excl,
                                    grep("^child_", nms, value=TRUE))
    outcome_component_excl <- setdiff(outcome_component_excl, yname)
  }
  
  double_lasso_fe_robust(
    data = data,
    yname = yname,
    dname = dname,
    fe = fe,
    exclude_vars = outcome_component_excl,
    lambda_rule = lambda_rule,
    post_treat_patterns = character(0)
  )
}

# USO
india_base_final <- read.csv("/Users/juanpablo/Desktop/india_base_final.csv")
india_base_final$Pole <- as.factor(india_base_final$Pole)

dl_adult <- run_double_lasso_activity(india_base_final, "adult_activity")
dl_child <- run_double_lasso_activity(india_base_final, "child_activity")

dl_adult$coeftest
dl_child$coeftest

dl_adult$selected_cols
dl_child$selected_cols

# IDs, instrumento, respond, etc.
dl_adult$excluded_auto
dl_child$excluded_auto

# DOUBLE LASSO Conservadora

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(sandwich)
  library(lmtest)
})

median_impute <- function(x) { m <- median(x, na.rm=TRUE); x[is.na(x)] <- m; x }
mode_impute <- function(x) {
  ux <- x[!is.na(x)]
  if (length(ux) == 0) return(x)
  m <- names(sort(table(ux), decreasing = TRUE))[1]
  x[is.na(x)] <- m
  x
}
demean_vec <- function(v, fe_id) v - ave(v, fe_id, FUN=function(x) mean(x, na.rm=TRUE))
demean_mat <- function(X, fe_id) apply(X, 2, demean_vec, fe_id = fe_id)

prep_X_matrix <- function(dfX, fe_id, impute_X=TRUE) {
  for (v in names(dfX)) {
    if (is.logical(dfX[[v]])) dfX[[v]] <- as.integer(dfX[[v]])
    if (is.character(dfX[[v]])) dfX[[v]] <- as.factor(dfX[[v]])
  }
  if (impute_X) {
    for (v in names(dfX)) {
      if (is.numeric(dfX[[v]]) || is.integer(dfX[[v]])) dfX[[v]] <- median_impute(dfX[[v]])
      if (is.factor(dfX[[v]]) || is.character(dfX[[v]])) dfX[[v]] <- as.factor(mode_impute(as.character(dfX[[v]])))
    }
  }
  fml <- as.formula(paste0("~ 0 + ", paste(names(dfX), collapse=" + ")))
  X <- model.matrix(fml, data=dfX)
  
  sds <- apply(X, 2, sd)
  X <- X[, which(is.finite(sds) & sds > 0), drop=FALSE]
  
  demean_mat(X, fe_id)
}

double_lasso_fe_withX <- function(data, yname, dname="treat", fe="Pole", xvars,
                                  lambda_rule=c("lambda.min","lambda.1se"),
                                  standardize=TRUE) {
  lambda_rule <- match.arg(lambda_rule)
  
  stopifnot(yname %in% names(data), dname %in% names(data), fe %in% names(data))
  
  df <- data %>%
    filter(!is.na(.data[[yname]]), !is.na(.data[[dname]]), !is.na(.data[[fe]]))
  
  fe_id <- as.factor(df[[fe]])
  
  if (is.logical(df[[dname]])) df[[dname]] <- as.integer(df[[dname]])
  if (is.factor(df[[dname]]) || is.character(df[[dname]])) df[[dname]] <- as.numeric(as.factor(df[[dname]])) - 1
  
  y <- demean_vec(df[[yname]], fe_id)
  d <- demean_vec(df[[dname]], fe_id)
  
  xvars <- setdiff(xvars, c(yname, dname, fe))
  xvars <- intersect(xvars, names(df))
  
  X <- prep_X_matrix(df[, xvars, drop=FALSE], fe_id)
  
  cv_y <- cv.glmnet(X, y, alpha=1, standardize=standardize)
  cv_d <- cv.glmnet(X, d, alpha=1, standardize=standardize)
  
  by <- as.vector(coef(cv_y, s=lambda_rule))[-1]
  bd <- as.vector(coef(cv_d, s=lambda_rule))[-1]
  
  sel <- sort(unique(c(which(by != 0), which(bd != 0))))
  selected_cols <- if (length(sel)==0) character(0) else colnames(X)[sel]
  
  if (length(sel)==0) {
    y_res <- y; d_res <- d
  } else {
    X_sel <- X[, sel, drop=FALSE]
    y_res <- resid(lm(y ~ X_sel))
    d_res <- resid(lm(d ~ X_sel))
  }
  
  final <- lm(y_res ~ d_res)
  V <- vcovCL(final, cluster=fe_id, type="HC1")
  ct <- coeftest(final, V)
  
  list(
    y=yname, d=dname, fe=fe,
    n=nrow(df),
    selected_cols=selected_cols,
    coeftest=ct
  )
}

# Cargar base

india_base_final <- read.csv("/Users/juanpablo/Desktop/india_base_final.csv")
india_base_final$Pole <- as.factor(india_base_final$Pole)

# Definir candidatas (A)

# Baseline pre-tratamiento (limpio)
X_pre <- c("gender","age","birthplace","religion","caste","education","homebuilt","village")

X_ses <- c("total_expenditure")

# Excluir
always_drop <- c("treat","forcing","Pole","respond")

# Construir X 

X_A <- unique(c(X_pre, X_ses))
X_A <- setdiff(X_A, always_drop)
X_A <- intersect(X_A, names(india_base_final))

#  Double LASSO 

dl_adult_A <- double_lasso_fe_withX(india_base_final, yname="adult_activity", xvars=X_A)
dl_child_A <- double_lasso_fe_withX(india_base_final, yname="child_activity", xvars=X_A)

dl_adult_A$coeftest
dl_child_A$coeftest

dl_adult_A$selected_cols
dl_child_A$selected_cols


# DOUBLE LASSO (Amplio)

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(sandwich)
  library(lmtest)
})

# base del código

median_impute <- function(x) { m <- median(x, na.rm=TRUE); x[is.na(x)] <- m; x }
mode_impute <- function(x) {
  ux <- x[!is.na(x)]
  if (length(ux) == 0) return(x)
  m <- names(sort(table(ux), decreasing = TRUE))[1]
  x[is.na(x)] <- m
  x
}
demean_vec <- function(v, fe_id) v - ave(v, fe_id, FUN=function(x) mean(x, na.rm=TRUE))
demean_mat <- function(X, fe_id) apply(X, 2, demean_vec, fe_id = fe_id)

prep_X_matrix <- function(dfX, fe_id, impute_X=TRUE) {
  for (v in names(dfX)) {
    if (is.logical(dfX[[v]])) dfX[[v]] <- as.integer(dfX[[v]])
    if (is.character(dfX[[v]])) dfX[[v]] <- as.factor(dfX[[v]])
  }
  if (impute_X) {
    for (v in names(dfX)) {
      if (is.numeric(dfX[[v]]) || is.integer(dfX[[v]])) dfX[[v]] <- median_impute(dfX[[v]])
      if (is.factor(dfX[[v]]) || is.character(dfX[[v]])) dfX[[v]] <- as.factor(mode_impute(as.character(dfX[[v]])))
    }
  }
  fml <- as.formula(paste0("~ 0 + ", paste(names(dfX), collapse=" + ")))
  X <- model.matrix(fml, data=dfX)
  sds <- apply(X, 2, sd)
  X <- X[, which(is.finite(sds) & sds > 0), drop=FALSE]
  demean_mat(X, fe_id)
}

double_lasso_fe_withX <- function(data, yname, dname="treat", fe="Pole", xvars,
                                  lambda_rule=c("lambda.min","lambda.1se"),
                                  standardize=TRUE) {
  lambda_rule <- match.arg(lambda_rule)
  
  stopifnot(yname %in% names(data), dname %in% names(data), fe %in% names(data))
  
  df <- data %>%
    filter(!is.na(.data[[yname]]), !is.na(.data[[dname]]), !is.na(.data[[fe]]))
  
  fe_id <- as.factor(df[[fe]])
  
  if (is.logical(df[[dname]])) df[[dname]] <- as.integer(df[[dname]])
  if (is.factor(df[[dname]]) || is.character(df[[dname]])) df[[dname]] <- as.numeric(as.factor(df[[dname]])) - 1
  
  y <- demean_vec(df[[yname]], fe_id)
  d <- demean_vec(df[[dname]], fe_id)
  
  forbidden <- c(yname, dname, fe, "adult_activity", "child_activity")
  xvars <- setdiff(xvars, forbidden)
  xvars <- intersect(xvars, names(df))
  
  X <- prep_X_matrix(df[, xvars, drop=FALSE], fe_id)
  
  cv_y <- cv.glmnet(X, y, alpha=1, standardize=standardize)
  cv_d <- cv.glmnet(X, d, alpha=1, standardize=standardize)
  
  by <- as.vector(coef(cv_y, s=lambda_rule))[-1]
  bd <- as.vector(coef(cv_d, s=lambda_rule))[-1]
  
  sel <- sort(unique(c(which(by != 0), which(bd != 0))))
  selected_cols <- if (length(sel)==0) character(0) else colnames(X)[sel]
  
  if (length(sel)==0) {
    y_res <- y; d_res <- d
  } else {
    X_sel <- X[, sel, drop=FALSE]
    y_res <- resid(lm(y ~ X_sel))
    d_res <- resid(lm(d ~ X_sel))
  }
  
  final <- lm(y_res ~ d_res)
  V <- vcovCL(final, cluster=fe_id, type="HC1")
  ct <- coeftest(final, V)
  
  list(y=yname, d=dname, fe=fe, n=nrow(df), selected_cols=selected_cols, coeftest=ct)
}

# Bases

india_base_final <- read.csv("/Users/juanpablo/Desktop/india_base_final.csv")
india_base_final$Pole <- as.factor(india_base_final$Pole)

# candidatas

X_pre <- c("gender","age","birthplace","religion","caste","education","homebuilt","distance","village")

# posibles

X_D_except_elecExp <- c(
  "elec_value","elec_hours","hours_avail","outages","low_voltage",
  "lighting_hours","adult_lighting","child_lighting",
  "kerosene_lamps","num_kerosene_lamps","kerosene_lamp_hours","kerosene_other","kerosene_expenditure",
  "appliances","appliances_use",
  "television","radio","refrigerator","mobile","fans","cooler","iron","inverter","water_pump",
  "ic_bulb","cfl_bulb","led_light","tube_light",
  "satisfaction","knowledge","aspirations","business_interest","income_increase",
  "yrsConnected"
)

X_ses <- c("total_expenditure")

always_drop <- c("treat","forcing","Pole","respond","electricity_expenditure",
                 "food_expenditure","education_expenditure")  

X_B <- unique(c(X_pre, X_ses, X_D_except_elecExp))
X_B <- setdiff(X_B, always_drop)
X_B <- intersect(X_B, names(india_base_final))

# Double LASSO

dl_adult_B <- double_lasso_fe_withX(india_base_final, yname="adult_activity", xvars=X_B)
dl_child_B <- double_lasso_fe_withX(india_base_final, yname="child_activity", xvars=X_B)

dl_adult_B$coeftest
dl_child_B$coeftest

dl_adult_B$selected_cols
dl_child_B$selected_cols