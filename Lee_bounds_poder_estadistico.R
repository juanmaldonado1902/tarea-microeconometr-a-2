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