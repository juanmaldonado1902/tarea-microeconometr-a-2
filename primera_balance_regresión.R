
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
