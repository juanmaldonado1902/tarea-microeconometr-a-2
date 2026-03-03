
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