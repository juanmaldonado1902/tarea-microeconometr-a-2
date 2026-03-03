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