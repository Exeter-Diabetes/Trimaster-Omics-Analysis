
# Initial setup ----

# load libraries
library(glmnet)
library(Metrics)
library(tidyverse)

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data <- set_up_dataset(differential = TRUE) %>%
  mutate(sex = as.numeric(sex)) # needed for BART (can't be a factor)

# proteomics variables
variables_proteomics <- data %>% select(contains("proteom")) %>% colnames()
# clinical features
clinical_features <- c("prebmi", "sex", "agetx", "preegfr", "t2dmduration")


data <- data %>%
  mutate_at(c(variables_proteomics, "prebmi", "agetx", "preegfr", "t2dmduration", "prehba1c"), function(x) scale(x))

# LASSO differential response analysis ----

## Order of analysis ----

# 1. Fit LASSO model with cross validation to select the right lambda
# 2. Refit model with chosen lambda
# 4. RMSE and R^2


## Proteomics ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_SGLT2_proteomics_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_SGLT2_proteomics_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_SGLT2_proteomics_cv_model$lambda.min)

# variable selection
LASSO_DPP4_vs_SGLT2_proteomics_final_model_var_select_top_20 <- coef(LASSO_DPP4_vs_SGLT2_proteomics_final_model) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "proteomics") %>%
  mutate(coef = abs(s0)) %>%
  filter(coef != 0 & proteomics != "(Intercept)" & proteomics != "prehba1c") %>%
  arrange(desc(coef)) %>%
  select(proteomics) %>%
  slice_head(n = 20) %>%
  unlist()

# Predictions
LASSO_DPP4_vs_SGLT2_proteomics_final_predictions <- predict(LASSO_DPP4_vs_SGLT2_proteomics_final_model, s = LASSO_DPP4_vs_SGLT2_proteomics_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_SGLT2_proteomics_info <- data.frame(
  lambda = LASSO_DPP4_vs_SGLT2_proteomics_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_SGLT2_proteomics_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_SGLT2_proteomics_final_predictions)
) %>%
  rename("rsq" = s1)



### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_TZD_proteomics_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_TZD_proteomics_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_TZD_proteomics_cv_model$lambda.min)

# variable selection
LASSO_DPP4_vs_TZD_proteomics_final_model_var_select_top_20 <- coef(LASSO_DPP4_vs_TZD_proteomics_final_model) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "proteomics") %>%
  mutate(coef = abs(s0)) %>%
  filter(coef != 0 & proteomics != "(Intercept)" & proteomics != "prehba1c") %>%
  arrange(desc(coef)) %>%
  select(proteomics) %>%
  slice_head(n = 20) %>%
  unlist()

# Predictions
LASSO_DPP4_vs_TZD_proteomics_final_predictions <- predict(LASSO_DPP4_vs_TZD_proteomics_final_model, s = LASSO_DPP4_vs_TZD_proteomics_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_TZD_proteomics_info <- data.frame(
  lambda = LASSO_DPP4_vs_TZD_proteomics_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_TZD_proteomics_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_TZD_proteomics_final_predictions)
) %>%
  rename("rsq" = s1)




### SGLT2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_SGLT2_vs_TZD_proteomics_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_SGLT2_vs_TZD_proteomics_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_SGLT2_vs_TZD_proteomics_cv_model$lambda.min)

# variable selection
LASSO_SGLT2_vs_TZD_proteomics_final_model_var_select_top_20 <- coef(LASSO_SGLT2_vs_TZD_proteomics_final_model) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "proteomics") %>%
  mutate(coef = abs(s0)) %>%
  filter(coef != 0 & proteomics != "(Intercept)" & proteomics != "prehba1c") %>%
  arrange(desc(coef)) %>%
  select(proteomics) %>%
  slice_head(n = 20) %>%
  unlist()

# Predictions
LASSO_SGLT2_vs_TZD_proteomics_final_predictions <- predict(LASSO_SGLT2_vs_TZD_proteomics_final_model, s = LASSO_SGLT2_vs_TZD_proteomics_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_SGLT2_vs_TZD_proteomics_info <- data.frame(
  lambda = LASSO_SGLT2_vs_TZD_proteomics_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_SGLT2_vs_TZD_proteomics_final_predictions),
  rsq = rsq(y, LASSO_SGLT2_vs_TZD_proteomics_final_predictions)
) %>%
  rename("rsq" = s1)



## Proteomics (Top 20) ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(LASSO_DPP4_vs_SGLT2_proteomics_final_model_var_select_top_20, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_SGLT2_proteomics_top_20_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_SGLT2_proteomics_top_20_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_SGLT2_proteomics_top_20_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_SGLT2_proteomics_top_20_final_predictions <- predict(LASSO_DPP4_vs_SGLT2_proteomics_top_20_final_model, s = LASSO_DPP4_vs_SGLT2_proteomics_top_20_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_SGLT2_proteomics_top_20_info <- data.frame(
  lambda = LASSO_DPP4_vs_SGLT2_proteomics_top_20_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_SGLT2_proteomics_top_20_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_SGLT2_proteomics_top_20_final_predictions)
) %>%
  rename("rsq" = s1)



### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(LASSO_DPP4_vs_TZD_proteomics_final_model_var_select_top_20, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_TZD_proteomics_top_20_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_TZD_proteomics_top_20_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_TZD_proteomics_top_20_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_TZD_proteomics_top_20_final_predictions <- predict(LASSO_DPP4_vs_TZD_proteomics_top_20_final_model, s = LASSO_DPP4_vs_TZD_proteomics_top_20_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_TZD_proteomics_top_20_info <- data.frame(
  lambda = LASSO_DPP4_vs_TZD_proteomics_top_20_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_TZD_proteomics_top_20_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_TZD_proteomics_top_20_final_predictions)
) %>%
  rename("rsq" = s1)



### SGLT2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(LASSO_SGLT2_vs_TZD_proteomics_final_model_var_select_top_20, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_SGLT2_vs_TZD_proteomics_top_20_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_SGLT2_vs_TZD_proteomics_top_20_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_SGLT2_vs_TZD_proteomics_top_20_cv_model$lambda.min)

# Predictions
LASSO_SGLT2_vs_TZD_proteomics_top_20_final_predictions <- predict(LASSO_SGLT2_vs_TZD_proteomics_top_20_final_model, s = LASSO_SGLT2_vs_TZD_proteomics_top_20_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_SGLT2_vs_TZD_proteomics_top_20_info <- data.frame(
  lambda = LASSO_SGLT2_vs_TZD_proteomics_top_20_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_SGLT2_vs_TZD_proteomics_top_20_final_predictions),
  rsq = rsq(y, LASSO_SGLT2_vs_TZD_proteomics_top_20_final_predictions)
) %>%
  rename("rsq" = s1)



## Clinical features ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_SGLT2_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_SGLT2_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_SGLT2_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_SGLT2_clinical_features_final_predictions <- predict(LASSO_DPP4_vs_SGLT2_clinical_features_final_model, s = LASSO_DPP4_vs_SGLT2_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_SGLT2_clinical_features_info <- data.frame(
  lambda = LASSO_DPP4_vs_SGLT2_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_SGLT2_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_SGLT2_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)


### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_TZD_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_TZD_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_TZD_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_TZD_clinical_features_final_predictions <- predict(LASSO_DPP4_vs_TZD_clinical_features_final_model, s = LASSO_DPP4_vs_TZD_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_TZD_clinical_features_info <- data.frame(
  lambda = LASSO_DPP4_vs_TZD_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_TZD_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_TZD_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)



### SGLT2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_SGLT2_vs_TZD_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_SGLT2_vs_TZD_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_SGLT2_vs_TZD_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_SGLT2_vs_TZD_clinical_features_final_predictions <- predict(LASSO_SGLT2_vs_TZD_clinical_features_final_model, s = LASSO_SGLT2_vs_TZD_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_SGLT2_vs_TZD_clinical_features_info <- data.frame(
  lambda = LASSO_SGLT2_vs_TZD_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_SGLT2_vs_TZD_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_SGLT2_vs_TZD_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)



## Proteomics + Clinical features ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_final_predictions <- predict(LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_final_model, s = LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_info <- data.frame(
  lambda = LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)



### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()



# CV for alpha
LASSO_DPP4_vs_TZD_proteomics_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_TZD_proteomics_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_TZD_proteomics_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_TZD_proteomics_clinical_features_final_predictions <- predict(LASSO_DPP4_vs_TZD_proteomics_clinical_features_final_model, s = LASSO_DPP4_vs_TZD_proteomics_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_TZD_proteomics_clinical_features_info <- data.frame(
  lambda = LASSO_DPP4_vs_TZD_proteomics_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_TZD_proteomics_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_TZD_proteomics_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)



### SGLT2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()



# CV for alpha
LASSO_SGLT2_vs_TZD_proteomics_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_SGLT2_vs_TZD_proteomics_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_SGLT2_vs_TZD_proteomics_clinical_features_final_predictions <- predict(LASSO_SGLT2_vs_TZD_proteomics_clinical_features_final_model, s = LASSO_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_SGLT2_vs_TZD_proteomics_clinical_features_info <- data.frame(
  lambda = LASSO_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_SGLT2_vs_TZD_proteomics_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_SGLT2_vs_TZD_proteomics_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)



## Proteomics (Top 20) + Clinical features ----


### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(LASSO_DPP4_vs_SGLT2_proteomics_final_model_var_select_top_20, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_final_predictions <- predict(LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_final_model, s = LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_info <- data.frame(
  lambda = LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)



### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(LASSO_DPP4_vs_TZD_proteomics_final_model_var_select_top_20, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_final_predictions <- predict(LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_final_model, s = LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_info <- data.frame(
  lambda = LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)




### SGLT2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(LASSO_SGLT2_vs_TZD_proteomics_final_model_var_select_top_20, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for alpha
LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$lambda.min)

# Predictions
LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_final_predictions <- predict(LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_final_model, s = LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$lambda.min, newx = X) %>% unlist()

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_info <- data.frame(
  lambda = LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$lambda.min,
  rmse = Metrics::rmse(y, LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_final_predictions),
  rsq = rsq(y, LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_final_predictions)
) %>%
  rename("rsq" = s1)




## Results summary ----

# functions used
rsq_summary_function <- function(data, comparison, models) {
  return(
    data %>% 
      select(rsq) %>%
      mutate(
        comparison = comparison,
        models = models
      )
  )
}


### R2 summary ----
rsq_summary <- rsq_summary_function(LASSO_DPP4_vs_SGLT2_proteomics_info, "DPP4 vs SGLT2", "Proteomics") %>%
  rbind(
    rsq_summary_function(LASSO_DPP4_vs_TZD_proteomics_info, "DPP4 vs TZD", "Proteomics"),
    rsq_summary_function(LASSO_SGLT2_vs_TZD_proteomics_info, "SGLT2 vs TZD", "Proteomics"),
    rsq_summary_function(LASSO_DPP4_vs_SGLT2_proteomics_top_20_info, "DPP4 vs SGLT2", "Proteomics (Top 20)"),
    rsq_summary_function(LASSO_DPP4_vs_TZD_proteomics_top_20_info, "DPP4 vs TZD", "Proteomics (Top 20)"),
    rsq_summary_function(LASSO_SGLT2_vs_TZD_proteomics_top_20_info, "SGLT2 vs TZD", "Proteomics (Top 20)"),
    rsq_summary_function(LASSO_DPP4_vs_SGLT2_clinical_features_info, "DPP4 vs SGLT2", "Clinical features"),
    rsq_summary_function(LASSO_DPP4_vs_TZD_clinical_features_info, "DPP4 vs TZD", "Clinical features"),
    rsq_summary_function(LASSO_SGLT2_vs_TZD_clinical_features_info, "SGLT2 vs TZD", "Clinical features"),
    rsq_summary_function(LASSO_DPP4_vs_SGLT2_proteomics_clinical_features_info, "DPP4 vs SGLT2", "Proteomics + Clinical features"),
    rsq_summary_function(LASSO_DPP4_vs_TZD_proteomics_clinical_features_info, "DPP4 vs TZD", "Proteomics + Clinical features"),
    rsq_summary_function(LASSO_SGLT2_vs_TZD_proteomics_clinical_features_info, "SGLT2 vs TZD", "Proteomics + Clinical features"),
    rsq_summary_function(LASSO_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_info, "DPP4 vs SGLT2", "Proteomics (Top 20) + Clinical features"),
    rsq_summary_function(LASSO_DPP4_vs_TZD_proteomics_top_20_clinical_features_info, "DPP4 vs TZD", "Proteomics (Top 20) + Clinical features"),
    rsq_summary_function(LASSO_SGLT2_vs_TZD_proteomics_top_20_clinical_features_info, "SGLT2 vs TZD", "Proteomics (Top 20) + Clinical features")
  )


plot_rsq_summary <- rsq_summary %>%
  mutate(
    models = factor(models, levels = rev(c("Clinical features", "Proteomics", "Proteomics (Top 20)", "Proteomics + Clinical features", "Proteomics (Top 20) + Clinical features")))
  ) %>%
  ggplot(aes(y = comparison, x = rsq, colour = models)) +
  geom_point(position = position_dodge(width = 0.5)) +
  labs(x = "R2", y = "Differential effects models") +
  guides(color = guide_legend("Variable combinations", reverse = TRUE, ncol = 1)) +
  theme_classic() +
  facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_line(),
    panel.grid.minor.x = element_line(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

pdf("Plots/05.LASSO_R2.pdf", width = 7, height = 6)
plot_rsq_summary
dev.off()





