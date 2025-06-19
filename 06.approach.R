
# Initial setup ----

# load libraries
library(glmnet)
library(Metrics)
library(tidyverse)

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data_response <- set_up_dataset(differential = FALSE) %>%
  mutate(sex = as.numeric(sex))
data_differential <- set_up_dataset(differential = TRUE)

# proteomics variables
variables_proteomics <- data_response %>% select(contains("proteom")) %>% colnames()
# clinical features
clinical_features <- c("prebmi", "sex", "agetx", "preegfr", "t2dmduration")

# standardisation
data_response <- data_response %>%
  mutate_at(c(variables_proteomics, "prebmi", "agetx", "preegfr", "t2dmduration", "prehba1c"), function(x) scale(x))


# Approach ----

# 1. LASSO with proteomics + clinical features
# 2. Differential response using clinical features / proteomics / clinical features + proteomics (top 20)



## LASSO: Proteomics + Clinical features ----

### DPP4 ----
X = data_response %>%
  filter(drugclass == "DPP4") %>%
  drop_na(posthba1cfinal, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data_response %>%
  filter(drugclass == "DPP4") %>%
  drop_na(posthba1cfinal, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(posthba1cfinal) %>%
  unlist()

# CV for alpha
LASSO_DPP4_proteomics_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 10)

# Refit model with chosen alpha
LASSO_DPP4_proteomics_clinical_features_cv_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_DPP4_proteomics_clinical_features_cv_model$lambda.min)

# variable selection
LASSO_DPP4_proteomics_clinical_features_var_selection <- coef(LASSO_DPP4_proteomics_clinical_features_cv_final_model) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "proteomics") %>%
  mutate(coef = abs(s0)) %>%
  filter(coef != 0 & proteomics != "(Intercept)" & !(proteomics %in% c("prehba1c", clinical_features))) %>%
  arrange(desc(coef)) %>%
  select(proteomics) %>%
  slice_head(n = 20) %>%
  unlist()


### SGLT2 ----
X = data_response %>%
  filter(drugclass == "SGLT2") %>%
  drop_na(posthba1cfinal, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data_response %>%
  filter(drugclass == "SGLT2") %>%
  drop_na(posthba1cfinal, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(posthba1cfinal) %>%
  unlist()

# CV for alpha
LASSO_SGLT2_proteomics_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_SGLT2_proteomics_clinical_features_cv_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_SGLT2_proteomics_clinical_features_cv_model$lambda.min)

# variable selection
LASSO_SGLT2_proteomics_clinical_features_var_selection <- coef(LASSO_SGLT2_proteomics_clinical_features_cv_final_model) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "proteomics") %>%
  mutate(coef = abs(s0)) %>%
  filter(coef != 0 & proteomics != "(Intercept)" & !(proteomics %in% c("prehba1c", clinical_features))) %>%
  arrange(desc(coef)) %>%
  select(proteomics) %>%
  slice_head(n = 20) %>%
  unlist()



### TZD ----
X = data_response %>%
  filter(drugclass == "TZD") %>%
  drop_na(posthba1cfinal, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c"))) %>%
  as.matrix()
y = data_response %>%
  filter(drugclass == "TZD") %>%
  drop_na(posthba1cfinal, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(posthba1cfinal) %>%
  unlist()

# CV for alpha
LASSO_TZD_proteomics_clinical_features_cv_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)

# Refit model with chosen alpha
LASSO_TZD_proteomics_clinical_features_cv_final_model <- glmnet(X, y, alpha = 1, lambda = LASSO_TZD_proteomics_clinical_features_cv_model$lambda.min)

# variable selection
LASSO_TZD_proteomics_clinical_features_var_selection <- coef(LASSO_TZD_proteomics_clinical_features_cv_final_model) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "proteomics") %>%
  mutate(coef = abs(s0)) %>%
  filter(coef != 0 & proteomics != "(Intercept)" & !(proteomics %in% c("prehba1c", clinical_features))) %>%
  arrange(desc(coef)) %>%
  select(proteomics) %>%
  slice_head(n = 20) %>%
  unlist()


## Differential effects ----

### Clinical features ----

#### DPP4 vs SGLT2 ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_DPP4_vs_SGLT2_clinical_features_model <- lm(y ~ ., data = cbind(X, y = y))

LM_DPP4_vs_SGLT2_clinical_features_predictions <- predict(LM_DPP4_vs_SGLT2_clinical_features_model, newdata = X)
  
# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_DPP4_vs_SGLT2_clinical_features_information <- data.frame(
  rmse = Metrics::rmse(y, LM_DPP4_vs_SGLT2_clinical_features_predictions),
  rsq = rsq(y, LM_DPP4_vs_SGLT2_clinical_features_predictions)
)


#### DPP4 vs TZD ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_DPP4_vs_TZD_clinical_features_model <- lm(y ~ ., data = cbind(X, y = y))

LM_DPP4_vs_TZD_clinical_features_predictions <- predict(LM_DPP4_vs_TZD_clinical_features_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_DPP4_vs_TZD_clinical_features_information <- data.frame(
  rmse = Metrics::rmse(y, LM_DPP4_vs_TZD_clinical_features_predictions),
  rsq = rsq(y, LM_DPP4_vs_TZD_clinical_features_predictions)
)


#### SGLT2 vs TZD ----
X = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_SGLT2_vs_TZD_clinical_features_model <- lm(y ~ ., data = cbind(X, y = y))

LM_SGLT2_vs_TZD_clinical_features_predictions <- predict(LM_SGLT2_vs_TZD_clinical_features_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_SGLT2_vs_TZD_clinical_features_information <- data.frame(
  rmse = Metrics::rmse(y, LM_SGLT2_vs_TZD_clinical_features_predictions),
  rsq = rsq(y, LM_SGLT2_vs_TZD_clinical_features_predictions)
)



### Proteomics ----

#### DPP4 vs SGLT2 ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(unique(c(LASSO_DPP4_proteomics_clinical_features_var_selection, LASSO_SGLT2_proteomics_clinical_features_var_selection, "prehba1c"))))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_DPP4_vs_SGLT2_proteomics_model <- lm(y ~ ., data = cbind(X, y = y))

LM_DPP4_vs_SGLT2_proteomics_predictions <- predict(LM_DPP4_vs_SGLT2_proteomics_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_DPP4_vs_SGLT2_proteomics_information <- data.frame(
  rmse = Metrics::rmse(y, LM_DPP4_vs_SGLT2_proteomics_predictions),
  rsq = rsq(y, LM_DPP4_vs_SGLT2_proteomics_predictions)
)


#### DPP4 vs TZD ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(unique(c(LASSO_DPP4_proteomics_clinical_features_var_selection, LASSO_TZD_proteomics_clinical_features_var_selection, "prehba1c"))))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_DPP4_vs_TZD_proteomics_model <- lm(y ~ ., data = cbind(X, y = y))

LM_DPP4_vs_TZD_proteomics_predictions <- predict(LM_DPP4_vs_TZD_proteomics_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_DPP4_vs_TZD_proteomics_information <- data.frame(
  rmse = Metrics::rmse(y, LM_DPP4_vs_TZD_proteomics_predictions),
  rsq = rsq(y, LM_DPP4_vs_TZD_proteomics_predictions)
)


#### SGLT2 vs TZD ----
X = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(unique(c(LASSO_SGLT2_proteomics_clinical_features_var_selection, LASSO_TZD_proteomics_clinical_features_var_selection, "prehba1c"))))
y = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_SGLT2_vs_TZD_proteomics_model <- lm(y ~ ., data = cbind(X, y = y))

LM_SGLT2_vs_TZD_proteomics_predictions <- predict(LM_SGLT2_vs_TZD_proteomics_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_SGLT2_vs_TZD_proteomics_information <- data.frame(
  rmse = Metrics::rmse(y, LM_SGLT2_vs_TZD_proteomics_predictions),
  rsq = rsq(y, LM_SGLT2_vs_TZD_proteomics_predictions)
)



### Proteomics + clinical features ----

#### DPP4 vs SGLT2 ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(unique(c(LASSO_DPP4_proteomics_clinical_features_var_selection, LASSO_SGLT2_proteomics_clinical_features_var_selection, clinical_features, "prehba1c"))))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_DPP4_vs_SGLT2_proteomics_clinical_features_model <- lm(y ~ ., data = cbind(X, y = y))

LM_DPP4_vs_SGLT2_proteomics_clinical_features_predictions <- predict(LM_DPP4_vs_SGLT2_proteomics_clinical_features_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_DPP4_vs_SGLT2_proteomics_clinical_features_information <- data.frame(
  rmse = Metrics::rmse(y, LM_DPP4_vs_SGLT2_proteomics_clinical_features_predictions),
  rsq = rsq(y, LM_DPP4_vs_SGLT2_proteomics_clinical_features_predictions)
)


#### DPP4 vs TZD ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(unique(c(LASSO_DPP4_proteomics_clinical_features_var_selection, LASSO_TZD_proteomics_clinical_features_var_selection, clinical_features, "prehba1c"))))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_DPP4_vs_TZD_proteomics_clinical_features_model <- lm(y ~ ., data = cbind(X, y = y))

LM_DPP4_vs_TZD_proteomics_clinical_features_predictions <- predict(LM_DPP4_vs_TZD_proteomics_clinical_features_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_DPP4_vs_TZD_proteomics_clinical_features_information <- data.frame(
  rmse = Metrics::rmse(y, LM_DPP4_vs_TZD_proteomics_clinical_features_predictions),
  rsq = rsq(y, LM_DPP4_vs_TZD_proteomics_clinical_features_predictions)
)


#### SGLT2 vs TZD ----
X = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(unique(c(LASSO_SGLT2_proteomics_clinical_features_var_selection, LASSO_TZD_proteomics_clinical_features_var_selection, clinical_features, "prehba1c"))))
y = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

LM_SGLT2_vs_TZD_proteomics_clinical_features_model <- lm(y ~ ., data = cbind(X, y = y))

LM_SGLT2_vs_TZD_proteomics_clinical_features_predictions <- predict(LM_SGLT2_vs_TZD_proteomics_clinical_features_model, newdata = X)

# Additional information
rsq <- function (x, y) cor(x, y) ^ 2

LM_SGLT2_vs_TZD_proteomics_clinical_features_information <- data.frame(
  rmse = Metrics::rmse(y, LM_SGLT2_vs_TZD_proteomics_clinical_features_predictions),
  rsq = rsq(y, LM_SGLT2_vs_TZD_proteomics_clinical_features_predictions)
)




# Results summary ----

# functions used
rsq_summary_function <- function(data, comparison, models) {
  return(
    data %>% 
      select(rsq) %>%
      mutate(comparison = comparison, models = models)
  )
}

### R2 summary ----
rsq_summary <- rsq_summary_function(LM_DPP4_vs_SGLT2_clinical_features_information, "DPP4 vs SGLT2", "Clinical features") %>%
  rbind(
    rsq_summary_function(LM_DPP4_vs_TZD_clinical_features_information, "DPP4 vs TZD", "Clinical features"),
    rsq_summary_function(LM_SGLT2_vs_TZD_clinical_features_information, "SGLT2 vs TZD", "Clinical features"),
    rsq_summary_function(LM_DPP4_vs_SGLT2_proteomics_information, "DPP4 vs SGLT2", "Proteomics"),
    rsq_summary_function(LM_DPP4_vs_TZD_proteomics_information, "DPP4 vs TZD", "Proteomics"),
    rsq_summary_function(LM_SGLT2_vs_TZD_proteomics_information, "SGLT2 vs TZD", "Proteomics"),
    rsq_summary_function(LM_DPP4_vs_SGLT2_proteomics_clinical_features_information, "DPP4 vs SGLT2", "Proteomics + Clinical features"),
    rsq_summary_function(LM_DPP4_vs_TZD_proteomics_clinical_features_information, "DPP4 vs TZD", "Proteomics + Clinical features"),
    rsq_summary_function(LM_SGLT2_vs_TZD_proteomics_clinical_features_information, "SGLT2 vs TZD", "Proteomics + Clinical features")
  )


plot_rsq_summary <- rsq_summary %>%
  mutate(
    models = factor(models, levels = rev(c("Clinical features", "Proteomics", "Proteomics + Clinical features")))
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


### Proteomics summary ----
human_protein_atlas <- read.delim("Proteomics Descriptions/human_protein_atlas.tsv")
unitprotkb_function <- readxl::read_excel("Proteomics Descriptions/uniprotkb_function.xlsx")

proteomics_summary <- variables_proteomics %>%
  as.data.frame() %>%
  set_names(c("Names")) %>%
  left_join(
    LASSO_SGLT2_proteomics_clinical_features_var_selection %>%
      as.data.frame() %>%
      set_names(c("Names")) %>%
      mutate(SGLT2 = 1), by = c("Names")
  ) %>%
  left_join(
    LASSO_DPP4_proteomics_clinical_features_var_selection %>%
      as.data.frame() %>%
      set_names(c("Names")) %>%
      mutate(DPP4 = 1), by = c("Names")
  ) %>%
  left_join(
    LASSO_TZD_proteomics_clinical_features_var_selection %>%
      as.data.frame() %>%
      set_names(c("Names")) %>%
      mutate(TZD = 1), by = c("Names")
  ) %>%
  mutate(row = 1:n()) %>%
  group_by(row) %>%
  mutate(Included = sum(c(SGLT2, DPP4, TZD), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Included > 0) %>%
  select(-c(row, Included)) %>%
  as.data.frame() %>%
  mutate(Names = toupper(gsub("proteomics_", "", Names))) %>%
  left_join(
    human_protein_atlas %>%
      rename("Names" = "Gene", "Gene description" = "Gene.description"), by = c("Names")
  ) %>%
  relocate("Gene description", .after = "Names") %>%
  left_join(
    unitprotkb_function %>%
      rename("Uniprot" = "Entry") %>%
      select(-c(`Gene Names`, `Entry Name`, `Protein names`)), by = c("Uniprot")
  ) %>%
  select(-c(`Uniprot`))

proteomics_summary


