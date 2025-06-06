
# Initial setup ----

# Example paper of bartMachine
# https://doi.org/10.3390/rs16183379

# load libraries
options(java.parameters = "-Xmx4000m")
library(bartMachine)
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

# BART differential response analysis ----

## Order of analysis ----

# 1. Fit BART model with cross validation to select the right hyperparameters
# 2. Check the fit of the model (convergence and pred vs obs)
# 3. Variable importance (based on splits)
# 4. RMSE and R^2

# Parameters to try
num_tree_cvs <- seq(50, 200, by = 50)
k_cvs <- seq(1, 9, by = 2) # 2
nu = seq(1, 9, by = 2) # 2
q <- seq(0.1, 0.9, by = 0.2) # 0.2
nu_q_cvs <- apply(expand.grid(nu, q), 1, as.numeric)
nu_q_cvs <- lapply(seq_len(ncol(nu_q_cvs)), function(i) nu_q_cvs[, i])


## Proteomics ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_cv_model <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_cv_model, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_SGLT2_proteomics_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_SGLT2.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_SGLT2_proteomics_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_SGLT2_proteomics_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance.rds")
}

# variable importance top 20
BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance_top_20 <- setdiff(names(BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance$avg_var_props), "prehba1c")[1:20]


# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_SGLT2_proteomics_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_SGLT2_proteomics_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_SGLT2_proteomics_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_SGLT2_proteomics_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_SGLT2_proteomics_cv_model$k,
  nu = BART_DPP4_vs_SGLT2_proteomics_cv_model$nu,
  q = BART_DPP4_vs_SGLT2_proteomics_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_SGLT2_proteomics_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_SGLT2_proteomics_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_SGLT2_proteomics_info, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_info.rds")

### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_cv_model <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_TZD_proteomics_cv_model, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_TZD_proteomics_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_TZD_proteomics_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_TZD_proteomics_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model_var_importance.rds")
}

# variable importance top 20
BART_DPP4_vs_TZD_proteomics_cv_model_var_importance_top_20 <- setdiff(names(BART_DPP4_vs_TZD_proteomics_cv_model_var_importance$avg_var_props), "prehba1c")[1:20]

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_TZD_proteomics_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_TZD_proteomics_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_TZD_proteomics_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_TZD_proteomics_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_TZD_proteomics_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_TZD_proteomics_cv_model$k,
  nu = BART_DPP4_vs_TZD_proteomics_cv_model$nu,
  q = BART_DPP4_vs_TZD_proteomics_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_TZD_proteomics_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_TZD_proteomics_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_TZD_proteomics_info, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_info.rds")

### SGTL2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_cv_model <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_SGLT2_vs_TZD_proteomics_cv_model, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_SGLT2_vs_TZD_proteomics_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_SGLT2_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_SGLT2_vs_TZD_proteomics_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance <- investigate_var_importance(BART_SGLT2_vs_TZD_proteomics_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance.rds")
}

# variable importance top 20
BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance_top_20 <- setdiff(names(BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance$avg_var_props), "prehba1c")[1:20]

# variable selection
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model_var_select.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_cv_model_var_select <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model_var_select.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_cv_model_var_select <- var_selection_by_permute_cv(BART_SGLT2_vs_TZD_proteomics_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_cv_model_var_select, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_SGLT2_vs_TZD_proteomics_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_SGLT2_vs_TZD_proteomics_info <- data.frame(
  num_tree_cvs = BART_SGLT2_vs_TZD_proteomics_cv_model$num_trees,
  k_cvs = BART_SGLT2_vs_TZD_proteomics_cv_model$k,
  nu = BART_SGLT2_vs_TZD_proteomics_cv_model$nu,
  q = BART_SGLT2_vs_TZD_proteomics_cv_model$q,
  rmse = Metrics::rmse(y, BART_SGLT2_vs_TZD_proteomics_cv_model$y_hat_train),
  rsq = rsq(y, BART_SGLT2_vs_TZD_proteomics_cv_model$y_hat_train)
)
saveRDS(BART_SGLT2_vs_TZD_proteomics_info, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_info.rds")



## Proteomics (Top 20) ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance_top_20, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_SGLT2.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_SGLT2_proteomics_top_20_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model$k,
  nu = BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model$nu,
  q = BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_SGLT2_proteomics_top_20_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_info, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_info.rds")



### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(BART_DPP4_vs_TZD_proteomics_cv_model_var_importance_top_20, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_top_20_cv_model <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_top_20_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_cv_model, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_TZD_proteomics_top_20_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_TZD_proteomics_top_20_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_TZD_proteomics_top_20_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_TZD_proteomics_top_20_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_TZD_proteomics_top_20_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_TZD_proteomics_top_20_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_TZD_proteomics_top_20_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_TZD_proteomics_top_20_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_TZD_proteomics_top_20_cv_model$k,
  nu = BART_DPP4_vs_TZD_proteomics_top_20_cv_model$nu,
  q = BART_DPP4_vs_TZD_proteomics_top_20_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_TZD_proteomics_top_20_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_TZD_proteomics_top_20_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_info, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_info.rds")




### SGLT2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance_top_20, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_top_20_cv_model <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_top_20_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_SGLT2_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_importance.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_importance <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_importance <- investigate_var_importance(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_importance, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_select.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_select <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_select.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_select <- var_selection_by_permute_cv(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_select, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_SGLT2_vs_TZD_proteomics_top_20_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_SGLT2_vs_TZD_proteomics_top_20_info <- data.frame(
  num_tree_cvs = BART_SGLT2_vs_TZD_proteomics_top_20_cv_model$num_trees,
  k_cvs = BART_SGLT2_vs_TZD_proteomics_top_20_cv_model$k,
  nu = BART_SGLT2_vs_TZD_proteomics_top_20_cv_model$nu,
  q = BART_SGLT2_vs_TZD_proteomics_top_20_cv_model$q,
  rmse = Metrics::rmse(y, BART_SGLT2_vs_TZD_proteomics_top_20_cv_model$y_hat_train),
  rsq = rsq(y, BART_SGLT2_vs_TZD_proteomics_top_20_cv_model$y_hat_train)
)
saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_info, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_info.rds")













## Clinical features ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_clinical_features_cv_model <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_SGLT2_clinical_features_cv_model, file = "Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_SGLT2_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_SGLT2.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_SGLT2_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_SGLT2_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_SGLT2_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_SGLT2_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_SGLT2_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_SGLT2_clinical_features_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_SGLT2_clinical_features_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_SGLT2_clinical_features_cv_model$k,
  nu = BART_DPP4_vs_SGLT2_clinical_features_cv_model$nu,
  q = BART_DPP4_vs_SGLT2_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_SGLT2_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_SGLT2_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_SGLT2_clinical_features_info, file = "Interim_files/BART_DPP4_vs_SGLT2_clinical_features_info.rds")

### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model.rds")) {
  # load file
  BART_DPP4_vs_TZD_clinical_features_cv_model <- readRDS("Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_TZD_clinical_features_cv_model, file = "Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_TZD_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_TZD_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_TZD_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_TZD_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_clinical_features_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_TZD_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_TZD_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_clinical_features_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_TZD_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_TZD_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_TZD_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_TZD_clinical_features_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_TZD_clinical_features_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_TZD_clinical_features_cv_model$k,
  nu = BART_DPP4_vs_TZD_clinical_features_cv_model$nu,
  q = BART_DPP4_vs_TZD_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_TZD_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_TZD_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_TZD_clinical_features_info, file = "Interim_files/BART_DPP4_vs_TZD_clinical_features_info.rds")

### SGTL2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model.rds")) {
  # load file
  BART_SGLT2_vs_TZD_clinical_features_cv_model <- readRDS("Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_SGLT2_vs_TZD_clinical_features_cv_model, file = "Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_SGLT2_vs_TZD_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_SGLT2_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_SGLT2_vs_TZD_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_SGLT2_vs_TZD_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_SGLT2_vs_TZD_clinical_features_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_clinical_features_cv_model_var_importance, file = "Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_SGLT2_vs_TZD_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_SGLT2_vs_TZD_clinical_features_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_clinical_features_cv_model_var_select, file = "Interim_files/BART_SGLT2_vs_TZD_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_SGLT2_vs_TZD_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_SGLT2_vs_TZD_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_SGLT2_vs_TZD_clinical_features_info <- data.frame(
  num_tree_cvs = BART_SGLT2_vs_TZD_clinical_features_cv_model$num_trees,
  k_cvs = BART_SGLT2_vs_TZD_clinical_features_cv_model$k,
  nu = BART_SGLT2_vs_TZD_clinical_features_cv_model$nu,
  q = BART_SGLT2_vs_TZD_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_SGLT2_vs_TZD_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_SGLT2_vs_TZD_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_SGLT2_vs_TZD_clinical_features_info, file = "Interim_files/BART_SGLT2_vs_TZD_clinical_features_info.rds")



## Proteomics + Clinical features ----

### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_SGLT2.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_SGLT2_proteomics_clinical_features_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$k,
  nu = BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$nu,
  q = BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_SGLT2_proteomics_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_SGLT2_proteomics_clinical_features_info, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_clinical_features_info.rds")

### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_TZD_proteomics_clinical_features_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model$k,
  nu = BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model$nu,
  q = BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_TZD_proteomics_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_TZD_proteomics_clinical_features_info, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_clinical_features_info.rds")

### SGTL2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_SGLT2_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_importance, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_select, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_SGLT2_vs_TZD_proteomics_clinical_features_info <- data.frame(
  num_tree_cvs = BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$num_trees,
  k_cvs = BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$k,
  nu = BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$nu,
  q = BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_SGLT2_vs_TZD_proteomics_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_SGLT2_vs_TZD_proteomics_clinical_features_info, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_clinical_features_info.rds")


## Proteomics (Top 20) + Clinical features ----


### DPP4 vs SGLT2 ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(BART_DPP4_vs_SGLT2_proteomics_cv_model_var_importance_top_20, clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_SGLT2.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$k,
  nu = BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$nu,
  q = BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_info, file = "Interim_files/BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_info.rds")



### DPP4 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(BART_DPP4_vs_TZD_proteomics_cv_model_var_importance_top_20, clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_DPP4_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model)
  
  saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_info <- data.frame(
  num_tree_cvs = BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$num_trees,
  k_cvs = BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$k,
  nu = BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$nu,
  q = BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_info, file = "Interim_files/BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_info.rds")



### SGLT2 vs TZD ----
X = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(BART_SGLT2_vs_TZD_proteomics_cv_model_var_importance_top_20, clinical_features, "prehba1c")))
y = data %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()


# CV for hyperparameters (takes a bit of time)
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    serialize = TRUE
  )
  saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model.rds")
}

# Check model convergence
check_bart_error_assumptions(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model, hetero_plot = "yhats")
pdf("Plots/04.error_assumptions_SGLT2_TZD.pdf", width = 6, height = 8)
plot_y_vs_yhat(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model)
dev.off()


# variable importance
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance <- investigate_var_importance(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance.rds")
}

# variable selection
if (file.exists("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select.rds")) {
  # load file
  BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select <- readRDS("Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select.rds")
} else {
  # run analysis
  BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select <- var_selection_by_permute_cv(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model)
  
  saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_select.rds")
}

# # partial dependency plots (run this for only the important variables?)
# pd_plot(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model, which(c(variables_proteomics, "prehba1c") == names(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model_var_importance$avg_var_props[1])))


# Additional information (hyperparameters and error)
rsq <- function (x, y) cor(x, y) ^ 2

BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_info <- data.frame(
  num_tree_cvs = BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$num_trees,
  k_cvs = BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$k,
  nu = BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$nu,
  q = BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$q,
  rmse = Metrics::rmse(y, BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$y_hat_train),
  rsq = rsq(y, BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_cv_model$y_hat_train)
)
saveRDS(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_info, file = "Interim_files/BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_info.rds")




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
rsq_summary <- rsq_summary_function(BART_DPP4_vs_SGLT2_proteomics_info, "DPP4 vs SGLT2", "Proteomics") %>%
  rbind(
    rsq_summary_function(BART_DPP4_vs_TZD_proteomics_info, "DPP4 vs TZD", "Proteomics"),
    rsq_summary_function(BART_SGLT2_vs_TZD_proteomics_info, "SGLT2 vs TZD", "Proteomics"),
    rsq_summary_function(BART_DPP4_vs_SGLT2_proteomics_top_20_info, "DPP4 vs SGLT2", "Proteomics (Top 20)"),
    rsq_summary_function(BART_DPP4_vs_TZD_proteomics_top_20_info, "DPP4 vs TZD", "Proteomics (Top 20)"),
    rsq_summary_function(BART_SGLT2_vs_TZD_proteomics_top_20_info, "SGLT2 vs TZD", "Proteomics (Top 20)"),
    rsq_summary_function(BART_DPP4_vs_SGLT2_clinical_features_info, "DPP4 vs SGLT2", "Clinical features"),
    rsq_summary_function(BART_DPP4_vs_TZD_clinical_features_info, "DPP4 vs TZD", "Clinical features"),
    rsq_summary_function(BART_SGLT2_vs_TZD_clinical_features_info, "SGLT2 vs TZD", "Clinical features"),
    rsq_summary_function(BART_DPP4_vs_SGLT2_proteomics_clinical_features_info, "DPP4 vs SGLT2", "Proteomics + Clinical features"),
    rsq_summary_function(BART_DPP4_vs_TZD_proteomics_clinical_features_info, "DPP4 vs TZD", "Proteomics + Clinical features"),
    rsq_summary_function(BART_SGLT2_vs_TZD_proteomics_clinical_features_info, "SGLT2 vs TZD", "Proteomics + Clinical features"),
    rsq_summary_function(BART_DPP4_vs_SGLT2_proteomics_top_20_clinical_features_info, "DPP4 vs SGLT2", "Proteomics (Top 20) + Clinical features"),
    rsq_summary_function(BART_DPP4_vs_TZD_proteomics_top_20_clinical_features_info, "DPP4 vs TZD", "Proteomics (Top 20) + Clinical features"),
    rsq_summary_function(BART_SGLT2_vs_TZD_proteomics_top_20_clinical_features_info, "SGLT2 vs TZD", "Proteomics (Top 20) + Clinical features")
  )


plot_rsq_summary <- rsq_summary %>%
  mutate(
    models = factor(models, levels = rev(c("Clinical features", "Proteomics", "Proteomics (Top 20)", "Proteomics + Clinical features", "Proteomics (Top 20) + Clinical features")))
  ) %>%
  ggplot(aes(y = comparison, x = rsq, colour = models)) +
  geom_point(position = position_dodge(width = 0.5)) +
  labs(x = "R2", y = "Differential effects models") +
  guides(color = guide_legend("Variable combinations", reverse = TRUE, ncol = 1)) +
  theme_minimal() +
  theme(legend.position = "bottom")

