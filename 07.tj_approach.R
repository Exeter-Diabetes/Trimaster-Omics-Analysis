
# Initial setup ----

# load libraries
options(java.parameters = "-Xmx20000m")
library(bartMachine)
library(Metrics)
library(tidyverse)

# load nonlinvarsel
# install.packages('foreach', dep = T)
# url <- 'http://www.rob-mcculloch.org/chm/nonlinvarsel_0.0.1.9001.tar.gz'
# download.file(url, destfile = 'temp')
# install.packages('temp', repos = NULL, type='source')
library(nonlinvarsel)

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data_differential <- set_up_dataset(differential = TRUE) %>%
  mutate(sex = as.numeric(sex))

# proteomics variables
variables_proteomics <- data_differential %>% select(contains("proteom")) %>% colnames()
# clinical features
clinical_features <- c("prebmi", "sex", "agetx", "preegfr", "t2dmduration")

# standardisation
data_differential <- data_differential %>%
  mutate_at(c(variables_proteomics, clinical_features, "prehba1c"), function(x) scale(x))

# Approach ----

# 1. Fit BART model with many trees and CV (to find the best model with folds)
# 2. Check var importance and var selection
# 3. Use `nonlinvarsel` package to select top 20 union (back and forw)
# 4. Refit BART model with clinical features / proteomics / clinical features + proteomics (check in sample vs out of sample)


# Parameters to try
num_tree_cvs <- seq(50, 200, by = 50)
k_cvs <- seq(1, 9, by = 4) # 2
nu = seq(1, 9, by = 4) # 2
q <- seq(0.1, 0.9, by = 0.4) # 0.2
nu_q_cvs <- apply(expand.grid(nu, q), 1, as.numeric)
nu_q_cvs <- lapply(seq_len(ncol(nu_q_cvs)), function(i) nu_q_cvs[, i])


# DPP4 vs SGLT2 ----

## Default BART models ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# fit default model
if (file.exists("Interim_files/BART_default_DPP4_vs_SGLT2_model.rds")) {
  BART_default_DPP4_vs_SGLT2_model <- readRDS("Interim_files/BART_default_DPP4_vs_SGLT2_model.rds")
} else {
  BART_default_DPP4_vs_SGLT2_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    num_burn_in = 5000,
    num_iterations_after_burn_in = 5000,
    serialize = TRUE
  )
  
  saveRDS(BART_default_DPP4_vs_SGLT2_model, "Interim_files/BART_default_DPP4_vs_SGLT2_model.rds")
}

# 
# library(BART)
# 
# BART_model_sparse <- wbart(
#   x.train = as.data.frame(X),
#   y.train = y,
#   sparse = TRUE,
#   nskip = 10000,
#   ndpost = 10000
# )
# vars_selected <- BART_model_sparse$varprob.mean %>% as.data.frame() %>% rownames_to_column()
# View(vars_selected)
# 
# library(MASS)
# model <- lm.ridge(y ~ ., data = cbind(X, y = y))
# 
# 
# # CV for alpha
# LASSO_DPP4_vs_SGLT2_proteomics_cv_model <- cv.glmnet(as.matrix(X), y, alpha = 1, nfolds = 5)
# LASSO_DPP4_vs_SGLT2_proteomics_cv_model
# coef(LASSO_DPP4_vs_SGLT2_proteomics_cv_model)
# 
# 
# library(Boruta)
# 
# selected_features_list <- list()
# 
# for (i in 1:20) {
#   boruta_run <- Boruta(X, y, doTrace = 1)
#   selected <- getSelectedAttributes(boruta_run, withTentative = FALSE)
#   # selected <- names(boruta_run$finalDecision[boruta_run$finalDecision != "Rejected"])
#   selected_features_list[[i]] <- selected
# }
# 
# # Count frequency
# all_selected <- unlist(selected_features_list)
# table_selected <- sort(table(all_selected), decreasing = TRUE)
# 
# # Features selected in >=80% of runs
# stable_features <- names(table_selected[table_selected >= 16])
# stable_features
# 

# check error assumptions
# plot_convergence_diagnostics(BART_default_DPP4_vs_SGLT2_model)
# check_bart_error_assumptions(BART_default_DPP4_vs_SGLT2_model, hetero_plot = "ys")
# plot_y_vs_yhat(BART_default_DPP4_vs_SGLT2_model, credible_intervals = TRUE)

## Variable selection ----

# variable selection from nonlinvarsel
if (file.exists("Interim_files/BART_default_DPP4_vs_SGLT2_model_nonlinvarsel.rds")) {
  BART_default_DPP4_vs_SGLT2_model_nonlinvarsel <- readRDS("Interim_files/BART_default_DPP4_vs_SGLT2_model_nonlinvarsel.rds")
} else {
  BART_default_DPP4_vs_SGLT2_model_nonlinvarsel_forward <- vsf(X, BART_default_DPP4_vs_SGLT2_model$y_hat_train)
  BART_default_DPP4_vs_SGLT2_model_nonlinvarsel_backward <- vsb(X, BART_default_DPP4_vs_SGLT2_model$y_hat_train)
  
  # combine together
  BART_default_DPP4_vs_SGLT2_model_nonlinvarsel <- print(BART_default_DPP4_vs_SGLT2_model_nonlinvarsel_forward) %>%
    as.data.frame() %>%
    set_names(c("Forward")) %>%
    rownames_to_column() %>%
    filter(str_detect(rowname, "proteomic")) %>%
    left_join(
      print(BART_default_DPP4_vs_SGLT2_model_nonlinvarsel_backward) %>%
        as.data.frame() %>%
        set_names(c("Backward")) %>%
        rownames_to_column() %>%
        filter(str_detect(rowname, "proteomic"))
    ) %>%
    arrange(Forward) %>%
    mutate(Foward_pos = 1:n()) %>%
    arrange(Backward) %>%
    mutate(Backward_pos = 1:n()) %>%
    mutate(Average_pos = (Foward_pos + Backward_pos) / 2) %>%
    arrange(Average_pos)
  
  saveRDS(BART_default_DPP4_vs_SGLT2_model_nonlinvarsel, "Interim_files/BART_default_DPP4_vs_SGLT2_model_nonlinvarsel.rds")
}

## Refit models ----
current_dataset <- data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c)

current_dataset <- current_dataset %>%
  slice(sample(1:nrow(current_dataset)))

n_cv_groups = 5
current_dataset$cv_groups <- do.call(rbind.data.frame, split(1:nrow(current_dataset), cut(1:nrow(current_dataset), breaks = n_cv_groups, labels = FALSE))) %>% t() %>% as.data.frame() %>% remove_rownames() %>% gather() %>% distinct() %>% arrange(value) %>% select(key) %>% as.data.frame() %>% unlist() %>% as.vector()


### Clinical features ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(clinical_features, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(clinical_features, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_DPP4_vs_SGLT2_clinical_features_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_SGLT2_clinical_features_info", interim_info)



### Proteomics ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # variable selection
  vars_selected <- BART_default_DPP4_vs_SGLT2_model_nonlinvarsel %>%
    filter(Foward_pos < 11 | Backward_pos < 11) %>%
    select(rowname) %>%
    unlist()
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_DPP4_vs_SGLT2_proteomics_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_SGLT2_proteomics_info", interim_info)



### Proteomics + Clinical features ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # variable selection
  vars_selected <- BART_default_DPP4_vs_SGLT2_model_nonlinvarsel %>%
    filter(Foward_pos < 11 | Backward_pos < 11) %>%
    select(rowname) %>%
    unlist()
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, clinical_features, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, clinical_features, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_DPP4_vs_SGLT2_proteomics_clinical_features_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_SGLT2_proteomics_clinical_features_info", interim_info)




# DPP4 vs TZD ----

## Default BART models ----
X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# fit default model
if (file.exists("Interim_files/BART_default_DPP4_vs_TZD_model.rds")) {
  BART_default_DPP4_vs_TZD_model <- readRDS("Interim_files/BART_default_DPP4_vs_TZD_model.rds")
} else {
  BART_default_DPP4_vs_TZD_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    num_burn_in = 5000,
    num_iterations_after_burn_in = 5000,
    serialize = TRUE
  )
  
  saveRDS(BART_default_DPP4_vs_TZD_model, "Interim_files/BART_default_DPP4_vs_TZD_model.rds")
}


# check error assumptions
# plot_convergence_diagnostics(BART_default_DPP4_vs_TZD_model)
# check_bart_error_assumptions(BART_default_DPP4_vs_TZD_model, hetero_plot = "ys")
# plot_y_vs_yhat(BART_default_DPP4_vs_TZD_model, credible_intervals = TRUE)

## Variable selection ----

# variable selection from nonlinvarsel
if (file.exists("Interim_files/BART_default_DPP4_vs_TZD_model_nonlinvarsel.rds")) {
  BART_default_DPP4_vs_TZD_model_nonlinvarsel <- readRDS("Interim_files/BART_default_DPP4_vs_TZD_model_nonlinvarsel.rds")
} else {
  BART_default_DPP4_vs_TZD_model_nonlinvarsel_forward <- vsf(X, BART_default_DPP4_vs_TZD_model$y_hat_train)
  BART_default_DPP4_vs_TZD_model_nonlinvarsel_backward <- vsb(X, BART_default_DPP4_vs_TZD_model$y_hat_train)
  
  # combine together
  BART_default_DPP4_vs_TZD_model_nonlinvarsel <- print(BART_default_DPP4_vs_TZD_model_nonlinvarsel_forward) %>%
    as.data.frame() %>%
    set_names(c("Forward")) %>%
    rownames_to_column() %>%
    filter(str_detect(rowname, "proteomic")) %>%
    left_join(
      print(BART_default_DPP4_vs_TZD_model_nonlinvarsel_backward) %>%
        as.data.frame() %>%
        set_names(c("Backward")) %>%
        rownames_to_column() %>%
        filter(str_detect(rowname, "proteomic"))
    ) %>%
    arrange(Forward) %>%
    mutate(Foward_pos = 1:n()) %>%
    arrange(Backward) %>%
    mutate(Backward_pos = 1:n()) %>%
    mutate(Average_pos = (Foward_pos + Backward_pos) / 2) %>%
    arrange(Average_pos)
  
  saveRDS(BART_default_DPP4_vs_TZD_model_nonlinvarsel, "Interim_files/BART_default_DPP4_vs_TZD_model_nonlinvarsel.rds")
}

## Refit models ----
current_dataset <- data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c)

current_dataset <- current_dataset %>%
  slice(sample(1:nrow(current_dataset)))

n_cv_groups = 5
current_dataset$cv_groups <- do.call(rbind.data.frame, split(1:nrow(current_dataset), cut(1:nrow(current_dataset), breaks = n_cv_groups, labels = FALSE))) %>% t() %>% as.data.frame() %>% remove_rownames() %>% gather() %>% distinct() %>% arrange(value) %>% select(key) %>% as.data.frame() %>% unlist() %>% as.vector()


### Clinical features ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(clinical_features, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(clinical_features, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_DPP4_vs_TZD_clinical_features_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_TZD_clinical_features_info", interim_info)



### Proteomics ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # variable selection
  vars_selected <- BART_default_DPP4_vs_TZD_model_nonlinvarsel %>%
    filter(Foward_pos < 11 | Backward_pos < 11) %>%
    select(rowname) %>%
    unlist()
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_DPP4_vs_TZD_proteomics_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_TZD_proteomics_info", interim_info)



### Proteomics + Clinical features ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # variable selection
  vars_selected <- BART_default_DPP4_vs_TZD_model_nonlinvarsel %>%
    filter(Foward_pos < 11 | Backward_pos < 11) %>%
    select(rowname) %>%
    unlist()
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, clinical_features, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, clinical_features, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_DPP4_vs_TZD_proteomics_clinical_features_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_TZD_proteomics_clinical_features_info", interim_info)



# SGLT2 vs TZD ----

## Default BART models ----
X = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

# fit default model
if (file.exists("Interim_files/BART_default_SGLT2_vs_TZD_model.rds")) {
  BART_default_SGLT2_vs_TZD_model <- readRDS("Interim_files/BART_default_SGLT2_vs_TZD_model.rds")
} else {
  BART_default_SGLT2_vs_TZD_model <- bartMachineCV(
    X = X,
    y = y,
    num_tree_cvs = num_tree_cvs,
    k_cvs = k_cvs,
    nu_q_cvs = nu_q_cvs,
    k_folds = 5,
    num_burn_in = 5000,
    num_iterations_after_burn_in = 5000,
    serialize = TRUE
  )
  
  saveRDS(BART_default_SGLT2_vs_TZD_model, "Interim_files/BART_default_SGLT2_vs_TZD_model.rds")
}


# check error assumptions
# plot_convergence_diagnostics(BART_default_SGLT2_vs_TZD_model)
# check_bart_error_assumptions(BART_default_SGLT2_vs_TZD_model, hetero_plot = "ys")
# plot_y_vs_yhat(BART_default_SGLT2_vs_TZD_model, credible_intervals = TRUE)

## Variable selection ----

# variable selection from nonlinvarsel
if (file.exists("Interim_files/BART_default_SGLT2_vs_TZD_model_nonlinvarsel.rds")) {
  BART_default_SGLT2_vs_TZD_model_nonlinvarsel <- readRDS("Interim_files/BART_default_SGLT2_vs_TZD_model_nonlinvarsel.rds")
} else {
  BART_default_SGLT2_vs_TZD_model_nonlinvarsel_forward <- vsf(X, BART_default_SGLT2_vs_TZD_model$y_hat_train)
  BART_default_SGLT2_vs_TZD_model_nonlinvarsel_backward <- vsb(X, BART_default_SGLT2_vs_TZD_model$y_hat_train)
  
  # combine together
  BART_default_SGLT2_vs_TZD_model_nonlinvarsel <- print(BART_default_SGLT2_vs_TZD_model_nonlinvarsel_forward) %>%
    as.data.frame() %>%
    set_names(c("Forward")) %>%
    rownames_to_column() %>%
    filter(str_detect(rowname, "proteomic")) %>%
    left_join(
      print(BART_default_SGLT2_vs_TZD_model_nonlinvarsel_backward) %>%
        as.data.frame() %>%
        set_names(c("Backward")) %>%
        rownames_to_column() %>%
        filter(str_detect(rowname, "proteomic"))
    ) %>%
    arrange(Forward) %>%
    mutate(Foward_pos = 1:n()) %>%
    arrange(Backward) %>%
    mutate(Backward_pos = 1:n()) %>%
    mutate(Average_pos = (Foward_pos + Backward_pos) / 2) %>%
    arrange(Average_pos)
  
  saveRDS(BART_default_SGLT2_vs_TZD_model_nonlinvarsel, "Interim_files/BART_default_SGLT2_vs_TZD_model_nonlinvarsel.rds")
}

## Refit models ----
current_dataset <- data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c)

current_dataset <- current_dataset %>%
  slice(sample(1:nrow(current_dataset)))

n_cv_groups = 5
current_dataset$cv_groups <- do.call(rbind.data.frame, split(1:nrow(current_dataset), cut(1:nrow(current_dataset), breaks = n_cv_groups, labels = FALSE))) %>% t() %>% as.data.frame() %>% remove_rownames() %>% gather() %>% distinct() %>% arrange(value) %>% select(key) %>% as.data.frame() %>% unlist() %>% as.vector()


### Clinical features ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(clinical_features, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(clinical_features, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_SGLT2_vs_TZD_clinical_features_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_SGLT2_vs_TZD_clinical_features_info", interim_info)



### Proteomics ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # variable selection
  vars_selected <- BART_default_SGLT2_vs_TZD_model_nonlinvarsel %>%
    filter(Foward_pos < 11 | Backward_pos < 11) %>%
    select(rowname) %>%
    unlist()
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_SGLT2_vs_TZD_proteomics_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_SGLT2_vs_TZD_proteomics_info", interim_info)



### Proteomics + Clinical features ----
interim_info <- NULL
for (cv_it in 1:5) {
  
  # variable selection
  vars_selected <- BART_default_SGLT2_vs_TZD_model_nonlinvarsel %>%
    filter(Foward_pos < 11 | Backward_pos < 11) %>%
    select(rowname) %>%
    unlist()
  
  # training
  X = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, clinical_features, "prehba1c")))
  y = current_dataset %>%
    filter(cv_groups != paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  # validation
  X_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(all_of(c(vars_selected, clinical_features, "prehba1c")))
  y_val = current_dataset %>%
    filter(cv_groups == paste0("V", cv_it)) %>%
    select(benefit) %>%
    unlist()
  
  # fit default model
  file_name <- paste0("BART_SGLT2_vs_TZD_proteomics_clinical_features_model_v", cv_it)
  if (file.exists(paste0("Interim_files/", file_name, ".rds"))) {
    interim_model <- readRDS(paste0("Interim_files/", file_name, ".rds"))
    assign(file_name, interim_model)
  } else {
    interim_model <- bartMachineCV(
      X = X,
      y = y,
      num_tree_cvs = num_tree_cvs,
      k_cvs = k_cvs,
      nu_q_cvs = nu_q_cvs,
      k_folds = 5,
      num_burn_in = 5000,
      num_iterations_after_burn_in = 5000,
      serialize = TRUE
    )
    assign(file_name, interim_model)
    
    saveRDS(interim_model, paste0("Interim_files/", file_name, ".rds"))
  }
  
  # check error assumptions
  # plot_convergence_diagnostics(interim_model)
  # check_bart_error_assumptions(interim_model, hetero_plot = "ys")
  # plot_y_vs_yhat(interim_model, credible_intervals = TRUE)
  
  # prediction
  interim_prediction <- predict(interim_model, new_data = X_val)
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, interim_model$y_hat_train),
      rsq_in = rsq(y, interim_model$y_hat_train),
      rmse_out = Metrics::rmse(y_val, interim_prediction),
      rsq_out = rsq(y_val, interim_prediction),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_SGLT2_vs_TZD_proteomics_clinical_features_info", interim_info)









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
rsq_summary <- rsq_summary_function(BART_DPP4_vs_SGLT2_clinical_features_info, "DPP4 vs SGLT2", "Clinical features") %>%
  rbind(
    rsq_summary_function(BART_DPP4_vs_TZD_clinical_features_info, "DPP4 vs TZD", "Clinical features"),
    rsq_summary_function(BART_SGLT2_vs_TZD_clinical_features_info, "SGLT2 vs TZD", "Clinical features"),
    rsq_summary_function(BART_DPP4_vs_SGLT2_proteomics_info, "DPP4 vs SGLT2", "Proteomics"),
    rsq_summary_function(BART_DPP4_vs_TZD_proteomics_info, "DPP4 vs TZD", "Proteomics"),
    rsq_summary_function(BART_SGLT2_vs_TZD_proteomics_info, "SGLT2 vs TZD", "Proteomics"),
    rsq_summary_function(BART_DPP4_vs_SGLT2_proteomics_clinical_features_info, "DPP4 vs SGLT2", "Proteomics + Clinical features"),
    rsq_summary_function(BART_DPP4_vs_TZD_proteomics_clinical_features_info, "DPP4 vs TZD", "Proteomics + Clinical features"),
    rsq_summary_function(BART_SGLT2_vs_TZD_proteomics_clinical_features_info, "SGLT2 vs TZD", "Proteomics + Clinical features")
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
    BART_default_DPP4_vs_SGLT2_model_nonlinvarsel %>%
      filter(Foward_pos < 11 | Backward_pos < 11) %>%
      select(Names = rowname) %>%
      mutate(DPP4_vs_SGLT2 = 1), by = c("Names")
  ) %>%
  left_join(
    BART_default_DPP4_vs_TZD_model_nonlinvarsel %>%
      filter(Foward_pos < 11 | Backward_pos < 11) %>%
      select(Names = rowname) %>%
      mutate(DPP4_vs_TZD = 1), by = c("Names")
  ) %>%
  left_join(
    BART_default_SGLT2_vs_TZD_model_nonlinvarsel %>%
      filter(Foward_pos < 11 | Backward_pos < 11) %>%
      select(Names = rowname) %>%
      mutate(SGLT2_vs_TZD = 1), by = c("Names")
  ) %>%
  mutate(row = 1:n()) %>%
  group_by(row) %>%
  mutate(Included = sum(c(DPP4_vs_SGLT2), na.rm = TRUE)) %>%
  # mutate(Included = sum(c(DPP4_vs_SGLT2, DPP4_vs_TZD, SGLT2_vs_TZD), na.rm = TRUE)) %>%
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

View(proteomics_summary)





