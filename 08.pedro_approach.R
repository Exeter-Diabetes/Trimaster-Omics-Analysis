
# Initial setup ----

# load libraries
library(Metrics)
library(tidyverse)
library(Boruta)

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

# 1. Fit the Boruta method on all proteomics and clinical features 100 times.
# 2. Take the features confirmed at least 80% of runs
# 3. Divide the data into 5 folds
# 4. Fit the default LM model in a 5 fold cross-validation and test the in-sample and out-sample R2.
# 4.1. Clinical features
# 4.2. Selected Proteomics
# 4.3. Clinical features + Proteomics


# DPP4 vs SGLT2 ----

X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

## Variable selection: Boruta ----
# number of iterations
boruta_it <- 1000

if (file.exists("Interim_files/DPP4_vs_SGLT2_table_selected.rds")) {
  DPP4_vs_SGLT2_table_selected.rds <- readRDS("Interim_files/DPP4_vs_SGLT2_table_selected.rds")
} else {
  
  selected_features_list <- list()
  
  for (i in 1:boruta_it) {
    print(i)
    boruta_run <- Boruta(X, y, doTrace = 0)
    selected <- getSelectedAttributes(boruta_run, withTentative = FALSE)
    selected_features_list[[i]] <- selected
  }
  
  # Count frequency
  all_selected <- unlist(selected_features_list)
  DPP4_vs_SGLT2_table_selected <- sort(table(all_selected), decreasing = TRUE)
  
  saveRDS(DPP4_vs_SGLT2_table_selected, file = "Interim_files/DPP4_vs_SGLT2_table_selected.rds")
  
}
# Features selected in >=80% of runs
DPP4_vs_SGLT2_variables_selected <- names(DPP4_vs_SGLT2_table_selected[DPP4_vs_SGLT2_table_selected >= round(0.8*boruta_it)])



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
for (cv_it in 1:n_cv_groups) {
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
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
  vars_selected <- grep("^proteomics", DPP4_vs_SGLT2_variables_selected, value = TRUE)
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
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
  vars_selected <- grep("^proteomics", DPP4_vs_SGLT2_variables_selected, value = TRUE)
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_SGLT2_proteomics_clinical_features_info", interim_info)




# DPP4 vs TZD ----

X = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_DPP4 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

## Variable selection: Boruta ----

if (file.exists("Interim_files/DPP4_vs_TZD_table_selected.rds")) {
  DPP4_vs_TZD_table_selected.rds <- readRDS("Interim_files/DPP4_vs_TZD_table_selected.rds")
} else {
  
  selected_features_list <- list()
  
  for (i in 1:boruta_it) {
    print(i)
    boruta_run <- Boruta(X, y, doTrace = 0)
    selected <- getSelectedAttributes(boruta_run, withTentative = FALSE)
    selected_features_list[[i]] <- selected
  }
  
  # Count frequency
  all_selected <- unlist(selected_features_list)
  DPP4_vs_TZD_table_selected <- sort(table(all_selected), decreasing = TRUE)
  
  saveRDS(DPP4_vs_TZD_table_selected, "Interim_files/DPP4_vs_TZD_table_selected.rds")
  
}
# Features selected in >=80% of runs
DPP4_vs_TZD_variables_selected <- names(DPP4_vs_TZD_table_selected[DPP4_vs_TZD_table_selected >= round(0.8*boruta_it)])



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
for (cv_it in 1:n_cv_groups) {
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
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
  vars_selected <- grep("^proteomics", DPP4_vs_TZD_variables_selected, value = TRUE)
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
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
  vars_selected <- grep("^proteomics", DPP4_vs_TZD_variables_selected, value = TRUE)
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
      cv = cv_it
    )
  )
  
}

# change name
assign("BART_DPP4_vs_TZD_proteomics_clinical_features_info", interim_info)




# SGLT2 vs TZD ----

X = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(all_of(c(variables_proteomics, clinical_features, "prehba1c")))
y = data_differential %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  drop_na(benefit, all_of(c(variables_proteomics, clinical_features)), prehba1c) %>%
  select(benefit) %>%
  unlist()

## Variable selection: Boruta ----

if (file.exists("Interim_files/SGLT2_vs_TZD_table_selected.rds")) {
  SGLT2_vs_TZD_table_selected.rds <- readRDS("Interim_files/SGLT2_vs_TZD_table_selected.rds")
} else {
  
  selected_features_list <- list()
  
  for (i in 1:boruta_it) {
    print(i)
    boruta_run <- Boruta(X, y, doTrace = 0)
    selected <- getSelectedAttributes(boruta_run, withTentative = FALSE)
    selected_features_list[[i]] <- selected
  }
  
  # Count frequency
  all_selected <- unlist(selected_features_list)
  SGLT2_vs_TZD_table_selected <- sort(table(all_selected), decreasing = TRUE)
  
  saveRDS(SGLT2_vs_TZD_table_selected, "Interim_files/SGLT2_vs_TZD_table_selected.rds")
  
}
# Features selected in >=80% of runs
SGLT2_vs_TZD_variables_selected <- names(SGLT2_vs_TZD_table_selected[SGLT2_vs_TZD_table_selected >= round(0.8*boruta_it)])



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
for (cv_it in 1:n_cv_groups) {
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
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
  vars_selected <- grep("^proteomics", SGLT2_vs_TZD_variables_selected, value = TRUE)
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
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
  vars_selected <- grep("^proteomics", SGLT2_vs_TZD_variables_selected, value = TRUE)
  
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
  linear_model <- lm(y ~ ., data = cbind(X, y = y))
  
  # Additional information (hyperparameters and error)
  rsq <- function (x, y) cor(x, y) ^ 2
  
  interim_info <- rbind(
    interim_info,
    data.frame(
      rmse_in = Metrics::rmse(y, predict(linear_model, newdata = X)),
      rsq_in = rsq(y, predict(linear_model, newdata = X)),
      rmse_out = Metrics::rmse(y_val, predict(linear_model, newdata = X_val)),
      rsq_out = rsq(y_val, predict(linear_model, newdata = X_val)),
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
      select(rsq_out) %>%
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
  ggplot(aes(y = comparison, x = rsq_out, colour = models)) +
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

# proteomics_summary <- variables_proteomics %>%
#   as.data.frame() %>%
#   set_names(c("Names")) %>%
#   left_join(
#     BART_default_DPP4_vs_SGLT2_model_nonlinvarsel %>%
#       filter(Foward_pos < 11 | Backward_pos < 11) %>%
#       select(Names = rowname) %>%
#       mutate(DPP4_vs_SGLT2 = 1), by = c("Names")
#   ) %>%
#   left_join(
#     BART_default_DPP4_vs_TZD_model_nonlinvarsel %>%
#       filter(Foward_pos < 11 | Backward_pos < 11) %>%
#       select(Names = rowname) %>%
#       mutate(DPP4_vs_TZD = 1), by = c("Names")
#   ) %>%
#   left_join(
#     BART_default_SGLT2_vs_TZD_model_nonlinvarsel %>%
#       filter(Foward_pos < 11 | Backward_pos < 11) %>%
#       select(Names = rowname) %>%
#       mutate(SGLT2_vs_TZD = 1), by = c("Names")
#   ) %>%
#   mutate(row = 1:n()) %>%
#   group_by(row) %>%
#   mutate(Included = sum(c(DPP4_vs_SGLT2), na.rm = TRUE)) %>%
#   # mutate(Included = sum(c(DPP4_vs_SGLT2, DPP4_vs_TZD, SGLT2_vs_TZD), na.rm = TRUE)) %>%
#   ungroup() %>%
#   filter(Included > 0) %>%
#   select(-c(row, Included)) %>%
#   as.data.frame() %>%
#   mutate(Names = toupper(gsub("proteomics_", "", Names))) %>%
#   left_join(
#     human_protein_atlas %>%
#       rename("Names" = "Gene", "Gene description" = "Gene.description"), by = c("Names")
#   ) %>%
#   relocate("Gene description", .after = "Names") %>%
#   left_join(
#     unitprotkb_function %>%
#       rename("Uniprot" = "Entry") %>%
#       select(-c(`Gene Names`, `Entry Name`, `Protein names`)), by = c("Uniprot")
#   ) %>%
#   select(-c(`Uniprot`))
# 
# View(proteomics_summary)





