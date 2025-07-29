
# Initial setup ----

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data <- set_up_dataset(differential = TRUE)

# proteomics variables
variables_proteomics <- data %>% select(contains("proteom")) %>% colnames()

# load libraries
library(FactoMineR)
library(factoextra)


# Data formatting ----

data_analysis <- data %>%
  # select only proteomics
  select(all_of(c("study_id", variables_proteomics))) %>%
  # drop all rows with missigness (mostly patients without drug data)
  drop_na(variables_proteomics) %>%
  # scale all values to be normalised
  mutate(across(all_of(variables_proteomics), ~ scale(.)[,1]))


# PCA ----

## Run PCA ----
pca_analysis <- FactoMineR::PCA(
  X = data_analysis %>%
    select(-c("study_id")),
  scale.unit = FALSE, # already done
  ncp = 10 # this does not affect pca, affects the components given 
)


## Explained variance ----
# Variance explained by different components
plot_varexplained <- factoextra::fviz_eig(pca_analysis, addlabels = TRUE, ncp = 10)

# Attribute importance for specific components (in this case 1 and 2)
plot_attributeimp <- factoextra::fviz_pca_var(
  pca_analysis, 
  col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE,
  axes = c(1, 2)
  )

# COS2 contribution to specific components
plot_cos2imp <- factoextra::fviz_cos2(
  pca_analysis,
  choice = "var",
  axes = c(1, 2)
)


# Extract PCA info ----

## Patient values ----
pca_individuals <- factoextra::get_pca_ind(pca_analysis)$coord %>%
  as.data.frame() %>%
  set_names(paste0("PC", 1:10)) %>% # change column names
  cbind(
    study_id = data_analysis$study_id # add study id column
  ) %>%
  relocate("study_id", .before = "PC1") # place it at the beginning of the data.frame

# save table
dir.create("Data")
save(pca_individuals, file = "Data/pca_individuals.Rdata")

## Variable values ----

### Extract PCA values for each variable (used to understand vars in PCA after regression)
pca_variable <- factoextra::get_pca_var(pca_analysis)

# save table
save(pca_variable, file = "Data/pca_variable.Rdata")




# # Fit regression for specific benefit
# ## With PCA
# lm(benefit ~ prehba1c + PCA1 + PCA2 + PCA3+ PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_analysis %>% left_join(data %>% select(study_id, prehba1c)) %>%left_join(pca_individuals))
# 
# ## With PCA and BMI
# lm(benefit ~ prehba1c + prebmi + preegfr + PCA1 + PCA2 + PCA3+ PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_analysis %>% left_join(data %>% select(study_id, prehba1c, prebmi, preegfr)) %>%left_join(pca_individuals))
# 
# ## do LASSO for regression to reduce components to 0 (check if it discards anything)
# library(glmnet)
# x = data_analysis %>% left_join(data %>% select(study_id, prehba1c, prebmi, preegfr)) %>%left_join(pca_individuals) %>%
#   select(-c("benefit", "study_id", variables_proteomics)) %>%
#   as.matrix()
# y = data_analysis %>% left_join(data %>% select(study_id, prehba1c, prebmi, preegfr)) %>%left_join(pca_individuals) %>%
#   select(c("benefit")) %>%
#   unlist()
# cv_model <- cv.glmnet(x, y, alpha = 1)
# best_model <- glmnet(x, y, alpha = 1, lambda = cv_model$lambda.min)
# 
# # Check correlations for each PCA
# pca_variable$cor[,7] %>% View()
# pca_variable$cor[,8] %>% View()
# 
# # Check highest (absolute) correlation between var and PCA
# cor_highest <- pca_variable$cor[order(abs(pca_variable$cor[,8]), decreasing = TRUE)[1:30],8]




