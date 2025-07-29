
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


