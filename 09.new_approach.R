# this approach is:
# 1. check which vars are associated (p<0.05) with differential effect in univariate
# 2. take the union of those selected for all comparisons
# 3. PCA on the remaining proteomics
# 4. see which PCs are associated post clinical features

# Initial setup ----

# load libraries
library(FactoMineR)
library(factoextra)
library(glmnet)

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data <- set_up_dataset(differential = TRUE)
load("Data/pca_individuals.Rdata") # pca_individuals
load("Data/pca_variable.Rdata") # pca_variable

variables_proteomics <- data %>% select(contains("proteom")) %>% colnames()
clinical_features <- c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
pc_features <- paste0("PC", 1:10)

# Proteomics summary
human_protein_atlas <- read.delim("Proteomics Descriptions/human_protein_atlas.tsv")
unitprotkb_function <- readxl::read_excel("Proteomics Descriptions/uniprotkb_function.xlsx")

## Data formatting ----
data_analysis_sglt2_tzd <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, variables_proteomics, "resphba1c_SGLT2", "resphba1c_TZD"))) %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  select(-c("resphba1c_SGLT2", "resphba1c_TZD")) %>%
  # drop missing data
  drop_na() %>%
  # scale features (but not study_id or benefit)
  mutate(
    across(setdiff(c(clinical_features, variables_proteomics), "sex"), 
           ~ as.numeric(scale(.x)))
  )



data_analysis_sglt2_dpp4 <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, variables_proteomics, "resphba1c_SGLT2", "resphba1c_DPP4"))) %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_DPP4) %>%
  select(-c("resphba1c_SGLT2", "resphba1c_DPP4")) %>%
  # drop missing data
  drop_na() %>%
  # scale features (but not study_id or benefit)
  mutate(
    across(setdiff(c(clinical_features, variables_proteomics), "sex"), 
           ~ as.numeric(scale(.x)))
  )



data_analysis_tzd_dpp4 <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, variables_proteomics, "resphba1c_TZD", "resphba1c_DPP4"))) %>%
  mutate(benefit = resphba1c_TZD - resphba1c_DPP4) %>%
  select(-c("resphba1c_TZD", "resphba1c_DPP4")) %>%
  # drop missing data
  drop_na() %>%
  # scale features (but not study_id or benefit)
  mutate(
    across(setdiff(c(clinical_features, variables_proteomics), "sex"), 
           ~ as.numeric(scale(.x)))
  )


# Unadjusted Regression analysis ----
sglt2_tzd_unadjusted <- regression_analysis_function(
  data = data_analysis_sglt2_tzd,
  var_outcome = "benefit", 
  var_proteomics = variables_proteomics
) %>%
  filter(p_value <= 0.05)

sglt2_dpp4_unadjusted <- regression_analysis_function(
  data = data_analysis_sglt2_dpp4,
  var_outcome = "benefit", 
  var_proteomics = variables_proteomics
) %>%
  filter(p_value <= 0.05)

tzd_dpp4_unadjusted <- regression_analysis_function(
  data = data_analysis_tzd_dpp4,
  var_outcome = "benefit", 
  var_proteomics = variables_proteomics
) %>%
  filter(p_value <= 0.05)


variables_proteomics_simpler <- unique(c(sglt2_tzd_unadjusted$proteomic, sglt2_dpp4_unadjusted$proteomic, tzd_dpp4_unadjusted$proteomic))


# PCA ----

## Run PCA ----
n_pca <- 20
pca_analysis <- FactoMineR::PCA(
  X = data %>%
    select(all_of(variables_proteomics_simpler)) %>%
    drop_na(),
  scale.unit = FALSE, # already done
  ncp = n_pca # this does not affect pca, affects the components given 
)


## Explained variance ----
# Variance explained by different components
plot_varexplained <- factoextra::fviz_eig(pca_analysis, 
                                          addlabels = TRUE, 
                                          ncp = 10) +
  theme_minimal(base_size = n_pca)


## Extract PCA values
pca_individuals <- factoextra::get_pca_ind(pca_analysis)$coord %>%
  as.data.frame() %>%
  set_names(paste0("PC", 1:n_pca)) %>% # change column names
  cbind(
    study_id = data %>% select(all_of(c("study_id", variables_proteomics_simpler))) %>% drop_na() %>% select(study_id) %>% unlist() # add study id column
  ) %>%
  relocate("study_id", .before = "PC1") # place it at the beginning of the data.frame

### Extract PCA values for each variable (used to understand vars in PCA after regression)
pca_variable <- factoextra::get_pca_var(pca_analysis)


# PC analysis adjusted by clinical features ----

## Data formatting ----
data_analysis_sglt2_tzd <- data_analysis_sglt2_tzd %>%
  left_join(
    pca_individuals, by = "study_id"
  )

data_analysis_sglt2_dpp4 <- data_analysis_sglt2_dpp4 %>%
  left_join(
    pca_individuals, by = "study_id"
  )

data_analysis_tzd_dpp4 <- data_analysis_tzd_dpp4 %>%
  left_join(
    pca_individuals, by = "study_id"
  )

## Regression
formula <- paste("benefit ~", paste(clinical_features, collapse = "+"), "+", paste(pc_features, collapse = "+"))

lm_sglt2_tzd <- lm(
  formula = formula,
  data = data_analysis_sglt2_tzd
)

lm_sglt2_dpp4 <- lm(
  formula = formula,
  data = data_analysis_sglt2_dpp4
)

lm_tzd_dpp4 <- lm(
  formula = formula,
  data = data_analysis_tzd_dpp4
)


pca_variable$cor[order(abs(pca_variable$cor[,1]), decreasing = TRUE)[1:10],1]
pca_variable$cor[order(abs(pca_variable$cor[,2]), decreasing = TRUE)[1:10],2]
pca_variable$cor[order(abs(pca_variable$cor[,8]), decreasing = TRUE)[1:10],8]


