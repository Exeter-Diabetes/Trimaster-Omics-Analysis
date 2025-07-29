
# Initial setup ----

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data <- set_up_dataset(differential = TRUE)
load("Data/pca_individuals.Rdata") # pca_individuals
load("Data/pca_variable.Rdata") # pca_variable

# proteomics variables
variables_proteomics <- data %>% select(contains("proteom")) %>% colnames()
clinical_features <- c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
pc_features <- paste0("PC", 1:10)


# Proteomics summary
human_protein_atlas <- read.delim("Proteomics Descriptions/human_protein_atlas.tsv")
unitprotkb_function <- readxl::read_excel("Proteomics Descriptions/uniprotkb_function.xlsx")

# load libraries
library(FactoMineR)
library(factoextra)
library(glmnet)

## Data formatting ----
data_analysis <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, "resphba1c_SGLT2", "resphba1c_TZD"))) %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  select(-c("resphba1c_SGLT2", "resphba1c_TZD")) %>%
  # join PCA to main dataset
  left_join(
    pca_individuals, by = c("study_id")
  ) %>%
  # drop missing data
  drop_na()

data_proteomics <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, variables_proteomics, "resphba1c_SGLT2", "resphba1c_TZD"))) %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD) %>%
  select(-c("resphba1c_SGLT2", "resphba1c_TZD")) %>%
  # drop missing data
  drop_na()
  

# Regression analysis ----


formula <- paste("benefit ~", paste(clinical_features, collapse = "+"), "+", paste(pc_features, collapse = "+"))

## LASSO ----
x = data_analysis %>% 
  select(all_of(c(clinical_features, pc_features))) %>%
  as.matrix()
y = data_analysis %>%
  select(c("benefit")) %>%
  unlist()
cv_model <- cv.glmnet(x, y, alpha = 1)
best_model <- glmnet(x, y, alpha = 1, lambda = cv_model$lambda.min)


## Linear regression ----
SGLT2_TZD_lm <- lm(
  formula = formula,
  data = data_analysis
)



# Checking PCs identified ----

## PC 4, 7 and 8 identified
pc_identified <- NULL
pc_identified[["PC4"]] <- pca_variable$cor[order(abs(pca_variable$cor[,4]), decreasing = TRUE)[1:10],4]
pc_identified[["PC7"]] <- pca_variable$cor[order(abs(pca_variable$cor[,7]), decreasing = TRUE)[1:10],7]
pc_identified[["PC8"]] <- pca_variable$cor[order(abs(pca_variable$cor[,8]), decreasing = TRUE)[1:10],8]


## Regression analysis ----
# (negative is SGLT2, positive is TZD)
## univariate analysis of each PC proteomics (adjusted for clinical features)
PC4_adj <- regression_analysis_function(
  data = data_proteomics,
  var_outcome = "benefit", 
  var_proteomics = names(pc_identified[["PC4"]]), 
  adj_vars = c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
)

PC7_adj <- regression_analysis_function(
  data = data_proteomics,
  var_outcome = "benefit", 
  var_proteomics = names(pc_identified[["PC7"]]), 
  adj_vars = c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
)

PC8_adj <- regression_analysis_function(
  data = data_proteomics,
  var_outcome = "benefit", 
  var_proteomics = names(pc_identified[["PC8"]]), 
  adj_vars = c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
)


output_proteomics <- rbind(
  PC4_adj %>%
    mutate(
      drug1 = "Towards SGLT2",
      drug2 = "Towards TZD",
      PCA = "PC4"
    ) %>%
    mutate(
      benefit = ifelse(coef <= 0, drug1, drug2),
      proteomic = toupper(gsub("proteomics_", "", proteomic))
    ) %>%
    select(-c(p_value_adj_bonf, p_value_adj_FDR, drug1, drug2)) %>%
    relocate(benefit, .after = proteomic) %>%
    rename(Names = proteomic) %>%
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
    select(-c(`Uniprot`)) %>%
    relocate(PCA, .before = "Names"),
  PC7_adj %>%
    mutate(
      drug1 = "Towards SGLT2",
      drug2 = "Towards TZD",
      PCA = "PC7"
    ) %>%
    mutate(
      benefit = ifelse(coef <= 0, drug1, drug2),
      proteomic = toupper(gsub("proteomics_", "", proteomic))
    ) %>%
    select(-c(p_value_adj_bonf, p_value_adj_FDR, drug1, drug2)) %>%
    relocate(benefit, .after = proteomic) %>%
    rename(Names = proteomic) %>%
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
    select(-c(`Uniprot`)) %>%
    relocate(PCA, .before = "Names"),
  PC8_adj %>%
    mutate(
      drug1 = "Towards SGLT2",
      drug2 = "Towards TZD",
      PCA = "PC8"
    ) %>%
    mutate(
      benefit = ifelse(coef <= 0, drug1, drug2),
      proteomic = toupper(gsub("proteomics_", "", proteomic))
    ) %>%
    select(-c(p_value_adj_bonf, p_value_adj_FDR, drug1, drug2)) %>%
    relocate(benefit, .after = proteomic) %>%
    rename(Names = proteomic) %>%
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
    select(-c(`Uniprot`)) %>%
    relocate(PCA, .before = "Names")
)




