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
clinical_features <- c("prehba1c")
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
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD,
        prehba1c = scale(prehba1c)) %>%
  select(-c("resphba1c_SGLT2", "resphba1c_TZD")) %>%
  # join PCA to main dataset
  left_join(
    pca_individuals, by = c("study_id"),
  ) %>%
  # drop missing data
  drop_na()




data_proteomics <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, variables_proteomics, "resphba1c_SGLT2", "resphba1c_TZD"))) %>%
  mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD, 
         prehba1c = scale(prehba1c)) %>%
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

summary(SGLT2_TZD_lm)

# Plot results
plot_forest(SGLT2_TZD_lm, 
            title = "SGLT2 vs TZD (Differential Response)", 
            save_path = "Plots/SGLT2_v_TZD_differential_response_unadjusted.png",
            width = 10, height = 5.5)


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
  adj_vars = c("prehba1c")
)


PC7_adj <- regression_analysis_function(
  data = data_proteomics,
  var_outcome = "benefit", 
  var_proteomics = names(pc_identified[["PC7"]]), 
  adj_vars = c("prehba1c")
)

PC8_adj <- regression_analysis_function(
  data = data_proteomics,
  var_outcome = "benefit", 
  var_proteomics = names(pc_identified[["PC8"]]), 
  adj_vars = c("prehba1c")
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




# Individual drug response ---------


## Data formatting ----
data_SGLT2 <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, variables_proteomics, "resphba1c_SGLT2"))) %>%
  # join PCA to main dataset
  left_join(
    pca_individuals, by = c("study_id")
  ) %>%
  # drop missing data
  drop_na() %>%
  mutate(
    resp = resphba1c_SGLT2,
    prehba1c = scale(prehba1c)
  )


# Get associations of individual proteins with drug resposne
results_protein_SGLT2 <- run_top_protein_regressions(
  data_drug = data_SGLT2,
  outcome_var = "resp",
  clinical_vars = clinical_features,
  output_proteomics = output_proteomics,
  drug_name = "SGLT2"
)


# --- TZD response ----
data_TZD <- data %>%
  # select only clinical features and outcomes
  select(all_of(c("study_id", clinical_features, variables_proteomics, "resphba1c_TZD"))) %>%
  # join PCA to main dataset
  left_join(
    pca_individuals, by = c("study_id")
  ) %>%
  # drop missing data
  drop_na() %>%
  mutate(
    resp = resphba1c_TZD, 
    prehba1c = scale(prehba1c)
  )


# Get associations of individual proteins with drug response
results_protein_TZD <- run_top_protein_regressions(
  data_drug = data_TZD,
  outcome_var = "resp",
  clinical_vars = clinical_features,
  output_proteomics = output_proteomics,
  drug_name = "TZD"
)


# Summarise response data
summary_sglt2_tzd <- generate_pc_summary_tables(
  output_proteomics = output_proteomics %>% filter(p_value < 0.05),
  results_drug1 = results_protein_SGLT2,
  results_drug2 = results_protein_TZD,
  drug1_label = "SGLT2",
  drug2_label = "TZD"
)


# Save summary tables
for (pc in names(summary_sglt2_tzd)) {
  write.csv(summary_sglt2_tzd[[pc]],
            file = paste0("Output/summary_table_SGLT2_v_TZD_unadjusted_", pc, ".csv"),
            row.names = FALSE)
}




