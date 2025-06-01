
# Initial setup ----

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data <- set_up_dataset()

# proteomics variables
variables_proteomics <- data %>% select(contains("proteom")) %>% colnames()


# Treatment response analysis ----

## Univariate analysis ----

## DPP4
DPP4_pvalues_univariate <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics
)

## SGLT2
SGLT2_pvalues_univariate <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics
)

## TZD
TZD_pvalues_univariate <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics
)

## Adjusted analysis ----

### HbA1c ----

## DPP4
DPP4_pvalues_adj_hba1c <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "prehba1c"
)

## SGLT2
SGLT2_pvalues_adj_hba1c <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "prehba1c"
)

## TZD
TZD_pvalues_adj_hba1c <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "prehba1c"
)


### Sex ----

## DPP4
DPP4_pvalues_adj_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "sex"
)

## SGLT2
SGLT2_pvalues_adj_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "sex"
)

## TZD
TZD_pvalues_adj_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "sex"
)


### BMI ----

## DPP4
DPP4_pvalues_adj_bmi <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "prebmi"
)

## SGLT2
SGLT2_pvalues_adj_bmi <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "prebmi"
)

## TZD
TZD_pvalues_adj_bmi <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = "prebmi"
)


### HbA1c + BMI ----

## DPP4
DPP4_pvalues_adj_hba1c_bmi <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi")
)

## SGLT2
SGLT2_pvalues_adj_hba1c_bmi <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi")
)

## TZD
TZD_pvalues_adj_hba1c_bmi <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi")
)



### HbA1c + Sex ----

## DPP4
DPP4_pvalues_adj_hba1c_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "sex")
)

## SGLT2
SGLT2_pvalues_adj_hba1c_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "sex")
)

## TZD
TZD_pvalues_adj_hba1c_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "sex")
)


### BMI + Sex ----

## DPP4
DPP4_pvalues_adj_bmi_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prebmi", "sex")
)

## SGLT2
SGLT2_pvalues_adj_bmi_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prebmi", "sex")
)

## TZD
TZD_pvalues_adj_bmi_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prebmi", "sex")
)



### HbA1c + BMI + Sex ----

## DPP4
DPP4_pvalues_adj_hba1c_bmi_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "DPP4"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi", "sex")
)

## SGLT2
SGLT2_pvalues_adj_hba1c_bmi_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "SGLT2"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi", "sex")
)

## TZD
TZD_pvalues_adj_hba1c_bmi_sex <- regression_analysis_function(
  data = data %>% filter(drugclass == "TZD"), var_outcome = "resphba1c", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi", "sex")
)





# Plots ----

## DPP4 ----
plot_DPP4_pvalues <- DPP4_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    DPP4_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    DPP4_pvalues_adj_sex %>%
      mutate(facet_title = "Adjusted Sex"),
    DPP4_pvalues_adj_bmi %>%
      mutate(facet_title = "Adjusted BMI"),
    DPP4_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    DPP4_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    DPP4_pvalues_adj_bmi_sex %>%
      mutate(facet_title = "Adjusted BMI + Sex"),
    DPP4_pvalues_adj_hba1c_bmi_sex %>%
      mutate(facet_title = "Adjusted HbA1c + BMI + Sex")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    pvalue_fill = ifelse(p_value <= 0.05, p_value, NA),
    pvalue_bonf_fill = ifelse(p_value_adj_bonf <= 0.05, p_value_adj_bonf, NA),
    pvalue_FDR_fill = ifelse(p_value_adj_FDR <= 0.05, p_value_adj_FDR, NA)
  ) %>%
  drop_na(proteomic_yaxis) %>%
  gather("key", "value", -c("proteomic", "coef", "p_value", "p_value_adj_bonf", "p_value_adj_FDR", "proteomic_yaxis", "facet_title")) %>%
  mutate(
    key = factor(key, levels = c("pvalue_fill", "pvalue_FDR_fill", "pvalue_bonf_fill"), labels = c("Unadjusted", "FDR", "Bonferroni")),
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted Sex", "Adjusted BMI", "Adjusted HbA1c + BMI", "Adjusted HbA1c + Sex", "Adjusted BMI + Sex", "Adjusted HbA1c + BMI + Sex"))
  ) %>%
  ggplot(aes(x = key, y = proteomic_yaxis, fill = value)) +
  geom_tile() +
  scale_fill_gradient("p-value", low = "#0072B2", high = "white", na.value = "white", limits = c(0, 0.05), breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), labels = c("0", "0.01", "0.02", "0.03", "0.04", "0.05")) +
  labs(x = "p-value", y = "Proteomics") +
  theme_classic() +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    legend.position = "bottom"
  )


plot_DPP4_coef <- DPP4_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    DPP4_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    DPP4_pvalues_adj_sex %>%
      mutate(facet_title = "Adjusted Sex"),
    DPP4_pvalues_adj_bmi %>%
      mutate(facet_title = "Adjusted BMI"),
    DPP4_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    DPP4_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    DPP4_pvalues_adj_bmi_sex %>%
      mutate(facet_title = "Adjusted BMI + Sex"),
    DPP4_pvalues_adj_hba1c_bmi_sex %>%
      mutate(facet_title = "Adjusted HbA1c + BMI + Sex")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    xaxis = "Coefficient"
  ) %>%
  drop_na(proteomic_yaxis) %>%
  mutate(
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted Sex", "Adjusted BMI", "Adjusted HbA1c + BMI", "Adjusted HbA1c + Sex", "Adjusted BMI + Sex", "Adjusted HbA1c + BMI + Sex"))
  ) %>%
  ggplot(aes(x = xaxis, y = proteomic_yaxis, fill = coef)) +
  geom_tile() +
  labs(y = "Proteomics", fill = "Direction of effect") +
  theme_classic() +
  scale_fill_gradient2(low = "green", high = "red", midpoint = 0) +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )
  


## SGLT2 ----
plot_SGLT2_pvalues <- SGLT2_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    SGLT2_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    SGLT2_pvalues_adj_sex %>%
      mutate(facet_title = "Adjusted Sex"),
    SGLT2_pvalues_adj_bmi %>%
      mutate(facet_title = "Adjusted BMI"),
    SGLT2_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    SGLT2_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    SGLT2_pvalues_adj_bmi_sex %>%
      mutate(facet_title = "Adjusted BMI + Sex"),
    SGLT2_pvalues_adj_hba1c_bmi_sex %>%
      mutate(facet_title = "Adjusted HbA1c + BMI + Sex")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    pvalue_fill = ifelse(p_value <= 0.05, p_value, NA),
    pvalue_bonf_fill = ifelse(p_value_adj_bonf <= 0.05, p_value_adj_bonf, NA),
    pvalue_FDR_fill = ifelse(p_value_adj_FDR <= 0.05, p_value_adj_FDR, NA)
  ) %>%
  drop_na(proteomic_yaxis) %>%
  gather("key", "value", -c("proteomic", "coef", "p_value", "p_value_adj_bonf", "p_value_adj_FDR", "proteomic_yaxis", "facet_title")) %>%
  mutate(
    key = factor(key, levels = c("pvalue_fill", "pvalue_FDR_fill", "pvalue_bonf_fill"), labels = c("Unadjusted", "FDR", "Bonferroni")),
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted Sex", "Adjusted BMI", "Adjusted HbA1c + BMI", "Adjusted HbA1c + Sex", "Adjusted BMI + Sex", "Adjusted HbA1c + BMI + Sex"))
  ) %>%
  ggplot(aes(x = key, y = proteomic_yaxis, fill = value)) +
  geom_tile() +
  scale_fill_gradient("p-value", low = "#E69F00", high = "white", na.value = "white", limits = c(0, 0.05), breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), labels = c("0", "0.01", "0.02", "0.03", "0.04", "0.05")) +
  labs(x = "p-value", y = "Proteomics") +
  theme_classic() +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    legend.position = "bottom"
  )


plot_SGLT2_coef <- SGLT2_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    SGLT2_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    SGLT2_pvalues_adj_sex %>%
      mutate(facet_title = "Adjusted Sex"),
    SGLT2_pvalues_adj_bmi %>%
      mutate(facet_title = "Adjusted BMI"),
    SGLT2_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    SGLT2_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    SGLT2_pvalues_adj_bmi_sex %>%
      mutate(facet_title = "Adjusted BMI + Sex"),
    SGLT2_pvalues_adj_hba1c_bmi_sex %>%
      mutate(facet_title = "Adjusted HbA1c + BMI + Sex")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    xaxis = "Coefficient"
  ) %>%
  drop_na(proteomic_yaxis) %>%
  mutate(
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted Sex", "Adjusted BMI", "Adjusted HbA1c + BMI", "Adjusted HbA1c + Sex", "Adjusted BMI + Sex", "Adjusted HbA1c + BMI + Sex"))
  ) %>%
  ggplot(aes(x = xaxis, y = proteomic_yaxis, fill = coef)) +
  geom_tile() +
  labs(y = "Proteomics", fill = "Direction of effect") +
  theme_classic() +
  scale_fill_gradient2(low = "green", high = "red", midpoint = 0) +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )



## TZD ----
plot_TZD_pvalues <- TZD_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    TZD_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    TZD_pvalues_adj_sex %>%
      mutate(facet_title = "Adjusted Sex"),
    TZD_pvalues_adj_bmi %>%
      mutate(facet_title = "Adjusted BMI"),
    TZD_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    TZD_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    TZD_pvalues_adj_bmi_sex %>%
      mutate(facet_title = "Adjusted BMI + Sex"),
    TZD_pvalues_adj_hba1c_bmi_sex %>%
      mutate(facet_title = "Adjusted HbA1c + BMI + Sex")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    pvalue_fill = ifelse(p_value <= 0.05, p_value, NA),
    pvalue_bonf_fill = ifelse(p_value_adj_bonf <= 0.05, p_value_adj_bonf, NA),
    pvalue_FDR_fill = ifelse(p_value_adj_FDR <= 0.05, p_value_adj_FDR, NA)
  ) %>%
  drop_na(proteomic_yaxis) %>%
  gather("key", "value", -c("proteomic", "coef", "p_value", "p_value_adj_bonf", "p_value_adj_FDR", "proteomic_yaxis", "facet_title")) %>%
  mutate(
    key = factor(key, levels = c("pvalue_fill", "pvalue_FDR_fill", "pvalue_bonf_fill"), labels = c("Unadjusted", "FDR", "Bonferroni")),
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted Sex", "Adjusted BMI", "Adjusted HbA1c + BMI", "Adjusted HbA1c + Sex", "Adjusted BMI + Sex", "Adjusted HbA1c + BMI + Sex"))
  ) %>%
  ggplot(aes(x = key, y = proteomic_yaxis, fill = value)) +
  geom_tile() +
  scale_fill_gradient("p-value", low = "#D55E00", high = "white", na.value = "white", limits = c(0, 0.05), breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), labels = c("0", "0.01", "0.02", "0.03", "0.04", "0.05")) +
  labs(x = "p-value", y = "Proteomics") +
  theme_classic() +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    legend.position = "bottom"
  )


plot_TZD_coef <- TZD_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    TZD_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    TZD_pvalues_adj_sex %>%
      mutate(facet_title = "Adjusted Sex"),
    TZD_pvalues_adj_bmi %>%
      mutate(facet_title = "Adjusted BMI"),
    TZD_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    TZD_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    TZD_pvalues_adj_bmi_sex %>%
      mutate(facet_title = "Adjusted BMI + Sex"),
    TZD_pvalues_adj_hba1c_bmi_sex %>%
      mutate(facet_title = "Adjusted HbA1c + BMI + Sex")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    xaxis = "Coefficient"
  ) %>%
  drop_na(proteomic_yaxis) %>%
  mutate(
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted Sex", "Adjusted BMI", "Adjusted HbA1c + BMI", "Adjusted HbA1c + Sex", "Adjusted BMI + Sex", "Adjusted HbA1c + BMI + Sex"))
  ) %>%
  ggplot(aes(x = xaxis, y = proteomic_yaxis, fill = coef)) +
  geom_tile() +
  labs(y = "Proteomics", fill = "Direction of effect") +
  theme_classic() +
  scale_fill_gradient2(low = "green", high = "red", midpoint = 0) +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )


# PDFs ----
pdf("Plots/02.response_analysis.pdf", width = 12, height = 8)
plot_DPP4_pvalues
plot_DPP4_coef
plot_SGLT2_pvalues
plot_SGLT2_coef
plot_TZD_pvalues
plot_TZD_coef
dev.off()



# # something
# SGLT2_pvalues_adj_sex %>% slice(1:3)
# 
# 
# TZD_pvalues_univariate %>% slice(1:5)
# TZD_pvalues_adj_hba1c %>% slice(1:5)
# TZD_pvalues_adj_sex %>% slice(1:5)
# TZD_pvalues_adj_hba1c_sex %>% slice(1:5)
# TZD_pvalues_adj_bmi %>% slice(1:10)
# # FABP4: Fatty acid-binding protein, adipocyte
# # CES1: Liver carboxylesterase 1
# # LEP: Leptin (protein hormone produced by fat cells)
# # GSTA1: Glutathione S-transferase A1 (glutathione?)
# # AKR1C4: Aldo-keto reductase family 1 member C4 (liver?)



