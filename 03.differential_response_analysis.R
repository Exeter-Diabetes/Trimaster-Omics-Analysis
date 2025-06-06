
# Initial setup ----

# load functions
source("00.functions.R")
source("01.cohort_definition.R")

# load data
data <- set_up_dataset(differential = TRUE)

# proteomics variables
variables_proteomics <- data %>% select(contains("proteom")) %>% colnames()


# Differential response analysis ----

## Univariate analysis ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_univariate <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_univariate <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_univariate <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics
)



## Adjusted analysis ----

### HbA1c ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_adj_hba1c <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = "prehba1c"
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_adj_hba1c <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = "prehba1c"
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_adj_hba1c <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = "prehba1c"
)

### Sex (+ HbA1c) ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_adj_hba1c_sex <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("sex", "prehba1c")
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_adj_hba1c_sex <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("sex", "prehba1c")
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_adj_hba1c_sex <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("sex", "prehba1c")
)


### BMI (+ HbA1c) ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_adj_hba1c_bmi <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prebmi", "prehba1c")
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_adj_hba1c_bmi <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prebmi", "prehba1c")
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_adj_hba1c_bmi <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prebmi", "prehba1c")
)


### eGFR (+ HbA1c) ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_adj_hba1c_egfr <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("preegfr", "prehba1c")
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_adj_hba1c_egfr <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("preegfr", "prehba1c")
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_adj_hba1c_egfr <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("preegfr", "prehba1c")
)


### Age at treatment (+ HbA1c) ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_adj_hba1c_agetx <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "agetx")
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_adj_hba1c_agetx <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "agetx")
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_adj_hba1c_agetx <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "agetx")
)


### Duration (+ HbA1c) ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_adj_hba1c_t2dmduration <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "t2dmduration")
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_adj_hba1c_t2dmduration <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "t2dmduration")
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_adj_hba1c_t2dmduration <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "t2dmduration")
)


### HbA1c + BMI + Sex + eGFR + agetx + duration ----

## DPP4 vs SGLT2 (negative is DPP4, positive is SGLT2)
DPP4_SGLT2_pvalues_adj <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_SGLT2),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
)

## DPP4 vs TZD (negative is DPP4, positive is TZD)
DPP4_TZD_pvalues_adj <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_DPP4 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
)

## SGLT2 vs TZD (negative is SGLT2, positive is TZD)
SGLT2_TZD_pvalues_adj <- regression_analysis_function(
  data = data %>%
    mutate(benefit = resphba1c_SGLT2 - resphba1c_TZD),
  var_outcome = "benefit", var_proteomics = variables_proteomics, adj_vars = c("prehba1c", "prebmi", "sex", "agetx", "preegfr", "t2dmduration")
)



# Plots ----

## DPP4 vs SGLT2 ----
plot_DPP4_SGLT2_pvalues <- DPP4_SGLT2_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    DPP4_SGLT2_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    DPP4_SGLT2_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    DPP4_SGLT2_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    DPP4_SGLT2_pvalues_adj_hba1c_egfr %>%
      mutate(facet_title = "Adjusted HbA1c + eGFR"),
    DPP4_SGLT2_pvalues_adj_hba1c_agetx %>%
      mutate(facet_title = "Adjusted HbA1c + Age at treatment"),
    DPP4_SGLT2_pvalues_adj_hba1c_t2dmduration %>%
      mutate(facet_title = "Adjusted HbA1c + Duration"),
    DPP4_SGLT2_pvalues_adj %>%
      mutate(facet_title = "Adjusted")
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
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted HbA1c + Sex", "Adjusted HbA1c + BMI", "Adjusted HbA1c + eGFR", "Adjusted HbA1c + Age at treatment", "Adjusted HbA1c + Duration", "Adjusted"))
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


plot_DPP4_SGLT2_coef <- DPP4_SGLT2_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    DPP4_SGLT2_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    DPP4_SGLT2_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    DPP4_SGLT2_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    DPP4_SGLT2_pvalues_adj_hba1c_egfr %>%
      mutate(facet_title = "Adjusted HbA1c + eGFR"),
    DPP4_SGLT2_pvalues_adj_hba1c_agetx %>%
      mutate(facet_title = "Adjusted HbA1c + Age at treatment"),
    DPP4_SGLT2_pvalues_adj_hba1c_t2dmduration %>%
      mutate(facet_title = "Adjusted HbA1c + Duration"),
    DPP4_SGLT2_pvalues_adj %>%
      mutate(facet_title = "Adjusted")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    xaxis = "Coefficient"
  ) %>%
  drop_na(proteomic_yaxis) %>%
  mutate(
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted HbA1c + Sex", "Adjusted HbA1c + BMI", "Adjusted HbA1c + eGFR", "Adjusted HbA1c + Age at treatment", "Adjusted HbA1c + Duration", "Adjusted"))
  ) %>%
  ggplot(aes(x = xaxis, y = proteomic_yaxis, fill = coef)) +
  geom_tile() +
  labs(y = "Proteomics", fill = "Direction of effect", subtitle = "Negative = DPP4, Positive = SGLT2") +
  theme_classic() +
  scale_fill_gradient2(low = "green", high = "red", midpoint = 0) +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )



## DPP4 vs TZD ----
plot_DPP4_TZD_pvalues <- DPP4_TZD_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    DPP4_TZD_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    DPP4_TZD_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    DPP4_TZD_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    DPP4_TZD_pvalues_adj_hba1c_egfr %>%
      mutate(facet_title = "Adjusted HbA1c + eGFR"),
    DPP4_TZD_pvalues_adj_hba1c_agetx %>%
      mutate(facet_title = "Adjusted HbA1c + Age at treatment"),
    DPP4_TZD_pvalues_adj_hba1c_t2dmduration %>%
      mutate(facet_title = "Adjusted HbA1c + Duration"),
    DPP4_TZD_pvalues_adj %>%
      mutate(facet_title = "Adjusted")
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
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted HbA1c + Sex", "Adjusted HbA1c + BMI", "Adjusted HbA1c + eGFR", "Adjusted HbA1c + Age at treatment", "Adjusted HbA1c + Duration", "Adjusted"))
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


plot_DPP4_TZD_coef <- DPP4_TZD_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    DPP4_TZD_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    DPP4_TZD_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    DPP4_TZD_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    DPP4_TZD_pvalues_adj_hba1c_egfr %>%
      mutate(facet_title = "Adjusted HbA1c + eGFR"),
    DPP4_TZD_pvalues_adj_hba1c_agetx %>%
      mutate(facet_title = "Adjusted HbA1c + Age at treatment"),
    DPP4_TZD_pvalues_adj_hba1c_t2dmduration %>%
      mutate(facet_title = "Adjusted HbA1c + Duration"),
    DPP4_TZD_pvalues_adj %>%
      mutate(facet_title = "Adjusted")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    xaxis = "Coefficient"
  ) %>%
  drop_na(proteomic_yaxis) %>%
  mutate(
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted HbA1c + Sex", "Adjusted HbA1c + BMI", "Adjusted HbA1c + eGFR", "Adjusted HbA1c + Age at treatment", "Adjusted HbA1c + Duration", "Adjusted"))
  ) %>%
  ggplot(aes(x = xaxis, y = proteomic_yaxis, fill = coef)) +
  geom_tile() +
  labs(y = "Proteomics", fill = "Direction of effect", subtitle = "Negative = DPP4, Positive = TZD") +
  theme_classic() +
  scale_fill_gradient2(low = "green", high = "red", midpoint = 0) +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )



## SGLT2 vs TZD ----
plot_SGLT2_TZD_pvalues <- SGLT2_TZD_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    SGLT2_TZD_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    SGLT2_TZD_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    SGLT2_TZD_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    SGLT2_TZD_pvalues_adj_hba1c_egfr %>%
      mutate(facet_title = "Adjusted HbA1c + eGFR"),
    SGLT2_TZD_pvalues_adj_hba1c_agetx %>%
      mutate(facet_title = "Adjusted HbA1c + Age at treatment"),
    SGLT2_TZD_pvalues_adj_hba1c_t2dmduration %>%
      mutate(facet_title = "Adjusted HbA1c + Duration"),
    SGLT2_TZD_pvalues_adj %>%
      mutate(facet_title = "Adjusted")
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
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted HbA1c + Sex", "Adjusted HbA1c + BMI", "Adjusted HbA1c + eGFR", "Adjusted HbA1c + Age at treatment", "Adjusted HbA1c + Duration", "Adjusted"))
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


plot_SGLT2_TZD_coef <- SGLT2_TZD_pvalues_univariate %>%
  mutate(facet_title = "Unadjusted") %>%
  rbind(
    SGLT2_TZD_pvalues_adj_hba1c %>%
      mutate(facet_title = "Adjusted HbA1c"),
    SGLT2_TZD_pvalues_adj_hba1c_sex %>%
      mutate(facet_title = "Adjusted HbA1c + Sex"),
    SGLT2_TZD_pvalues_adj_hba1c_bmi %>%
      mutate(facet_title = "Adjusted HbA1c + BMI"),
    SGLT2_TZD_pvalues_adj_hba1c_egfr %>%
      mutate(facet_title = "Adjusted HbA1c + eGFR"),
    SGLT2_TZD_pvalues_adj_hba1c_agetx %>%
      mutate(facet_title = "Adjusted HbA1c + Age at treatment"),
    SGLT2_TZD_pvalues_adj_hba1c_t2dmduration %>%
      mutate(facet_title = "Adjusted HbA1c + Duration"),
    SGLT2_TZD_pvalues_adj %>%
      mutate(facet_title = "Adjusted")
  ) %>%
  mutate(
    proteomic_yaxis = ifelse(p_value <= 0.05, gsub("proteomics_", "", proteomic), NA),
    xaxis = "Coefficient"
  ) %>%
  drop_na(proteomic_yaxis) %>%
  mutate(
    facet_title = factor(facet_title, levels = c("Unadjusted", "Adjusted HbA1c", "Adjusted HbA1c + Sex", "Adjusted HbA1c + BMI", "Adjusted HbA1c + eGFR", "Adjusted HbA1c + Age at treatment", "Adjusted HbA1c + Duration", "Adjusted"))
  ) %>%
  ggplot(aes(x = xaxis, y = proteomic_yaxis, fill = coef)) +
  geom_tile() +
  labs(y = "Proteomics", fill = "Direction of effect", subtitle = "Negative = SGLT2, Positive = TZD") +
  theme_classic() +
  scale_fill_gradient2(low = "green", high = "red", midpoint = 0) +
  guides(fill = guide_colorbar(barwidth = unit(10, "cm"))) +
  facet_wrap(~facet_title, nrow = 2) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )



# PDFs ----
pdf("Plots/03.differential_response_analysis.pdf", width = 12, height = 8)
plot_DPP4_SGLT2_pvalues
plot_DPP4_SGLT2_coef
plot_DPP4_TZD_pvalues
plot_DPP4_TZD_coef
plot_SGLT2_TZD_pvalues
plot_SGLT2_TZD_coef
dev.off()



# # something
# SGLT2_TZD_pvalues_univariate %>% slice(1:3)
# SGLT2_TZD_pvalues_adj_hba1c %>% slice(1:3)
# # FABP4: Fatty acid-binding protein, adipocyte
# # LEP: Leptin (protein hormone produced by fat cells)







