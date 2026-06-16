# # # # # # # # # # # # # # # # # # #
# In this file, we do a univariate analysis of proteomics versus drug response
# # # # # # # # # # # # # # # # # # # 

regression_analysis_function <- function(data, var_outcome = NULL, var_proteomics = NULL, drug, adj_vars = NULL) {
  
  # check if var_outcome isn't null
  if (is.null(var_outcome)) stop("Supply var_outcome.")
  if (!all(var_outcome %in% colnames(data))) stop("var_outcome not in the data")
  # check if var_proteomics isn't null
  if (is.null(var_proteomics)) stop("Supply var_proteomics.")
  if (!all(var_proteomics %in% colnames(data))) stop("Not all var_proteomics in data")
  # check whether adj_vars are in the data
  if (!is.null(adj_vars)) if(!all(adj_vars %in% colnames(data))) stop("Not all adj_vars in data")
  
  # output file
  interim_dataset <- NULL
  
  # iterate by each var_proteomics
  for (protein in var_proteomics) {
    
    # formula for iterations (adj_vars only added when not null because of paste)
    formula <- paste(paste(var_outcome, "~"), paste(c(protein, adj_vars), collapse = "+"))
    
    # model for iterations
    model <- lm(formula = as.formula(formula), data = data)
    
    # get confidence intervals
    ci_vals <- confint(model, parm = protein, level = 0.95)
    
    interim_dataset <- rbind(
      interim_dataset,
      data.frame(
        proteomic = protein,
        coef = summary(model)$coefficients[protein, 1],
        ci_low = ci_vals[1],
        ci_high = ci_vals[2],
        p_value = summary(model)$coefficients[protein, 4]
      )
    )
    
  }
  
  # multiple testing adjustment
  interim_dataset <- interim_dataset %>%
    mutate(
      p_value_adj_bonf = p.adjust(interim_dataset$p_value, method = "bonferroni"),
      p_value_adj_FDR = p.adjust(interim_dataset$p_value, method = "BH")
    ) %>%
    arrange(p_value)
  

  return(interim_dataset)
  
}

  

  
  
run_top_protein_regressions <- function(data_drug, outcome_var, clinical_vars, output_proteomics, drug_name) {
    
    # PCs that were significantly associated with differential drug response
    pcs_to_process <- unique(output_proteomics$PCA)
    results_list <- list()
    
    for (pc in pcs_to_process) {
      cat("Processing", pc, "for", drug_name, "\n")
      
      # Get top 10 proteins for that PC
      top_proteins <- output_proteomics %>%
        filter(PCA == pc) %>%
        pull(Names) %>%
        tolower() %>%
        paste0("proteomics_", .)
      
      # Check if any proteins are missing from the dataset
      missing_proteins <- top_proteins[!top_proteins %in% colnames(data_drug)]
      if (length(missing_proteins) > 0) {
        cat("Missing proteins for", pc, ":", paste(missing_proteins, collapse = ", "), "\n")
      }
      

      # Run regression (if any valid proteins)
      if (length(top_proteins) > 0) {
        regression_results <- regression_analysis_function(
          data = data_drug,
          var_outcome = outcome_var,
          var_proteomics = top_proteins,
          adj_vars = clinical_vars,
          drug = drug_name
        )
        
        regression_results$drug <- drug_name
        regression_results$PCA <- pc
        results_list[[pc]] <- regression_results
      }
    }
    
    # Combine results for each PC
    results_df <- bind_rows(results_list)
    return(results_df)
}



plot_forest <- function(model, 
                        title = title, 
                        save_path = NULL, 
                        width = 8, 
                        height = 6) {
  
  # Tidy model results
  model_results <- broom::tidy(model, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    # Keep only PCs
    filter(grepl("^PC", term)) %>%
    mutate(
      term = factor(term, levels = rev(term)),
      significance = ifelse(p.value < 0.05, "p < 0.05", "Non significant")
    )
  
  # Create plot
  p <- ggplot(model_results, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = significance)) +
    geom_pointrange(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    coord_flip() +
    scale_color_manual(values = c("p < 0.05" = "firebrick", "Non significant" = "black")) +
    labs(
      title = title,
      x = "",
      y = "Effect Estimate (95% CI)",
      color = "Significance"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
  
  # Save plot if path is provided
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = 300)
  }
  
  return(p)
}



# Function to generate summary tables for proteins in PCs significantly associated with differential drug response showing:
# 1) How each protein contributes to the differential drug response and in which direction 
# 2) How each protein is associated with response to each individual drug
generate_pc_summary_tables <- function(output_proteomics,
                                       results_drug1,
                                       results_drug2,
                                       drug1_label = "Drug1",
                                       drug2_label = "Drug2") {
  format_protein <- function(x) {
    paste0("proteomics_", tolower(x))
  }
  
  pcs <- unique(output_proteomics$PCA)
  summary_tables <- list()
  
  for (pc in pcs) {
    
    # Subset output for this PC (keeps coef, ci_low, ci_high from PC*_adj)
    pc_df <- output_proteomics %>%
      dplyr::filter(PCA == pc) %>%
      dplyr::mutate(protein_join = format_protein(Names))  # for joining
    
    # Response results for drug 1 (keep coef + CI)
    drug1_pc <- results_drug1 %>%
      dplyr::filter(PCA == pc) %>%
      dplyr::select(proteomic, coef, ci_low, ci_high) %>%
      dplyr::rename(!!paste0(drug1_label, "_Response") := coef,
                    !!paste0(drug1_label, "_CI_low")  := ci_low,
                    !!paste0(drug1_label, "_CI_high") := ci_high)
    
    # Response results for drug 2 (keep coef + CI)
    drug2_pc <- results_drug2 %>%
      dplyr::filter(PCA == pc) %>%
      dplyr::select(proteomic, coef, ci_low, ci_high) %>%
      dplyr::rename(!!paste0(drug2_label, "_Response") := coef,
                    !!paste0(drug2_label, "_CI_low")  := ci_low,
                    !!paste0(drug2_label, "_CI_high") := ci_high)
    
    # Join differential and individual drug response
    summary_table <- pc_df %>%
      dplyr::left_join(drug1_pc,  by = c("protein_join" = "proteomic")) %>%
      dplyr::left_join(drug2_pc,  by = c("protein_join" = "proteomic")) %>%
      dplyr::select(
        `Gene description`,
        benefit,
        coef, ci_low, ci_high,                           # differential (PC-level) coef + CI
        !!rlang::sym(paste0(drug1_label, "_Response")),
        !!rlang::sym(paste0(drug1_label, "_CI_low")),
        !!rlang::sym(paste0(drug1_label, "_CI_high")),
        !!rlang::sym(paste0(drug2_label, "_Response")),
        !!rlang::sym(paste0(drug2_label, "_CI_low")),
        !!rlang::sym(paste0(drug2_label, "_CI_high"))
      ) %>%
      dplyr::rename(
        `Drug Benefit Direction`      = benefit,
        `Differential Drug Response`  = coef,
        `Differential_CI_low`         = ci_low,
        `Differential_CI_high`        = ci_high
      ) %>%
      dplyr::mutate(
        `Differential Drug Response` = round(`Differential Drug Response`, 2),
        dplyr::across(dplyr::matches("CI_|_Response$"), ~ round(.x, 2))
      )
    
    summary_tables[[pc]] <- summary_table
  }
  
  return(summary_tables)
}


  
