# # # # # # # # # # # # # # # # # # #
# In this file, we do a univariate analysis of proteomics versus drug response
# # # # # # # # # # # # # # # # # # # 

#Functions ----
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
    
    # append new rows
    interim_dataset <- rbind(
      interim_dataset,
      data.frame(
        proteomic = protein,
        coef = summary(model)$coefficients[protein,1],
        p_value = summary(model)$coefficients[protein,4]
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
  
  # output object
  return(interim_dataset)
  
}

