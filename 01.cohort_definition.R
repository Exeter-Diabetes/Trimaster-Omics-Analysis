# # # # # # # # # # # # # # # # # # # 
# Function for loading initial dataset
# # # # # # # # # # # # # # # # # # # 

set_up_dataset <- function(rem_outliers = FALSE, differential = FALSE) {
  
  # Check inputs ----
  
  # check rem_outliers
  if (!is.logical(rem_outliers)) stop("rem_outliers needs to be TRUE or FALSE")
  # check differential
  if (!is.logical(differential)) stop("differential needs to be TRUE or FALSE")
  
  # load libraries
  require(tidyverse)
  require(readr)
  require(data.table)
  
  
  # Laura's code ----
  # LauraAnalysis\analysis_cluster\study_codes\00_main_analysis_clustervsmodel
  
  
  smoke_data_merge_orig <- read_csv("../../ALLAnalysisReadyData/TM_AnalaysiscombinedData.csv")
  load("../../LauraAnalysis/analysis_cluster/study_data_output/data_wcluster.Rdata") # name: data_wcluster
  
  
  smoke_data_merge      <- data.frame(study_id = as.vector(smoke_data_merge_orig$study_id), smoke = smoke_data_merge_orig$v1_Smoker)
  # table(is.na(smoke_data_merge$smoke))
  
  data <- merge(data_wcluster, smoke_data_merge, by = "study_id", all.x = TRUE)
  
  rm(data_wcluster)
  rm(smoke_data_merge_orig)
  rm(smoke_data_merge)
  #rm(data_wsmoking)
  
  
  #
  ## Load in genetic cluster data and merge
  #
  
  load("../../LauraAnalysis/analysis_cluster/study_data_output//genetic_cluser_dataRdata") # genetic_cluser_data
  
  # table(genetic_cluser_data$study_id%in%data$study_id)
  
  data_merged <- merge(data, genetic_cluser_data, by = "study_id", all.x = T)
  data        <- data_merged
  
  rm(data_merged)
  rm(genetic_cluser_data)
  
  
  # Data preparation
  
  # Exclude participants from the 1_SAID cluster in the analysis
  data_clusterexcluded <- data[!data$clusterID_name == "1_SAID", ]
  
  # generate a few variables needed later
  data_clusterexcluded$drugline <- data_clusterexcluded$vs_SU + data_clusterexcluded$vs_MFN + 1
  # data_clusterexcluded$vs_DIET_only[is.na(data_clusterexcluded$drugline)]
  
  data_clusterexcluded$ncurrtx  <- data_clusterexcluded$drugline
  
  variables_to_include <- c("study_id", "clusterID_name", "validhba1cA", "validhba1cB", "validhba1cC", 
                            "TreatmentOrder", "vs_HbA1c_result", "vs_Age_at_diagnosis", "vs_Gender", 
                            "vs_eGFR", "vs_BMI", "vs_ALT_result", "CHOLmmolLV1", "HDLmmolLV1", "CREATumolLV1",
                            "v1_Blood_Pressure_1sys", "v1_Blood_Pressure_1dia", 
                            "v2v1datediff", "v3v2datediff", "v4v3datediff", "smoke",
                            "vs_Ethnic_Group", "ageatscreening", "drugline", "ncurrtx", "HOMA_B", "HOMA_IR", 
                            "Drug1rank", "Drug1rank2", "Drug2rank", "Drug2rank2", "Drug3rank", "Drug3rank2",
                            "incp1", "incp2", "incp3", 
                            "t2d_betacell_pineg", "t2d_betacell_pipos", "t2d_bodyfat", "t2d_lipodystrophy", 
                            "t2d_liverlipid", "t2d_metabollicsyn", "t2d_obesity", "t2d_residualglycaemic")
  #"tol1", "tol2", "tol3")
  
  # reshape dataset 
  data_long_temp            <- melt(setDT(data_clusterexcluded[, variables_to_include]), 
                                    id.vars = c("study_id", "clusterID_name", "TreatmentOrder", "vs_HbA1c_result",
                                                "vs_Age_at_diagnosis", "vs_Gender", 
                                                "vs_eGFR", "vs_BMI", "vs_ALT_result", "CHOLmmolLV1", "HDLmmolLV1", "CREATumolLV1",
                                                "v1_Blood_Pressure_1sys", "v1_Blood_Pressure_1dia", 
                                                "v2v1datediff", "v3v2datediff", "v4v3datediff",  "smoke",
                                                "vs_Ethnic_Group", "ageatscreening", "drugline", "ncurrtx", 
                                                "HOMA_B", "HOMA_IR", 
                                                "Drug1rank", "Drug1rank2", "Drug2rank", "Drug2rank2", "Drug3rank", "Drug3rank2",
                                                "incp1", "incp2", "incp3", 
                                                "t2d_betacell_pineg", "t2d_betacell_pipos", "t2d_bodyfat", 
                                                "t2d_lipodystrophy", "t2d_liverlipid", "t2d_metabollicsyn", "t2d_obesity", "t2d_residualglycaemic"),
                                    #"tol1", "tol2", "tol3"), 
                                    variable.name = c("hba1c"))
  
  # removed all participants that have not valid hba1c recorded (will not have a study_id in long data format)
  data_long            <- data_long_temp[!is.na(data_long_temp$study_id), ]                  
  
  colnames(data_long)  <- c("study_id", "clusterID_name", "TreatmentOrder", "prehba1c", "Age_at_diagnosis", "sex",
                            "preegfr", "prebmi", "prealt", "pretotalcholesterol", "prehdl", "precreat", "presys", "predia", 
                            "v2v1datediff", "v3v2datediff", "v4v3datediff", "smoke", 
                            "ethnicity", "ageatscreening", "drugline", "ncurrtx", "HOMA_B", "HOMA_IR",
                            "Drug1rank", "Drug1rank2", "Drug2rank", "Drug2rank2", "Drug3rank", "Drug3rank2", 
                            "tol1", "tol2", "tol3", 
                            "t2d_betacell_pineg", "t2d_betacell_pipos", "t2d_bodyfat", 
                            "t2d_lipodystrophy", "t2d_liverlipid", "t2d_metabollicsyn", "t2d_obesity", "t2d_residualglycaemic", 
                            "hba1c_tx", "posthba1cfinal")
  
  # treatment variable
  data_long$tx_taken                                 <- rep(NA, times = dim(data_long)[1])
  data_long$tx_taken[data_long$hba1c_tx == "validhba1cA"] <- "A"
  data_long$tx_taken[data_long$hba1c_tx == "validhba1cB"] <- "B"
  data_long$tx_taken[data_long$hba1c_tx == "validhba1cC"] <- "C"
  
  data_long$drugclass                                      <- rep(NA, times = dim(data_long)[1])
  data_long$drugclass[data_long$hba1c_tx == "validhba1cA"] <- "TZD"
  data_long$drugclass[data_long$hba1c_tx == "validhba1cB"] <- "DPP4"
  data_long$drugclass[data_long$hba1c_tx == "validhba1cC"] <- "SGLT2"
  
  # period variable 
  data_long$period                            <- rep(NA, times = dim(data_long)[1])
  data_long$period[data_long$tx_taken == "A"] <- unlist(gregexpr("A", data_long$TreatmentOrder[data_long$tx_taken == "A"]))
  data_long$period[data_long$tx_taken == "B"] <- unlist(gregexpr("B", data_long$TreatmentOrder[data_long$tx_taken == "B"]))
  data_long$period[data_long$tx_taken == "C"] <- unlist(gregexpr("C", data_long$TreatmentOrder[data_long$tx_taken == "C"]))
  
  # create drugclass as factor with all 5 treatments
  data_long$drugclass                         <- factor(data_long$drugclass, levels = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"))
  
  # create drugline as factor
  data_long$drugline                          <- factor(data_long$drugline, levels = c("2", "3", "4", "5+"))
  
  # create ncurrtx as factor
  data_long$ncurrtx                           <- data_long$ncurrtx - 1            # as per Johns definition
  data_long$ncurrtx                           <- factor(data_long$ncurrtx, levels = c("1", "2", "3", "4+"))
  
  
  # change coding of ethnicity variable to fit with model
  # table(data_long$ethnicity)
  data_long$ethnicity[data_long$ethnicity == "A"] <- "White"
  data_long$ethnicity[data_long$ethnicity == "B"] <- "White"
  data_long$ethnicity[data_long$ethnicity == "C"] <- "White"
  data_long$ethnicity[data_long$ethnicity == "E"] <- "Mixed"
  data_long$ethnicity[data_long$ethnicity == "H"] <- "South Asian"
  data_long$ethnicity[data_long$ethnicity == "J"] <- "South Asian"
  data_long$ethnicity[data_long$ethnicity == "K"] <- "South Asian"
  data_long$ethnicity[data_long$ethnicity == "L"] <- "South Asian"
  data_long$ethnicity[data_long$ethnicity == "N"] <- "Black"
  data_long$ethnicity[data_long$ethnicity == "S"] <- "Other"
  data_long$ethnicity[data_long$ethnicity == "Z"] <- "Missing"
  
  data_long$ethnicity                             <- factor(data_long$ethnicity, levels = c("White", "South Asian", "Black", "Mixed", "Other", "Missing"))
  
  # deprivation index variable
  data_long$imd5                            <- rep("3", times = dim(data_long)[1])
  data_long$imd5                            <- factor(data_long$imd5, levels = c("1 (least)", "2", "3", "4", "5 (most)"))
  
  # age at treatment variable
  data_long$agetx                           <- data_long$ageatscreening  
  
  # sex variable needs to be Female/Male, female is baseline
  data_long$sex[data_long$sex == "M"]        <- "Male"
  data_long$sex[data_long$sex == "F"]        <- "Female"
  
  # str(data_long$sex)
  data_long$sex                              <- as.factor(data_long$sex)
  
  # T2D duration variable
  data_long$t2dmduration                     <- data_long$agetx - data_long$Age_at_diagnosis
  
  # summary(data_long$t2dmduration)
  
  # recode smoke variable
  data_long$smoke[data_long$smoke == 0]   <- "Non-smoker"
  data_long$smoke[data_long$smoke == 1]   <- "Active smoker"
  data_long$smoke[is.na(data_long$smoke)] <- "Not recorded"
  # table(data_long$smoke)
  data_long$smoke                         <- factor(data_long$smoke, levels = c("Non-smoker", "Active smoker", "Ex-smoker", "Not recorded"))
  # str(data_long$smoke)
  # table(data_long$smoke)
  
  # hba1cmonth variable
  data_long$hba1cmonth                    <- rep(NA, times = dim(data_long)[1])
  
  study_id_unique                         <- unique(data_long$study_id)
  for(i in 1: length(study_id_unique)){
    
    subset_participant <- subset(data_long, data_long$study_id == study_id_unique[i])
    
    data_long$hba1cmonth[data_long$study_id%in%study_id_unique[i]] <- as.vector(unlist(subset_participant[1, c("v2v1datediff", "v3v2datediff", "v4v3datediff")]))[subset_participant$period]
    
    # print(i)
  }
  
  data_long$hba1cmonth                     <- data_long$hba1cmonth /30.4
  
  
  # just checking:
  #subset(data_long, data_long$study_id == "TM010002")
  #subset(data_long, data_long$study_id == "TM030300")
  
  # include all periods that have outcome HbA1c recorded rather than all patients that have received all three treatments
  #data_long_temp <- data_long
  #data_long      <- data_long_temp[!is.na(data_long_temp$posthba1cfinal), ]
  
  data_long_cc      <- data_long[!is.na(data_long$posthba1cfinal), ]
  # table(table(data_long_cc$study_id))
  # 50 participants have only one valid hba1c recorded
  # 96 participants have two valid hba1c recorded
  # 309 participants have 3 3 valid hba1c recorded
  
  # I use the data of all completely recorded periods, even from participants who only have one hba1c recorded 
  # to decide on the optimal and worst tx (train the model)
  # but in the CD analysis I will only include the ones that have actually taken their concordant AND discordant treatment 
  
  
  # # Check for missing values
  # contains_any_na = sapply(data_long_cc, function(x) any(is.na(x)))
  # names(data_long_cc)[contains_any_na]

  
  # Above is Laura's code taken from LauraAnalysis\analysis_cluster\study_codes\00_main_analysis_clustervsmodel
  
  # Create response variables ----
  
  clean_dataset <- data_long_cc
  
  # create response variable
  clean_dataset <- clean_dataset %>%
    mutate(resphba1c = posthba1cfinal - prehba1c)
  
  # Omics information ----

  # find names of columns with omics
  require(readxl)
  omics_variables <- read_excel("../../Sample Analysis Data Dictionary.xlsx", sheet = "proteomics")
  colnames(omics_variables) <- c("string", "extra")
  omics_variables <- omics_variables %>%
    select(string) %>%
    slice(-1:-3) %>%
    unlist() %>% 
    str_trim() %>%                    # Remove leading/trailing spaces
    str_split("\\s+") %>%            # Split by one or more spaces
    map_chr(~ .x[1])               # Extract the first word
  
  
  omics_dataset <- read_csv("../../ALLAnalysisReadyData/TM_AnalaysiscombinedData.csv") %>%
    select(contains(c("study_id", omics_variables)))

  # rename variables to have proteomics before hand
  omics_dataset <- omics_dataset %>%
    rename_with(
      .fn = ~ if_else(
        map_lgl(.x, ~ any(str_detect(.x, omics_variables))),
        paste0("proteomics_", .x),
        .x
      )
    )
  
  # variables needed
  omics_dataset <- omics_dataset %>%
    select(contains(c("study_id", "_mean"))) %>%
    rename_with(~ str_remove(.x, "_mean"))
  
  # remove columns with full missingness
  omics_dataset <- omics_dataset %>%
    select(where(~ !all(is.na(.))))
  
  
  if (isTRUE(rem_outliers)) {
      
    ## Potential cleaning step: removing outliers (https://pmc.ncbi.nlm.nih.gov/articles/PMC11405273/)
    
    ## need to stardadise
    omics_dataset <- omics_dataset %>%
      mutate(across(contains("proteomics_"), ~ ifelse(abs(. - mean(., na.rm = TRUE)) > 5 * sd(., na.rm = TRUE), NA, .)))
    
  
  }
  
  # join with main dataset
  dataset_row_drug_initiation <- clean_dataset %>%
    left_join(
      omics_dataset, by = c("study_id")
    )
  
  ### ### ### ### ### ### ###
  
  # Differential response ----
  
  # check if this has been passed (if not, then pass above)
  if (!isTRUE(differential)) {
    
    # return output
    output_dataset <- dataset_row_drug_initiation
    return(output_dataset)
    
  } else {
    
    # pivot dataset wider
    dataset_row_patient <- dataset_row_drug_initiation %>%
      select(-c(hba1c_tx, tx_taken, period)) %>%
      group_by(study_id) %>%
      pivot_wider(
        names_from = c("drugclass"),
        values_from = c("posthba1cfinal", "resphba1c", "hba1cmonth")
      ) %>%
      ungroup()
    
    # return output
    output_dataset <- dataset_row_patient
    return(output_dataset)
    
  }
  
}
