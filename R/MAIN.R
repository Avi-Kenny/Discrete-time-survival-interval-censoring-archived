# Title: "Modeling HIV seroconversion dates"
# Author: Avi Kenny
# Date: 2020-09-30



#################.
##### Setup #####
#################.

# Set working directory
# If running this locally in RStudio, skip this section and click "Session" >>
#     "Set working directory" >> "To Source File Location"
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Jim Hughes/Project - Stepped wedge/z.stepped.wedge/R")
} else {
  setwd("z.stepped.wedge/R")
}

# Set code blocks to run
{
  run_main <- TRUE
  run_results <- FALSE
}



###########################.
##### Simulation code #####
###########################.

if (run_main) {
  
  library(simba)
  library(magrittr)
  
  run_on_cluster(
    
    first = {
      
      # Set up and configure simulation object
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = 10, # !!!!!
        parallel = "cluster",
        packages = c("dplyr", "mice", "survival")
      )
      
      # Specify seroconversion conditional probabilities (discrete hazards)
      # We assume an individual "of age X" has lived exactly X years
      # For an individual of age X, the number in the corresponding age bin
      #     represents the probability that that individual seroconverted
      #     between age X-1 and age X, given that they did not seroconvert by
      #     age X-1
      # !!!!! Change the actual numbers later based on AHRI cohort data
      p_sero_year <- list(
        male = list(
          "1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,
          "26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005
        ),
        female = list(
          "1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,
          "26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0
        )
      )
      
      sim %<>% add_constant(
        m_probs = construct_probs(p_sero_year),
        m = 5, # Number of MI replicates
        num_patients = 50, # !!!!!
        start_year = 2000,
        end_year = 2020
      )
      
      # Add functions to simulation object
      source("generate_dataset.R")
      source("one_simulation.R")
      sim %<>% add_creator(generate_dataset)
      sim %<>% add_script(one_simulation)
      
      # Set levels
      sim %<>% set_levels(
        hr_hiv = 1.8,
        hr_art = 1.2
      )
      
    },
    
    main = { sim %<>% run("one_simulation") },
    
    last = {
      
      # Calculate bias, variance inflation, power
      
      sim %>% summary()
      
    },
    
    cluster_config = list(
      sim_var = "sim",
      js = "slurm",
      dir = "/home/akenny/HIV.multiple.imputation"
    )
    
  )
  
}



###########################.
##### Process results #####
###########################.

if (run_results) {
  
  library(ggplot2)
  
  # !!!!! TO DO
  
}



###################.
##### Archive #####
###################.

if (FALSE) {
  
  # Misc MICE code
  {
    # Run analysis on each imputed dataset
    for (j in 1:m) {
      d_imputed <- mice::complete(imputation_object, j)
    }
    
    # Get degrees of freedom from a single analysis
    dfcom <- summary(cox_model_ideal)$waldtest[[2]]
    
    # Store list as mice::mira object and pool results
    imputation_results <- as.mira(analysis_list)
    pooled_model <- pool(imputation_results, dfcom = dfcom)
    
    # Get coefficient of hiv_status(3)
    coeff_ideal <- summary(cox_model_ideal)$coefficients[2,1]
    coeff_mi <- pooled_model$pooled[2,1]
  }
  
  # Convert yearly probabilities to monthly probabilities
  {
    convert_to_monthly_prob <- function(p) { 1 - (1-p)^(1/12) }
    psero <- list(
      mtct = psero_year$mtct,
      male = lapply(psero_year$male, convert_to_monthly_prob),
      female = lapply(psero_year$female, convert_to_monthly_prob)
    )
  }
  
}
