# Set configuration
source("R/config.R", local=T)

# # !!!!! TEMP DEBUGGING
# cfg$parallelize <- F
# cfg$model_sex <- "Female"
# cfg$samp_size <- 0 # !!!!!

# Load SimEngine + functions
{
  library(SimEngine)
  source("R/helpers.R", local=T)
  source("R/one_simulation.R", local=T)
  source("R/generate_data.R", local=T)
  source("R/likelihood.R", local=T)
}

source("R/models.R", local=T)

if (cfg$run_analysis) {

  source("R/process_data.R", local=T)
  source("R/analysis.R", local=T)
  
} else if (cfg$run_sims) {
  
  # Set level sets
  source("R/levels.R", local=T)
  
  # Run simulation
  source("R/run.R", local=T)
  
}

if (cfg$run_process_results) {
  
  # Tables and figures
  source("R/process_results.R", local=T)
  # source("R/process_results2.R", local=T) # !!!!!
  
}

if (cfg$run_dqa) {
  
  # Data quality assurance
  source("R/dqa.R", local=T)
  
}
