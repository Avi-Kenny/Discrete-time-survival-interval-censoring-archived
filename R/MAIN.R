# Main config
# Sys.setenv(SLURM_CPUS_PER_TASK=5) # For running locally
# Sys.setenv(avi_samp_size=20000) # For running locally
# Sys.setenv(avi_maxit=3000) # For running locally
# Sys.setenv(avi_r=2) # For running locally
# Sys.setenv(avi_reltol=1e-5) # For running locally
cfg <- list(
  run_sims = F, # !!!!!
  run_analysis = T, # !!!!!
  run_process = F,
  sim_level_set = "level_set_1",
  sim_run_or_update = "run",
  sim_num = 1000,
  sim_parallel = F,
  # sim_n_cores = 50,
  sim_n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
  sim_stop_at_error = F,
  model_version = 18
)

# Secondary config
source("R/config.R", local=T)

# Load SimEngine + functions
{
  library(SimEngine)
  source("R/helpers.R", local=T)
  source("R/one_simulation.R", local=T)
  source("R/generate_data.R", local=T)
  source("R/likelihood.R", local=T)
}

if (cfg$run_analysis) {
  
  set.seed(1)
  source("R/analysis.R", local=T)
  
} else if (cfg$run_sims) {
  
  # Set level sets
  source("R/levels.R", local=T)
  
  # Run simulation
  source("R/run.R", local=T)
  
} else if (cfg$run_process) {
  
  # Tables and figures
  source("R/process.R", local=T)
  
}
