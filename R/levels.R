# # Set global constants
# C <- list(alpha_1=0.5)

# Set simulation levels
if (cfg$run_sims && Sys.getenv("sim_run") %in% c("first", "")) {
  
  # L <- list(n=300,max_time=70,params=list(a_x=log(0.005), a_y=log(0.003), a_v=log(0.7),a_z=log(0.004),g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),g_v=c(log(1.2),log(1.1)), g_z=c(log(1.2),log(1.1)),beta_x=log(1.5), beta_z=log(0.6)))
  
  level_sets <- list()
  
  # Simulation 1: basic
  # Figures: ...
  level_sets[["level_set_1"]] <- list(
    n = 500,
    max_time = 20,
    model_version = 7,
    params = list(
      "param set 1" = lapply(list(
        a_x = 0.05,
        g_x = c(1.3,1.2),
        t_x = 1,
        a_s = 0.05,
        g_s = c(2,1.5),
        t_s = 1,
        beta_x = 1.5,
        a_y = 0.03,
        g_y = c(1.2,1.1),
        t_y = 1,
        a_v = 0.7,
        g_v = c(1.2,1.1)
      ), log)
    )
  )
  
  level_set <- level_sets[[cfg$sim_level_set]]
  
  # if (cfg$sim_level_set=="asdf") { cfg$keep = c(1:3,7:9,16:18,22:24) }
  
}
