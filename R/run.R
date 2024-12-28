if (cfg$sim_run_or_update=="run") {
  
  run_on_cluster(
    
    first = {
      
      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$sim_num,
        parallel = cfg$parallelize,
        n_cores = cfg$sim_n_cores,
        stop_at_error = cfg$sim_stop_at_error,
        # batch_levels = c("surv_true"),
        # return_batch_id = T,
        # seed = 123, # !!!!!
        packages = cfg$pkgs
      )
      sim <- do.call(set_levels, c(list(sim), level_set))
      # if (!is.null(cfg$keep)) { sim %<>% set_levels(.keep=cfg$keep) }
      
      # Simulation script
      sim %<>% set_script(one_simulation)
      
    },
    
    main = { sim %<>% run() },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
} else if (cfg$sim_run_or_update=="update") {
  
  update_sim_on_cluster(
    
    first = {
      sim <- readRDS(paste0(cluster_config$dir,"/sim.rds"))
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = { sim %<>% update_sim() },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
}
