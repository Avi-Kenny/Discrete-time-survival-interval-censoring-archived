
##############################################.
##### Generate dataset and fit Cox model #####
##############################################.

system.time({
  
  # Runtimes (seconds)
  # 2 (n=500), 5 (n=1000), 18 (n=2000), 65 (n=4000)
  
  # Generate data
  dat <- generate_data(
    n = 1000,
    max_time = 100,
    params = list(
      gamma = c(log(1.3),log(1.002)),
      xi = c(log(1.2),log(1.001)),
      beta_x = log(1.5)
    )
  )
  
  # Run Cox PH model
  model <- coxph(
    Surv(t_start, t_end, y) ~ z_sex + z_age + x + cluster(id),
    data = dat
  )
  print("True params: c(1.2, 1.001, 1.5)")
  print(summary(model)$conf.int)
  
})



