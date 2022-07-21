
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
      g_x = c(log(1.3),log(1.002)),
      g_y = c(log(1.2),log(1.001)),
      beta = log(1.5)
    )
  )
  
  # Run Cox PH model
  model <- coxph(
    Surv(t_start, t_end, y) ~ w_sex + w_age + x + cluster(id),
    data = dat
  )
  print("True params: c(1.2, 1.001, 1.5)")
  print(summary(model)$conf.int)
  
})
