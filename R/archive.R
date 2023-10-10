##############################.
##### VIZ: Others (temp) #####
##############################.

if (F) {
  
  # !!!!!
  {
    # True params: a_x=0.005, a_y=0.003, a_v=0.1/0.7, g_x=c(1.3,1.002), g_y=c(1.2,1.001), g_v=c(1.2,1.001), beta_x=1.5
    true_val <- 1.5
    x_window <- 0.2
    # v1 <- "cox_g_y2_est"
    # v2 <- "lik_F_g_y2_est"
    # v3 <- "lik_M_g_y2_est"
    v1 <- "cox_beta_x_est"
    v2 <- "lik_F_beta_x_est"
    v3 <- "lik_M_beta_x_est"
    r1 <- filter(sim$results, params=="10pct testing")
    r2 <- filter(sim$results, params=="70pct testing")
    # r1 <- filter(sim$results, n==1000 & max_time==100)
    x <- c(r1[[v1]], r1[[v2]], r1[[v3]],
           r2[[v1]], r2[[v2]], r2[[v3]])
    stats <- c("10pct, Cox", "10pct, Lik F", "10pct, Lik M",
               "70pct, Cox", "70pct, Lik F", "70pct, Lik M")
    df_plot <- data.frame(
      x = x,
      y = rep(0, length(x)),
      which = rep(factor(stats, levels=stats), each=round(length(x)/length(stats)))
    )
    ggplot(df_plot, aes(x=x, y=y, color=which)) +
      geom_vline(xintercept=true_val, alpha=0.5, linetype="dashed") +
      geom_jitter(width=0, height=1, alpha=0.3, size=3) +
      facet_wrap(~which, ncol=1, strip.position="left") + # scales="free_x"
      labs(y=NULL) +
      ylim(-2,2) +
      xlim(true_val-x_window,true_val+x_window) +
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            strip.text.y.left = element_text(angle=0),
            legend.position="none")
    
  }
  
  # r <- sim$results
  r500 <- filter(sim$results, n==500)
  r1000 <- filter(sim$results, n==1000)
  r2000 <- filter(sim$results, n==2000)
  ln <- c(nrow(r500), nrow(r1000), nrow(r2000))
  plot_data <- data.frame(
    x = c(r500$lik_beta_x_est, r1000$lik_beta_x_est, r2000$lik_beta_x_est),
    y = c(r500$cox_beta_x_est, r1000$cox_beta_x_est, r2000$cox_beta_x_est),
    n = c(rep(500, ln[1]), rep(1000, ln[2]), rep(2000, ln[3]))
  )
  ggplot(plot_data, aes(x=x, y=y)) +
    geom_abline(slope=1, intercept=0, color="orange") +
    geom_point(alpha=0.3) +
    geom_point(data=data.frame(x=1.5,y=1.5), color="forestgreen", size=3, alpha=0.7) +
    facet_wrap(~n) +
    labs(x="Likelihood model", y="Cox model")
  
  sim %>% summarize(
    coverage = list(
      list(name="cov_cox", truth=1.5, lower="cox_beta_x_ci_lo", upper="cox_beta_x_ci_hi"),
      list(name="cov_lik", truth=1.5, lower="lik_beta_x_ci_lo", upper="lik_beta_x_ci_hi")
    )
  )
  
  
  
  # # Read in simulation object
  # sim <- readRDS("../simba.out/sim_20210615.simba")
  
  # # Transform results
  # sim$results %<>% mutate(
  #   # est_hr_hiv = exp(est_hiv),
  #   # est_hr_art = exp(est_art),
  #   hr_hiv_lab = paste("HIV+ HR:",hr_hiv),
  #   hr_art_lab = paste("ART+ HR:",hr_art),
  #   method = ifelse(method=="ideal","Ideal",ifelse(method=="mi","MI","")),
  #   log_hr_hiv = log(hr_hiv),
  #   log_hr_art = log(hr_art)
  # )
  
  # # Plot of estimates (HIV)
  # # Export: 6" X 3"
  # ggplot(sim$results, aes(x=method, y=est_hiv, color=method)) +
  #   geom_point(alpha=0.2, size=2) +
  #   geom_hline(aes(yintercept=log(hr_hiv)), linetype="dotted") +
  #   facet_wrap(~hr_hiv_lab, ncol=2) +
  #   theme(legend.position="none")
  
  # # Plot of estimates (ART)
  # # Export: 6" X 3"
  # ggplot(sim$results, aes(x=method, y=est_art, color=method)) +
  #   geom_point(alpha=0.2, size=2) +
  #   geom_hline(aes(yintercept=log(hr_art)), linetype="dotted") +
  #   facet_wrap(~hr_art_lab, ncol=4) +
  #   theme(legend.position="none")
  
  
  
  
  
  # # HIV graph
  # # Export PDF 3x8
  # ggplot(sim$results, aes(x=method, y=est_hr_hiv, color=method)) +
  #   geom_point(alpha=0.2, size=2) +
  #   geom_hline(aes(yintercept=hr_hiv), linetype="dotted") +
  #   facet_wrap(~hr_hiv_lab, ncol=4) +
  #   theme(legend.position="none") +
  #   labs(title="Point estimates (50 simulation replicates per level)",
  #        x="Method", y="Estimated hazard ratio (HIV+ART-)")
  
  # # Coverage plots
  # summ <- sim %>% summary(
  #   coverage = list(
  #     name="cov_hiv", estimate="est_hiv", truth="log_hr_hiv", se="se_hiv"
  #   )
  # ) %>% mutate(
  #   hr_hiv_lab = paste("HIV+ HR:",hr_hiv),
  #   hr_art_lab = paste("ART+ HR:",hr_art),
  # )
  
  # ggplot(summ, aes(x=method, y=cov_hiv, color=method)) +
  #   geom_point(size=3) +
  #   geom_hline(aes(yintercept=0.95), linetype="dotted") +
  #   facet_wrap(~hr_hiv_lab, ncol=4) +
  #   theme(legend.position="none") +
  #   labs(title="Coverage (50 simulation replicates per level)",
  #        x="Method", y="95% CI Coverage (HIV+ART-)")
  
}



##############################.
##### Run dataset checks #####
##############################.

if (F) {
  
  # Setup
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(magrittr)
  library(data.table)
  library(survival)
  source("generate_dataset.R")
  source("perform_imputation.R")
  source("transform_dataset.R")
  source("run_analysis.R")
  source("helpers.R")
  p_sero_year <- convert_p_sero(list(male = list("1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,"26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005),female = list("1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,"26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0)))
  
  # Create "original" dataset
  dat_orig <- generate_dataset(
    num_patients = 10000,
    start_year = 2000,
    end_year = 2020,
    hazard_ratios = list("hiv"=1.7, "art"=1.4),
    p_sero_year = p_sero_year,
    p_death_year = p_death_year(1),
    u_mult = list("sero"=1, "death"=1)
  )
  
  # Create "imputed" dataset
  dat_imp <- perform_imputation(dat_orig, p_sero_year)
  
  # Create transformed datasets
  dat_cp_orig <- transform_dataset(dat_orig)
  dat_cp_imp <- transform_dataset(dat_imp)
  
  # !!!!! MICE TESTING: START
  # !!!!! TO DO
  # !!!!! MICE TESTING: END
  
  # Check 1: Compare overall seroconversion rates
  print(1-sum(is.na(dat_orig$sero_year))/nrow(dat_orig))
  print(1-sum(is.na(dat_imp$sero_year))/nrow(dat_imp))
  dat_orig %>% filter(sero_year<2000) %>% nrow()
  dat_imp %>% filter(sero_year<2000) %>% nrow()
  dat_orig %>% filter(sero_year>=2000) %>% nrow()
  dat_imp %>% filter(sero_year>=2000) %>% nrow()
  
  # Check 2: Check to see if distribution of seroconversion years is similar
  # !!!!! orig is getting much higher pre-2000 seroconversion rate
  ggplot(data.frame(
    sero_year = c(dat_orig$sero_year,dat_imp$sero_year),
    which = rep(c("original","imputed"),each=nrow(dat_orig))
  ), aes(x=sero_year, group=which, fill=factor(which))) +
    geom_histogram(color="white", bins=100) +
    geom_vline(xintercept = 2000, linetype="dotted") +
    facet_wrap(~which, ncol=2)
  
  # Check 3: Look at seroconversion years by "testing case"
  d3_orig <- dat_orig %>% group_by(case, sero_year) %>% summarize(count=n())
  d3_imp <- dat_imp %>% group_by(case, sero_year) %>% summarize(count=n())
  d3_orig$which <- "orig"
  d3_imp$which <- "imp"
  ggplot(
    rbind(d3_orig,d3_imp),
    aes(x=sero_year, y=count, group=which, fill=factor(case))
  ) +
    geom_bar(stat="identity", width=0.4) +
    geom_vline(xintercept = 2000, linetype="dotted") +
    facet_wrap(~case+which, ncol=2)
  
  # Check 4: Plot cascade status over time
  d4_orig <- dat_cp_orig %>% mutate(
    case = case_when(
      case==1 ~ "1. No testing data",
      case==2 ~ "2. Last test was neg",
      case==3 ~ "3. Neg test then pos test",
      case==4 ~ "4. First test was pos"
    )
  ) %>% group_by(start_year, casc_status, case) %>% summarize(num=n())
  ggplot(d4_orig, aes(x=start_year, y=num, color=casc_status)) +
    geom_line() + facet_wrap(~case, ncol=2)
  d4_imp <- dat_cp_imp %>% mutate(
    case = case_when(
      case==1 ~ "1. No testing data",
      case==2 ~ "2. Last test was neg",
      case==3 ~ "3. Neg test then pos test",
      case==4 ~ "4. First test was pos"
    )
  ) %>% group_by(start_year, casc_status, case) %>% summarize(num=n())
  ggplot(d4_imp, aes(x=start_year, y=num, color=casc_status)) +
    geom_line() + facet_wrap(~case, ncol=2)
  
  # Check 5: examine cascade status by case
  d5 <- dat_cp_orig %>% group_by(case, casc_status) %>%
    summarize(count=n())
  d5_imp <- dat_cp_imp %>% group_by(case, casc_status) %>% summarize(count=n())
  d5$imp <- d5_imp$count
  d5 %<>% rename(orig="count")
  d5b <- d5 %>% group_by(casc_status) %>% summarize(orig=sum(orig),imp=sum(imp))
  print(d5)
  print(d5b)
  print(d5 %>% filter(casc_status=="HIV+ART-") %>% mutate(diff=imp-orig))
  
  # Check 6: run analysis on both datasets; should yield comparable SEs
  run_analysis(dat_cp_orig, method="ideal")
  run_analysis(dat_cp_imp, method="ideal")
  
  # !!!!!!
  run_analysis(dat_cp_orig, method="censor")
  run_analysis(dat_cp_imp, method="censor")
  
  # Check 7: Examine case counts
  xtabs(~case, data=dat_orig)
  
  # Check 8: Examine death counts
  xtabs(~died, data=dat_orig)
  
  # Check 9: Compare death counts by casc_status
  # !!!!! This illustrates the source of the bias
  xtabs(~casc_status, data=filter(dat_cp_orig, died==1))
  xtabs(~casc_status, data=filter(dat_cp_imp, died==1))
  
  # Check 10: Compare death counts by case+casc_status
  xtabs(~case+casc_status, data=filter(dat_cp_orig, died==1))
  xtabs(~case+casc_status, data=filter(dat_cp_imp, died==1))
  
  # Check 11: examine datasets manually for abnormalities
  # write.table(dat_orig, file="dat_orig.csv", sep=",", row.names=FALSE)
  # write.table(dat_imp, file="dat_imp.csv", sep=",", row.names=FALSE)
  # write.table(dat_cp_orig, file="dat_cp_orig.csv", sep=",", row.names=FALSE)
  # write.table(dat_cp_imp, file="dat_cp_imp.csv", sep=",", row.names=FALSE)
  
}



###########################.
##### Mini-simulation #####
###########################.

if (F) {
  
  n <- 100000
  
  # x is our binary exposure variable
  prob_x <- 0.2
  x <- rbinom(n=n, size=1, prob=prob_x)
  
  # y is our binary outcome, correlated with x
  # rho_xy, the correlation between x and y, is the outcome of interest
  rho_xy <- 0.8
  y <- ifelse(runif(n)<rho_xy, x, rbinom(n=n, size=1, prob=prob_x))
  print(cor(x,y))
  
  est_corr <- c()
  rho_zx_vec <- seq(0,1,0.1)
  for (rho_zx in rho_zx_vec) {
    
    # z is a covariate correlated directly with x
    # Within this loop, we test different values of rho_zx, the correlation
    #     between x and z
    z <- ifelse(runif(n)<rho_zx, x, rbinom(n=n, size=1, prob=prob_x))
    
    # Impose missingness by deleting (100*k)% of the x values
    k <- 0.8
    x_trunc <- c(x[round(1:(n*(1-k)))],rep(NA,(n*k)))
    
    # Perform multiple imputation of x based on z; then estimate rho_xy
    est_cor_mi <- c()
    for (m in 1:10) {
      z_miss <- z[round((n*(1-k)+1):n)]
      x_imp <- ifelse(runif(n*k)<rho_zx, z_miss,
                      rbinom(n=(n*k), size=1, prob=prob_x))
      x_new <- c(x[round(1:(n*(1-k)))],x_imp)
      est_cor_mi <- c(est_cor_mi,cor(x_new,y))
    }
    est_corr <- c(est_corr, mean(est_cor_mi))
    
  }
  
  library(ggplot2)
  ggplot(
    data.frame(x=rho_zx_vec,y=est_corr),
    aes(x=x,y=y)) +
    geom_point() +
    geom_hline(yintercept=0.8, linetype="dotted") +
    labs(x="Correlation between x and z", y="Estimated rho_xy")
  
}



###############################################.
##### Old one_simulation() code (partial) #####
###############################################.

if (F) {
  
  dat2 <- impose_missingness(dat) # !!!!!


  # !!!!! Testing
  # C <- list(num_patients=10, start_year=2000, end_year=2001, m=5)
  # C <- list(num_patients=5000, start_year=2000, end_year=2002, m=5)
  # L <- list(method="mi", hr_hiv=1.4, hr_art=0.7)

  # Generate baseline data
  # !!!!! This is generated as a "true cohort" rather than an "open cohort"
  dat_baseline <- generate_data_baseline(
    num_patients = C$num_patients,
    start_year = C$start_year
  )

  # Set parameters
  params <- list(
    alpha0=-4,  alpha1=0.1,  alpha2=0.05,  alpha3=0, # alpha3=0.2
    beta0=-4,   beta1=0.1,   beta2=0.05,   beta3=0,  # beta3=0.2
    eta0=-4,    eta1=0.1,    eta2=0.05,    eta3=0,   # eta3=0.2
    gamma0=-4,  gamma1=0.1,  gamma2=0.05,  gamma3=0, # gamma3=0.2
    psi1=L$hr_hiv,
    psi2=L$hr_art
  )

  # Generate event data
  # !!!!! For now, all patients are HIV- at baseline
  dat_events <- apply(dat_baseline, MARGIN=1, function(r) {
    generate_data_events(
      id = r[["id"]],
      b_age = r[["b_age"]],
      sex = r[["sex"]],
      u = r[["u"]],
      start_year = C$start_year,
      end_year = C$end_year,
      baseline_status = NA,
      params = params
    )
  })
  attr(dat_events, "end_year") <- C$end_year

  # Take m samples from the posterior
  # theta_m <- posterior_param_sample(fit=fit, size=C$m)
  # !!!!! Temp: START
  {
    psi1_psample <- rnorm(C$m, mean=L$hr_hiv, sd=0.1)
    psi2_psample <- rnorm(C$m, mean=L$hr_art, sd=0.1)
    theta_m <- list()
    for (i in 1:C$m) {
      theta_m[[i]] <- params
      theta_m[[i]]$psi1 <- psi1_psample[i]
      theta_m[[i]]$psi2 <- psi2_psample[i]
    }
  }
  # !!!!! Temp: END

  if (L$method=="ideal") {
    # Transform data and run Cox PH analysis
    dat_cp <- transform_dataset(
      dat_baseline = dat_baseline,
      dat_events = dat_events
    )
    results <- run_analysis(dat_cp=dat_cp)
    # print(exp(results$est_hiv)); print(exp(results$est_art)); # !!!!!
  }

  if (L$method=="mi") {

    # Perform MI on second dataset
    dat_events_mi <- list()
    for (i in 1:C$m) {
      dat_events_mi[[i]] <- perform_imputation(
        dat_baseline = dat_baseline,
        dat_events = dat_events,
        theta_m = theta_m[[i]]
      )
    }

    results_mi <- list()
    for (i in 1:C$m) {
      # Transform data and run Cox PH analysis
      dat_cp <- transform_dataset(
        dat_baseline = dat_baseline,
        dat_events = dat_events_mi[[i]]
      )
      results_mi[[i]] <- run_analysis(dat_cp=dat_cp)
      # print(paste("Replicate:",i)) # !!!!!
      # print(paste("HIV:", exp(results_mi[[i]]$est_hiv))) # !!!!!
      # print(paste("ART:", exp(results_mi[[i]]$est_art))) # !!!!!
    }

    # Combine MI estimates using "Rubin's rules"
    v <- function(results_mi, attr) {
      sapply(results_mi, function(r) { r[[attr]] })
    }
    est_hiv <- v(results_mi, "est_hiv")
    est_art <- v(results_mi, "est_art")
    var_hiv <- (v(results_mi, "se_hiv"))^2
    var_art <- (v(results_mi, "se_art"))^2
    results <- list(
      est_hiv = mean(est_hiv),
      se_hiv = sqrt( mean(var_hiv) + (1+(1/C$m))*var(est_hiv) ),
      est_art = mean(est_art),
      se_art = sqrt( mean(var_art) + (1+(1/C$m))*var(est_art) )
    )

  }

  return (res)
  
}



#################################################.
##### Old posterior_param_sample() function #####
#################################################.

if (F) {
  
  #' Take a posterior sample of parameter vector theta
  #' 
  #' @param fit Stan model fit returned by fit_stan()
  #' @param size Number of samples to take (m); one per imputation
  #' @return A list of parameter samples
  posterior_param_sample <- function(fit, size) {
    
    n.iter <- nrow(fit[[1]])
    chains <- sample(c(1,2), size=size, replace=TRUE)
    rows <- sample(c(1:n.iter), size=size, replace=FALSE)
    alpha0 <- psi1 <- psi2 <- c()
    for (i in 1:C$m) {
      alpha0 <- c(alpha0, fit[[chains[i]]][[rows[i],"alpha0"]])
      psi1 <- c(psi1, fit[[chains[i]]][[rows[i],"psi1"]])
      psi2 <- c(psi2, fit[[chains[i]]][[rows[i],"psi2"]])
    }
    
    return(list(
      alpha0=alpha0, psi1=psi1, psi2=psi2
    ))
    
  }
  
}



#######################################.
##### Old run_analysis() function #####
#######################################.

if (F) {
  
  #' Run Cox PH analysis
  #'
  #' @param dat_cp A dataset returned by transform_dataset()
  #' @param options Placeholder; currently unused
  #' @return A list containing the following:
  #'     est_hiv: point estimate of HIV+ART- exposure coefficient
  #'     se_hiv: standard error of HIV+ART- exposure coefficient
  #'     est_art: point estimate of HIV+ART+ exposure coefficient
  #'     se_art: standard error of HIV+ART+ exposure coefficient
  
  run_analysis <- function(dat_cp, options=list()) {
    
    # Create "censor" dataset
    if (!is.null(options$method) && options$method=="censor") {
      
      # Add an ID row
      dat_cp <- cbind("obs_id"=c(1:nrow(dat_cp)),dat_cp)
      
      # Exclude patients with no testing data
      dat_cp %<>% filter(case!=1)
      
      # Exclude observation time prior to the first test
      dat_cp %<>% filter(start_year>=first_test)
      
      # Exclude observation time after the last negative test for case 2
      dat_cp %<>% filter(
        !(replace_na(case==2 & start_year>last_neg_test,FALSE))
      )
      
      # !!!!! Check censoring manually
      
    }
    
    # Fit time-varying Cox model ("ideal")
    # !!!!! Also try with robust SEs
    fit <- coxph(
      Surv(start_time, end_time, y) ~ factor(casc_status) + age + sex +
        cluster(id),
      data = dat_cp
    )
    summ <- summary(fit)$coefficients
    
    # # !!!!! Fit a logistic discrete survival model
    # fit2 <- glm(
    #   y ~ factor(casc_status) + age + sex,
    #   data = dat_cp,
    #   # start = c(-0.1,-0.1,-0.1,-0.1,-0.1),
    #   family = "binomial"
    #   # family = binomial(link="log")
    # )
    # summ <- summary(fit2)$coefficients
    
    results <- list(
      est_hiv = summ["factor(casc_status)HIV+ART-","coef"],
      se_hiv = summ["factor(casc_status)HIV+ART-","se(coef)"],
      est_art = summ["factor(casc_status)HIV+ART+","coef"],
      se_art = summ["factor(casc_status)HIV+ART+","se(coef)"]
      # est_hiv = summ["factor(casc_status)HIV+ART-","Estimate"],
      # se_hiv = summ["factor(casc_status)HIV+ART-","Std. Error"],
      # est_art = summ["factor(casc_status)HIV+ART+","Estimate"],
      # se_art = summ["factor(casc_status)HIV+ART+","Std. Error"]
    )
    
    return(results)
    
  }
  
}



############################################.
##### Old transform_dataset() function #####
############################################.

if (F) {
  
  #' Transform dataset into "counting process" format
  #'
  #' @param dat_baseline A dataset returned by generate_data_baseline()
  #' @param dat_events A dataset returned by generate_data_events()
  #' @return A dataset in "counting process" format for Cox PH analysis
  
  transform_dataset <- function(dat_baseline, dat_events) {
    
    # Extract variables
    start_year <- attr(dat_baseline, "start_year")
    end_year <- attr(dat_events, "end_year")
    I <- nrow(dat_baseline)
    
    # Data transformation
    dat_baseline$start_time <- 0
    dat_baseline$end_time <- sapply(dat_events, function(d) { d$T_i })
    dat_baseline$y_copy <- sapply(dat_events, function(d) { d$y[d$T_i] })
    
    # Put dataset into "counting process" format
    dat_cp <- survSplit(
      formula = Surv(start_time, end_time, y_copy) ~.,
      data = dat_baseline,
      cut = c(1:(12*(end_year-start_year)))
    )
    
    # Convert dat_events to a dataframe and attach to dat_cp
    # cbind is functioning as an inner join since both dataframes are sorted
    df_ev <- as.data.frame(rbindlist(dat_events))
    df_ev %<>% filter(y!=9)
    df_ev %<>% subset(select=-id)
    dat_cp %<>% cbind(df_ev)
    
    # Create exposure variable
    dat_cp %<>% mutate(
      casc_status = case_when(
        x==0 ~ "HIV-",
        x==1 & z==0 ~ "HIV+ART-",
        x==1 & z==1 ~ "HIV+ART+"
      ),
      age = b_age + start_time/12
    )
    
    return(dat_cp)
    
  }
  
}



#############################################.
##### Old perform_imputation() function #####
#############################################.

if (F) {
  
  #' Perform imputation on a dataset with missingness
  #'
  #' @param dat_baseline A dataset returned by generate_data_baseline()
  #' @param dat_events A dataset returned by generate_data_events()
  #' @param theta_m An posterior draw of the parameters
  #' @return dat_events, but with missing values in X imputed
  
  perform_imputation <- function(dat_baseline, dat_events, theta_m) {
    
    p <- theta_m
    
    # !!!!! Make sure we are memoising within perform_imputation
    
    # # !!!!! TESTING
    # dat_events_backup <- dat_events
    # dat_events <- dat_events[1:3]
    
    # Perform imputation for each patient
    x_imputed <- lapply(dat_events, function(de) {
      
      db_i <- dat_baseline[de$i,]
      
      # S_iX is the set of X values with positive posterior probability
      # The actual value assigned represents the number of zeros in X
      S_iX <- case_when(
        de$case == 1 ~ list(c(0:de$T_i)),
        de$case == 2 ~ list(c(de$last_neg_test:de$T_i)),
        de$case == 3 ~ list(c(de$last_neg_test:(de$first_pos_test-1))),
        de$case == 4 ~ list(c(0:(de$first_pos_test-1)))
      )[[1]]
      
      # Calculate component discrete hazards
      # Note: p and db_i are accessed globally
      p_it <- memoise(function(t) {
        expit(
          p$alpha0 + p$alpha1*db_i$sex + p$alpha2*(db_i$b_age+(t-1)/12) +
            p$alpha3*db_i$u
        )
      })
      q_it <- memoise(function(t,x) {
        min(0.99999, expit(
          p$gamma0 + p$gamma1*db_i$sex + p$gamma2*(db_i$b_age+(t-1)/12) +
            p$gamma3*db_i$u
        ) * exp(
          log(p$psi1)*x*(1-de$z[t]) +
            log(p$psi2)*x*de$z[t]
        ))
      })
      
      # In this block, we assign a probability to each possible value of S_iX
      # Note: d is accessed globally
      # Note: make sure these probabilities line up with those in
      #     generate_data_events.R and fit_stan.R
      probs <- sapply(S_iX, function(x) {
        
        # !!!!! Need to QA this
        
        if (de$case<=2) {
          
          if (x==de$last_neg_test) {
            P_X <- p_it(x+1)
          } else if (x %in% c((de$last_neg_test+1):(de$T_i-1))) {
            P_X_part <- prod(sapply(c((de$last_neg_test+1):x), function(s) {
              (1 - p_it(s))
            }))
            P_X <- P_X_part * p_it(x+1)
          } else if (x==de$T_i) {
            P_X <- prod(sapply(c((de$last_neg_test+1):de$T_i), function(s) {
              (1 - p_it(s))
            }))
          } else {
            stop("x is out of range; debug")
          }
          
        } else if (de$case>=3) {
          
          sum_p <- sum(sapply(c((de$last_neg_test+1):de$first_pos_test), p_it))
          P_X <- p_it(x+1) / sum_p
          
        }
        
        x_vec <- c(rep(0,x),rep(1,de$T_i-x))
        P_Y_part <- prod(sapply(c(1:(de$T_i-1)), function(s) {
          1 - q_it(s,x_vec[s])
        }))
        q_T_i <- q_it(de$T_i,x_vec[de$T_i])
        P_Y <- P_Y_part * ifelse(de$y[de$T_i]==1, q_T_i, 1-q_T_i)
        
        return(P_X*P_Y)
        
      })
      probs <- probs / sum(probs)
      
      if (round(sum(probs),6)!=1) {
        stop("S_iX probabilities don't sum to one; debug")
      }
      mult <- which(as.integer(rmultinom(n=1, size=1, prob=probs))==1)
      x_i_sample <- S_iX[mult]
      
      return (c(rep(0,x_i_sample),
                rep(1,de$T_i-x_i_sample),
                rep(9, sum(de$y==9))))
      
    })
    
    # Merge imputations back into dat_events
    dat_imputed <- dat_events
    for (i in 1:length(dat_events)) {
      dat_imputed[[i]]$x <- x_imputed[[i]]
    }
    
    return(dat_imputed)
    
  }
  
}



#########################.
##### Old MAIN code #####
#########################.

if (F) {
  
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
  
  # Dataset checks
  {
    
    # Check 4: Look at cascade status by year and "testing case"
    d4_orig <- dat_cp_orig %>% group_by(case, start_year) %>%
      summarize(num_hivpos_artneg=sum(casc_status=="HIV+ART-"))  
    d4_imp <- dat_cp_imp %>% group_by(case, start_year) %>%
      summarize(num_hivpos_artneg=sum(casc_status=="HIV+ART-"))  
    d4_orig$which <- "orig"
    d4_imp$which <- "imp"
    ggplot(
      rbind(d4_orig,d4_imp),
      aes(x=start_year, y=num_hivpos_artneg, group=which, fill=factor(case))
    ) +
      geom_bar(stat="identity", width=0.4) +
      facet_grid(rows=vars(case), cols=vars(which), scales="free_y")
    
  }
  
}



############################################.
##### Old dataset generating functions #####
############################################.

if (F) {
  
  #' Generate baseline data
  #'
  #' @param num_patients Number of patients in cohort
  #' @param start_year Start of cohort (Jan 1st, start_year)
  #' @return A dataframe, one row per patient, containing the following fields:
  #'     - id: patient ID variable
  #'     - b_age: baseline age (in completed years at start_year)
  #'     - birth_year: birth year (Jan 1st of year)
  #'     - sex: sex (0=female, 1=male)
  #'     - u: unmeasured "health behavior" variable
  #' @notes
  #'     - All dates are in CMC (Century Month Code) format, which is the number
  #'       of months since January 1st, 1900
  
  generate_data_baseline <- function(
    num_patients, start_year
  ) {
    
    # Generate baseline variables: patient_id, age, birth_year, sex, u
    # All patients assumed to be born on Jan 1st
    # Age represents number of completed years
    # All baseline variables represent values at Jan 1st of start_year
    # !!!!! Change the age distribution
    # !!!!! Add baseline status/testing data
    {
      id <- c(1:num_patients)
      b_age <- sample(1:80, size=num_patients, replace=TRUE)
      # birth_year <- start_year - b_age
      sex <- sample(c(0,1), size=num_patients, replace=TRUE)
      u <- rnorm(n=num_patients)
    }
    
    dat <- data.frame(id=id, b_age=b_age, sex=sex, u=u) # birth_year=birth_year
    attr(dat, "start_year") <- start_year
    attr(dat, "num_patients") <- num_patients
    
    return (dat)
    
  }
  
  
  #' Generate cohort events (seroconversion, ART initiation, testing, death)
  #'
  #' @param id Patient ID
  #' @param b_age Age of patient at start_year
  #' @param sex Sex of patient (0=female,1=male)
  #' @param u Latent health behavior variable
  #' @param start_year Start of cohort (Jan 1st, start_year)
  #' @param end_year End of cohort (Jan 1st, end_year)
  #' @param baseline_status Baseline cascade status; one of c("HIV-","HIV+ART-",
  #'     "HIV+ART+"); currently unused !!!!!
  #' @param params A list of Markov model parameters (alpha, beta, ...)
  #' @return A list containing the following:
  #'     - v: vector of testing indicators
  #'     - x: vector of serostatus indicators
  #'     - y: vector of outcome indicators
  #'     - z: vector of ART status indicators
  #'     - J: vector of outcome indicators
  #' @notes Much of this code mirrors code in fit_stan.R; ensure the two are in
  #'     sync with one another
  
  generate_data_events <- function(
    id, b_age, sex, u, start_year, end_year, baseline_status, params
  ) {
    
    p <- params
    
    # Set baseline variables
    x <- v <- z <- y <- c()
    x_last <- z_last <- 0
    
    # Sample events
    # Note: this code mirrors the MCMC code
    for (t in 1:(12*(end_year-start_year))) {
      
      if (length(y)==0 || max(y, na.rm=TRUE)==0) {
        
        # Seroconversion
        p_sero <- ifelse(x_last==1, 1, expit(
          p$alpha0 + p$alpha1*sex + p$alpha2*(b_age+(t-1)/12) + p$alpha3*u
        ))
        x <- c(x, rbinom(n=1, size=1, prob=p_sero))
        
        # Testing
        # !!!!! Add a condition s.t. patient doesn't get tested after POS test
        p_test <- expit(
          p$beta0 + p$beta1*sex + p$beta2*(b_age+(t-1)/12) + p$beta3*u
        )
        v <- c(v, rbinom(n=1, size=1, prob=p_test))
        
        # ART
        p_art <- ifelse(z_last==1, 1,
                        ifelse(x[length(x)]==0 || v[length(v)]==0, 0, expit(
                          p$eta0 + p$eta1*sex + p$eta2*(b_age+(t-1)/12) + p$eta3*u
                        ))
        )
        z <- c(z, rbinom(n=1, size=1, prob=p_art))
        
        # Outcome
        p_y <- min(0.99999, expit(
          p$gamma0 + p$gamma1*sex + p$gamma2*(b_age+(t-1)/12) + p$gamma3*u
        ) * exp(
          log(p$psi1)*x[length(x)]*(1-z[length(z)]) +
            log(p$psi2)*x[length(x)]*z[length(z)]
        ))
        if (p_y==0.99999) { warning("p_y>=0.99999") }
        y <- c(y, rbinom(n=1, size=1, prob=p_y))
        
        x_last <- x[length(x)]
        z_last <- z[length(z)]
        
      } else {
        
        # NA values coded as 9 (for Stan)
        v <- c(v, 9)
        x <- c(x, 9)
        y <- c(y, 9)
        z <- c(z, 9)
        
      }
      
    }
    
    # Add "testing case" to dataset
    #     Case 1: no testing data
    #     Case 2: most recent test was negative
    #     Case 3: negative test followed by a positive test
    #     Case 4: first test was positive
    T_i <- sum(v!=9)
    s <- as.integer(sum(v[1:T_i])>0)
    test_first <- s * min( (1:T_i) + T_i*(1-v[1:T_i]) )
    test_last <- s * max( (1:T_i)*v[1:T_i] )
    case <- ifelse(s==0, 1, ifelse(
      x[test_last]==0, 2, ifelse(
        x[test_first]==1, 4, 3
      )
    ))
    
    # Add last_neg_test and first_pos_test
    last_neg_test <- 0
    first_pos_test <- 0
    if (case==2) {
      last_neg_test <- test_last
    }
    if (case==3) {
      last_neg_test <- max( (1:T_i)*(v[1:T_i])*(1-x[1:T_i]) )
      first_pos_test <- min( (1:T_i) + T_i*(1-(x[1:T_i])*(v[1:T_i])) )
    }
    if (case==4) {
      first_pos_test <- test_first
    }
    
    # Calculate delta and delta*x
    miss <- rep(9,sum(x==9))
    if (case==1) {
      delta <- c(rep(0,T_i), miss)
    }
    if (case==2) {
      delta <- c(rep(1,last_neg_test), rep(0,T_i-last_neg_test), miss)
    }
    if (case==3) {
      delta <- c(rep(1,last_neg_test),
                 rep(0,first_pos_test-last_neg_test-1),
                 rep(1,T_i-first_pos_test+1),
                 miss)
    }
    if (case==4) {
      delta <- c(rep(0,first_pos_test-1), rep(1,T_i-first_pos_test+1), miss)
    }
    deltax <- delta*x
    deltax <- ifelse(deltax==81,9,deltax)
    
    # !!!!! Condense code when porting to Stan
    
    return(list(id=id, v=v, x=x, y=y, z=z, T_i=T_i, last_neg_test=last_neg_test,
                first_pos_test=first_pos_test, test_first=test_first,
                test_last=test_last, case=case, delta=delta, deltax=deltax))
    
  }
  
}



################################.
##### Old code from MAIN.R #####
################################.

if (F) {
  
  # Specify seroconversion conditional probabilities (discrete hazards)
  # For an individual of (exactly) age X-1, the number in the corresponding
  #     age bin represents the probability that the individual will
  #     seroconvert sometime in their Xth full year of life, given that they
  #     did not seroconvert by their (X-1)th full year of life.
  # !!!!! Change the actual numbers later based on AHRI cohort data
  p_sero_year <- convert_p_sero(list(
    male = list(
      "1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,
      "26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005
    ),
    female = list(
      "1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,
      "26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0
    )
  ))
  
}



##############################################.
##### Old code from perform_imputation() #####
##############################################.

if (F) {
  
  dataset$sero_year <- apply(
    # dataset$sero_year3 <- apply( # !!!!!
    X = dataset,
    MARGIN = 1,
    end_year = attributes(dataset)$end_year,
    FUN = function(r, end_year) {
      
      for (var in c("sex","died","death_year","last_neg_test","first_pos_test",
                    "birth_year", "case")) {
        assign(var, as.numeric(r[[var]]))
      }
      
      sex_mf <- ifelse(sex, "male", "female")
      end_year <- ifelse(died==1, death_year+1, end_year)
      
      # Case 1: no testing data
      if (case==1) {
        age_end <- end_year - birth_year
        
        # probs <- m_probs[[sex_mf]][[age_end]]
        # !!!!! TESTING
        if (died==0) {
          probs <- m_probs[[sex_mf]][[age_end]]
        } else {
          probs <- m_probs2[[sex_mf]][[age_end]]
        }
        
        pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
        if (pos==length(probs)) {
          return (NA)
        } else {
          return (birth_year+pos-1)
        }
      }
      
      # Case 2: most recent test was negative
      # Iteratively sample using discrete hazards
      if (case==2) {
        age_start <- last_neg_test - birth_year + 2
        age_end <- end_year - birth_year
        sero_year <- NA
        age <- age_start
        while (is.na(sero_year) && age<=age_end) {
          mult <- ifelse(died,hr_hiv_est,1) # !!!!!
          if (runif(1)<(p_sero_year[[sex_mf]][age])*mult) { # !!!!!
            sero_year <- birth_year + age - 1
          }
          age <- age + 1
        }
        return(sero_year)
      }
      
      # Case 3: negative test followed by a positive test
      # Sample from a multinomial with probs proportional to discrete hazards
      if (case==3) {
        age_start <- last_neg_test - birth_year + 2
        age_end <- first_pos_test - birth_year + 1
        probs <- p_sero_year[[sex_mf]][c(age_start:age_end)]
        pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
        return(last_neg_test + pos)
      }
      
      # Case 4: first test was positive
      # Sample from a multinomial with probs proportional to discrete hazards
      if (case==4) {
        age_start <- 1
        age_end <- first_pos_test - birth_year + 1
        probs <- p_sero_year[[sex_mf]][c(age_start:age_end)]
        pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
        return(birth_year + pos - 1)
      }
      
    }
  )
  
}



#####################.
##### Version 3 #####
#####################.

# Test JAGS data
if (F) {
  
  I <- 20
  dat <- list(
    I = I,
    J = c(),
    sex = sample(c(0,1), size=I, replace=TRUE),
    b_age = sample(c(30:50), size=I, replace=TRUE),
    v = matrix(NA, nrow=I, ncol=15),
    x = cbind(rep(0,I), matrix(NA, nrow=I, ncol=15)),
    y = matrix(NA, nrow=I, ncol=15),
    z = cbind(rep(0,I), matrix(0, nrow=I, ncol=15))
  )
  alpha0 <- -5;  alpha1 <- 0.3;  alpha2 <- 0.1;  alpha3 <- 0.3
  gamma0 <- -6;  gamma1 <- 0.1;  gamma2 <- 0.05;  gamma3 <- 0.1
  psi1 <- 0.4
  psi2 <- 0.2
  for (i in 1:I) {
    x_last <- 0
    for (j in 1:15) {
      
      p_sero <- expit(
        alpha0 + alpha1*dat$sex[i] + alpha2*(dat$b_age[i]+j-1) + alpha3*dat$u[i]
      )
      dat$x[i,j] <- rbinom(1, 1, ifelse(x_last==1, 1, p_sero))
      x_last <- dat$x[i,j]
      
      p_outcome <- min(1, expit(
        gamma0 + gamma1*dat$sex[i] + gamma2*(dat$b_age[i]+j-1) + gamma3*dat$u[i]
      ) * exp(
        log(psi1)*dat$x[i,j]*(1-dat$z[i,j]) + log(psi2)*dat$x[i,j]*dat$z[i,j]
      ))
      dat$y[i,j] <- rbinom(1, 1, p_outcome)
      
      if (dat$y[i,j]==1 || j==15) {
        dat$J <- c(dat$J, j)
        break
      }
      
    }
  }
}



#####################.
##### Version 2 #####
#####################.

if (F) {
  
  #' Format converter for p_sero_year
  #' @param p_sero_year number
  #' @return p_sero_year, but with individual years instead of buckets
  #' 
  convert_p_sero <- function(p_sero_year) {
    
    new_list <- list()
    for (sex in c("male", "female")) {
      p <- p_sero_year[[sex]]
      new_probs <- c(
        rep(p[["1"]],1), rep(p[["2-10"]],9), rep(p[["11-15"]],5),
        rep(p[["16-20"]],5), rep(p[["21-25"]],5), rep(p[["26-30"]],5),
        rep(p[["31-35"]],5), rep(p[["36-40"]],5), rep(p[["41-45"]],5),
        rep(p[["46-50"]],5), rep(0, 50)
      )
      new_list[[sex]] <- new_probs
    }
    
    return (new_list)
    
  }
  
  #' Return discrete hazards of death, by age
  #' @return A vector of discrete hazards, indexed by age
  #' 
  p_death_year <- function(mult) {
    
    # !!!!! Need to update these numbers
    probs <- c(
      rep(0.01, 9), # 1-9
      rep(0.002, 10), # 10-19
      rep(0.002, 10), # 20-29
      rep(0.004, 10), # 30-39
      rep(0.004, 10), # 40-49
      rep(0.01, 10), # 50-59
      rep(0.01, 10), # 60-69
      rep(0.03, 10), # 70-79
      rep(0.03, 10), # 80-89
      rep(0.1, 10), # 90-99
      rep(0.5, 10) # 100-109
    )
    
    return (mult*probs)
    
  }
  
  #' Construct multinomial probabilities
  #'
  #' @param p_sero_year A list of List of monthly seroconversion probabilities
  #'     (see simulation constants)
  #' @return A list of multinomial probabilities. The index of the list is the age
  #'     of an individual in the dataset. For an individual of (exactly) age 33,
  #'     the list value at index 33 will be a vector of length 34. This is a
  #'     vector of multinomial probabilities, where first entry is the prob that
  #'     the individual seroconverted by age 1 (MTCT), the second entry is the
  #'     prob that the individual seroconverted by age 2, etc. The last entry is
  #'     the prob that the individual never seroconverted.
  
  construct_m_probs <- function(p_sero_year) {
    
    pm <- p_sero_year$male
    pf <- p_sero_year$female
    
    p_sero_m <- list(c(pm[1],1-pm[1]))
    p_sero_f <- list(c(pf[1],1-pf[1]))
    
    for (age in 2:100) {
      
      # Set next item to previous item (without last prob)
      p_sero_m[[age]] <- p_sero_m[[age-1]][1:(age-1)]
      p_sero_f[[age]] <- p_sero_f[[age-1]][1:(age-1)]
      
      # Calculate and set next prob
      next_prob_m <- prod(1-pm[1:(age-1)]) * pm[age]
      next_prob_f <- prod(1-pf[1:(age-1)]) * pf[age]
      p_sero_m[[age]] <- c(p_sero_m[[age]], next_prob_m)
      p_sero_f[[age]] <- c(p_sero_f[[age]], next_prob_f)
      
      # Set prob of not seroconverting
      p_sero_m[[age]] <- c(p_sero_m[[age]], 1-sum(p_sero_m[[age]]))
      p_sero_f[[age]] <- c(p_sero_f[[age]], 1-sum(p_sero_f[[age]]))
      
    }
    
    return(list(
      male = p_sero_m,
      female = p_sero_f
    ))
    
  }
  
  # Generate probabilities for baseline serostatus
  m_probs <- construct_m_probs(p_sero_year)
  
  # Generate baseline variables: testing_prob
  # These will be the indices of the testing prob vector c(0,0.1,0.2,0.5,1)
  {
    multinom_probs_u0 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
    multinom_probs_u1 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
    d_testing_prob_u0 <- rmultinom(n=num_patients,size=1,prob=multinom_probs_u0)
    d_testing_prob_u1 <- rmultinom(n=num_patients,size=1,prob=multinom_probs_u1)
    d_testing_prob_u0 <- apply(d_testing_prob_u0, 2, function(x){which(x==1)})
    d_testing_prob_u1 <- apply(d_testing_prob_u1, 2, function(x){which(x==1)})
    d_testing_prob_index <- (1-d_u)*d_testing_prob_u0 + d_u*d_testing_prob_u1
    # testing_probs <- c(0.1,0.1,0.1,0.1,0.1)
    testing_probs <- c(0,0.1,0.2,0.5,1)
  }
  
  # Generate baseline seroconversion years
  # !!!!! Incorporate u ?
  d_sero_year <- sapply(d_patient_id, function(i) {
    sex_mf <- ifelse(d_sex[i], "male", "female")
    probs <- m_probs[[sex_mf]][[d_age[i]]]
    pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
    if (pos==length(probs)) {
      return (NA)
    } else {
      return (d_birth_year[i]+pos-1)
    }
  })
  
  # Old function parameters
  #' @param p_sero_year List of yearly discrete hazards of seroconversion, by sex
  #' @param p_death_year List of yearly discrete hazards of death
  #' @param u_mult List of hazard multipliers
  
  # Generate baseline seroconversion years
  # !!!!! Incorporate u ?
  d_sero_year <- sapply(d_patient_id, function(i) {
    sex_mf <- ifelse(d_sex[i], "male", "female")
    probs <- m_probs[[sex_mf]][[d_age[i]]]
    pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
    if (pos==length(probs)) {
      return (NA)
    } else {
      return (d_birth_year[i]+pos-1)
    }
  })
  
}



#####################.
##### Version 1 #####
#####################.

if (F) {
  
  # FN: Impute seroconversion dates
  # This function takes in a dataset and returns a vector of imputed seroconversion dates
  
  # !!!!! TO DO !!!!!
  # `s` values are placeholders; needs to reflect actual imputation model
  # Delete "local variables" that are not needed
  
  impute_sero_dates <- function(df, end_date) {
    
    # Set up vector to hold imputed dates
    s_vector <- c()
    
    # Generate imputed dates
    for (i in 1:nrow(df)) {
      
      # Set local variables
      dob <- df[i,"dob"]
      sex <- df[i,"sex"]
      alive <- df[i,"alive"]
      dod <- df[i,"dod"]
      last_test_neg <- df[i,"last_test_neg"]
      first_test_pos <- df[i,"first_test_pos"]
      case <- df[i,"case"]
      
      switch(case,
             
             # Case 1: First HIV test was pos
             "1" = {
               s <- round(dob+((first_test_pos-dob)/2))
             },
             
             # Case 2: 1+ neg tests followed by a pos test
             "2" = {
               s <- round(last_test_neg+((first_test_pos-last_test_neg)/2))
             },
             
             # Case 3: 1+ neg tests and no pos test (and is alive)
             "3" = {
               if (runif(1)<0.5) {
                 s <- NA
               } else {
                 s <- round(dob+((end_date-dob)/2))
               }
             },
             
             # Case 4: 1+ neg tests and no pos test (and is dead)
             "4" = {
               if (runif(1)<0.5) {
                 s <- NA
               } else {
                 s <- round(dob+((dod-dob)/2))
               }
             },
             
             # Case 5: no testing data
             "5" = {
               if (runif(1)<0.5) {
                 s <- NA
               } else {
                 s <- round(dob+((ifelse(is.na(dod),end_date,dod)-dob)/2))
               }
             }
             
      )
      
      s_vector <- c(s_vector,s)
      
    }
    
    # Return vector of imputed dates
    return (s_vector)
    
  }
  
  
  
  
  # FN: Transform dataset
  # This function takes a dataset and (for patients who seroconverted) splits each row into multiple rows, each one corresponding to a different HIV status; this is to facilitate analysis via a Cox PH model with a time-varying exposure
  # `df` is either a data frame returned by create_dataset_ideal() or an imputed dataset accessed via mice::complete(imputation_object, i), where imputation_object is returned by create_imputed_datasets()
  
  transform_dataset <- function(df, end_date) {
    
    # Create new data frame
    new_df <- data.frame(
      "patient_id" = integer(),
      "dob" = integer(),
      "sex" = integer(),
      "alive" = integer(),
      "dod" = integer(),
      "last_test_neg" = integer(),
      "first_test_pos" = integer(),
      "art_init" = integer(),
      "s" = integer(),
      "case" = integer(),
      "hiv_status" = integer(),
      "start_time" = integer(),
      "end_time" = integer(),
      "had_event" = integer()
    )
    
    # Loop through data frame rows and split into multiple rows
    for (i in 1:nrow(df)) {
      
      # Set local variables
      dob <- df[i,"dob"]
      dod <- df[i,"dod"]
      s <- df[i,"s"]
      art_init <- df[i,"art_init"]
      alive <- df[i,"alive"]
      
      # Patients who never seroconverted
      if (is.na(s)) {
        
        new_row_1 <- df[i,]
        new_row_1[1,"hiv_status"] <- 1
        new_row_1[1,"start_time"] <- dob
        new_row_1[1,"end_time"] <- ifelse(is.na(dod),end_date,dod)
        new_row_1[1,"had_event"] <- ifelse(alive==0,1,0)
        
        new_df[nrow(new_df)+1,] <- new_row_1
        
      }
      
      # Patients who seroconverted but never initiated ART
      if (!is.na(s) & is.na(art_init)) {
        
        new_row_1 <- df[i,]
        new_row_1[1,"hiv_status"] <- 1
        new_row_1[1,"start_time"] <- dob
        new_row_1[1,"end_time"] <- s
        new_row_1[1,"had_event"] <- 0
        
        new_row_2 <- df[i,]
        new_row_2[1,"hiv_status"] <- 2
        new_row_2[1,"start_time"] <- s
        new_row_2[1,"end_time"] <- ifelse(is.na(dod),end_date,dod)
        new_row_2[1,"had_event"] <- ifelse(alive==0,1,0)
        
        new_df[nrow(new_df)+1,] <- new_row_1
        new_df[nrow(new_df)+1,] <- new_row_2
        
      }
      
      # Patients who seroconverted and initiated ART
      if (!is.na(s) & !is.na(art_init)) {
        
        new_row_1 <- df[i,]
        new_row_1[1,"hiv_status"] <- 1
        new_row_1[1,"start_time"] <- dob
        new_row_1[1,"end_time"] <- s
        new_row_1[1,"had_event"] <- 0
        
        new_row_2 <- df[i,]
        new_row_2[1,"hiv_status"] <- 2
        new_row_2[1,"start_time"] <- s
        new_row_2[1,"end_time"] <- art_init
        new_row_2[1,"had_event"] <- 0
        
        new_row_3 <- df[i,]
        new_row_3[1,"hiv_status"] <- 3
        new_row_3[1,"start_time"] <- art_init
        new_row_3[1,"end_time"] <- ifelse(is.na(dod),end_date,dod)
        new_row_3[1,"had_event"] <- ifelse(alive==0,1,0)
        
        new_df[nrow(new_df)+1,] <- new_row_1
        new_df[nrow(new_df)+1,] <- new_row_2
        new_df[nrow(new_df)+1,] <- new_row_3
        
      }
      
    }
    
    # Remove rows where start_time == end_time
    new_df %<>% filter(start_time != end_time)
    
    # Return transformed data frame
    return (new_df)
    
  }
  
  
  
  # FN: Create imputed datasets
  # This function multiply imputes m datasets by leveraging the impute_sero_dates() function
  
  create_imputed_datasets <- function(d_reality, m) {
    
    # Useful links:
    # https://stats.stackexchange.com/questions/78632/multiple-imputation-for-missing-values
    # https://github.com/stefvanbuuren/mice/blob/master/R/mice.R
    # https://github.com/stefvanbuuren/mice/blob/master/R/sampler.R
    # https://github.com/stefvanbuuren/mice/blob/master/R/mice.impute.mean.R
    
    # Create custom imputation method
    # data2 is a copy of the original dataset
    # Function needs to be created in the global environment to work
    mice.impute.hivmi <- function(y, ry, x, wy = NULL, data2, ...) {
      # Returns a list of the imputed values
      return(impute_sero_dates(data2, (2019-1900)*12))
    }
    
    # Create `where` matrix
    mtx_where <- is.na(d_reality)
    s_position <- match("s",names(d_reality))
    mtx_where[,-s_position] <- FALSE
    
    # Create `predictorMatrix`
    id_position <- match("patient_id",names(d_reality))
    dim <- length(names(d_reality))
    mtx_predictor <- matrix(0, nrow=dim, ncol=dim)
    mtx_predictor[s_position,id_position] <- 1
    
    # Conduct imputation
    # The `where` argument specifies that only `s` should be imputed
    # The `predictorMatrix` argument ensures that all `s` values are imputed; in reality it is bypassed by the custom `mice.impute.hivmi` method.
    # The `remove.constant` argument prevents `s` from being removed from the "variables to be imputed" list (note that the `remove.constant` argument is not documented; I had to look in the MICE package source code to find it)
    # The custom `data2` argument is passed to the `mice.impute.hivmi` method, allowing us to access the full dataset when performing the imputation procedure
    imp <- mice(
      data = d_reality,
      method = "hivmi",
      m = m,
      print = F,
      maxit = 1,
      where = mtx_where,
      predictorMatrix = mtx_predictor,
      remove.constant = FALSE,
      data2 = d_reality
    )
    
    # Return imputation object
    return (imp)
    
  }
  
}
