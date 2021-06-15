#' Perform imputation on a dataset with missingness
#'
#' @param dat_baseline A dataset returned by generate_data_baseline()
#' @param dat_events A dataset returned by generate_data_events()
#' @param theta_m An posterior draw of the parameters
#' @return dat_events, but with missing values in X imputed

perform_imputation <- function(dat_baseline, dat_events, theta_m) {
  
  p <- theta_m
  db <- dat_baseline
  
  # !!!!! Make sure we are memoising within perform_imputation
  
  # # !!!!! TESTING
  # dat_events_backup <- dat_events
  # dat_events <- dat_events[1:3]
  
  # Perform imputation for each patient
  x_imputed <- lapply(dat_events, function(d) {
    
    # S_iX is the set of X values with positive posterior probability
    # The actual value assigned represents the number of zeros in X
    S_iX <- case_when(
      d$case == 1 ~ c(0:d$T_i),
      d$case == 2 ~ c(d$last_neg_test:d$T_i),
      d$case == 3 ~ c(d$last_neg_test:(d$first_pos_test-1)),
      d$case == 4 ~ c(0:(d$first_pos_test-1))
    )
    
    # Calculate component seroconversion discrete hazards
    p_it <- memoise(function(t) {
      expit(
        p$alpha0 + p$alpha1*db$sex + p$alpha2*(db$b_age+(t-1)/12) +
          p$alpha3*db$u
      )
    })
    
    # In this block, we assign a probability to each possible value of S_iX
    # Note: d is accessed globally
    # Note: make sure these probabilities line up with those in
    #     generate_data_events.R and fit_stan.R
    probs <- sapply(S_iX, function(x) {
      
      if (case<=2) {
        
        # !!!!! CONTINUE
        
        if (x==d$last_neg_test) {
          P_X <- p_it(x+1)
        } else if (x %in% c((d$last_neg_test+1):(d$T_i-1))) {
          P_X_part <- prod(sapply(c((d$last_neg_test+1):x), function(s) {
            (1 - p_it(s))
          }))
          P_X <- P_X_part * p_it(x+1)
        } else if (x==d$T_i) {
          P_X_part <- prod(sapply(c((d$last_neg_test+1):d$T_i), function(s) {
            (1 - p_it(s))
          }))
        } else {
          stop("x is out of range; debug")
        }
        
      } else if (case>=3) {
        
        P_X <- 999
        
      }
      
      P_Y <- 999
      
      return(P_X*P_Y)
      
    })
    probs <- probs / sum(probs)
    
    if (sum(probs)!=1) { stop("S_iX probabilities don't sum to one; debug") } # !!!!!
    mult <- rmultinom(n=1, size=1, prob=probs)
    x_i_sample <- S_iX[mult]
    
    # Break off NAs (9s)
    
    
  })
  
  # Merge imputations back into dat_events
  dat_imputed <- dat_events
  for (i in 1:length(dat_events)) {
    dat_events[[i]]$x <- x_imputed[[i]]
  }
  
  return(dat_imputed)
  
  
  
  
  
  
  
  
  
  # !!!!! Recycle all code below
  
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
