
###################.
##### Recycle #####
###################.

# QA checks (consider recycling)
if (FALSE) {
  
  # This should be zero
  dataset %>% filter(baseline_status!="HIV-" & !is.na(last_neg_test)) %>%
    nrow()
  
  # !!!!! Add more QA checks
  
}



################################.
##### Old code from MAIN.R #####
################################.

if (FALSE) {
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

if (FALSE) {
  
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
{
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



#####################.
##### Version 1 #####
#####################.

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
