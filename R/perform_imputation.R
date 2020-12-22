#' Perform imputation on a dataset with missingness
#'
#' @param dataset A dataset returned by generate_dataset()
#' @param p_sero_year A list of seroconversion probabilities returned by
#'     convert_p_sero()
#' @return The same dataset but with missing seroconversion dates imputed

perform_imputation <- function(dataset, p_sero_year) {
  
  # !!!!! The fundamental flaw is that the imputation model MUST include the outcome
  # !!!!! These are placeholders
  hr_hiv_est <- ifelse(exists("L"), L$hr_hiv, 1.7) # !!!!!
  
  # Remove unknown seroconversion info
  # dataset$sero_year2 <- dataset$sero_year # !!!!!
  dataset$sero_year <- NA
  
  m_probs <- construct_m_probs(p_sero_year)
  m_probs2 <- construct_m_probs(lapply(p_sero_year, function(x){x*hr_hiv_est})) # !!!!!
  
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
  
  # # !!!!! If condition holds, give imputed version
  # dataset %<>% mutate(
  #   sero_year = ifelse(case==2, sero_year3, sero_year2)
  # )
  # # !!!!!
  
  return(dataset)
  
}
