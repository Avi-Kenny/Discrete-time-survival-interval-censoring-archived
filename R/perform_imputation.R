#' Perform imputation on a dataset with missingness
#'
#' @param dataset A dataset returned by generate_dataset()
#' @param p_sero_year A list of seroconversion probabilities returned by
#'     convert_p_sero()
#' @return The same dataset but with missing seroconversion dates imputed

perform_imputation <- function(dataset, p_sero_year) {
  
  m_probs <- construct_m_probs(p_sero_year)
  
  dataset$sero_year <- sapply(c(1:nrow(dataset)), function(i) {
    
    d <- dataset[i,]
    sex_mf <- ifelse(d$sex, "male", "female")
    end_year <- ifelse(d$died==1, d$death_year, attributes(dataset)$end_year)
    
    # !!!!! When generating data, make sure everyone falls into the proper category (particularly for baseline data)
    # !!!!! Then, modify run_analysis() accordingly
    
    # Case 1: no testing data
    if (is.na(d$last_neg_test) & is.na(d$first_pos_test)) {
      age_end <- end_year - d$birth_year
      probs <- m_probs[[sex_mf]][[age_end]]
      sero_age <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
      if (sero_age==length(probs)) {
        return(NA)
      } else {
        return(end_year-(age_end-sero_age))
      }
    }
    
    # Case 2: most recent test was negative
    # Iteratively sample using discrete hazards
    if (!is.na(d$last_neg_test) & is.na(d$first_pos_test)) {
      age_start <- d$last_neg_test - d$birth_year + 1
      age_end <- end_year - d$birth_year
      sero_year <- NA
      age <- age_start
      while (is.na(sero_year) && age<=age_end) {
        if (runif(1)<p_sero_year[[sex_mf]][age]) {
          sero_year <- d$birth_year + age
        }
        age <- age + 1
      }
      return(sero_year)
    }
    
    # Case 3: negative test followed by a positive test
    # Sample from a multinomial with probs proportional to discrete hazards
    if (!is.na(d$last_neg_test) & !is.na(d$first_pos_test)) {
      age_start <- d$last_neg_test - d$birth_year
      age_end <- d$first_pos_test - d$birth_year
      probs <- p_sero_year[[sex_mf]][c(age_start:age_end)]
      pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
      sero_age <- c(age_start:age_end)[pos]
      return(d$first_pos_test-(age_end-sero_age))
    }
    
    # Case 4: first test was positive
    # Sample from a multinomial with probs proportional to discrete hazards
    if (is.na(d$last_neg_test) & !is.na(d$first_pos_test)) {
      age_end <- d$first_pos_test - d$birth_year
      probs <- m_probs[[sex_mf]][[age_end]]
      probs <- probs[1:(length(probs)-1)]
      sero_age <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
      return(d$first_pos_test-(age_end-sero_age))
    }
    
  })
  
  return(dataset)
  
}
