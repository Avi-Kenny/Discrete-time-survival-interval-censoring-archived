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
