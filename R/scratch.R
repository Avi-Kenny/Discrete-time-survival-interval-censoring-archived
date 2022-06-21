
# Generate an incomplete dataset from a complete dataset
# Probability of being assigned to a "case" is independent of all other vars
impose_missingness <- function(dat) {
  
  # Generate the "testing case" variable
  #     Case 1: no testing data
  #     Case 2: most recent test was negative
  #     Case 3: negative test followed by a positive test
  #     Case 4: first test was positive
  n <- attr(dat, "n")
  # case <- sample(c(1:4), size=n, replace=T)
  dat$case <- rep(NA, nrow(dat))
  dat$test_reg <- rep(NA, nrow(dat))
  
  # Loop through patients
  for (i in c(1:n)) {
    
    # Extract variables from dataset
    rows <- which(dat$id==i)
    x_i <- dat[rows,"x"]
    n_obs_i <- length(rows)
    max_time <- attr(dat, "max_time")
    
    # Assign to one of three testing regimens:
    #   Reg 1: 0 tests (prob 0.1)
    #   Reg 2: ceil(max_time/8) tests (prob 0.45)
    #   Reg 3: ceil(max_time/4) tests (prob 0.45)
    test_reg <- sample(c(1:3), size=1, prob=c(0.1,0.45,0.45))
    if (test_reg==2) {
      num_tests <- ceiling(max_time/8)
    } else if (test_reg==3) {
      num_tests <- ceiling(max_time/4)
    } else {
      num_tests <- 0
    }
    tests <- sort(sample(c(1:max_time), size=num_tests))
    tests <- tests[tests<=n_obs_i]
    if (length(tests)==0) {
      case <- 1 # No testing data
      x_i <- rep(NA, n_obs_i)
    } else {
      if (x_i[tests[1]]==1) {
        case <- 4 # First test POS
        x_i[c(1:round(tests[1]-1))] <- NA
      } else if (x_i[tests[length(tests)]]==0) {
        case <- 2 # most recent test NEG
        if (n_obs_i>tests[length(tests)]) {
          x_i[c(1:n_obs_i)>tests[length(tests)]] <- NA
        }
      } else {
        case <- 3 # NEG test then POS test
        ind_neg <- max(tests[x_i[tests]==0])
        ind_pos <- min(tests[x_i[tests]==1])
        if (ind_pos-ind_neg>1) {
          x_i[c(round(ind_neg+1):round(ind_pos-1))] <- NA
        }
      }
    }
    
    dat[rows,"x"] <- x_i
    dat[rows,"test_reg"] <- test_reg
    dat[rows,"case"] <- case
    
    # # Set new X variable value based on case
    # if (case[i]==1) {
    #   x_new <- rep(NA, n_obs_i)
    # } else if (case[i]==2) {
    #   
    #   x_new <- 999 # !!!!!
    # } else if (case[i]==3) {
    #   x_new <- 999 # !!!!!
    # } else if (case[i]==4) {
    #   x_new <- 999 # !!!!!
    # }
    
    # # Add case variable to dataframe
    # dat[rows,"case"] <- case[i]
    
  }
  
  return (dat)
  
}

