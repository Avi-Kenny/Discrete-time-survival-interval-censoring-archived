#################################################.
##### Old model code archived on 2024-11-02 #####
#################################################.

if (F) {
  
  # Create a natural cubic spline basis
  if (which=="age (0-100), 4DF") {
    grid <- scale_age(seq(0,100, length.out=500))
    k <- scale_age(c(0,25,50,75,100))
  } else if (which=="age (13,20,30,60,90)") {
    grid <- scale_age(seq(13,90, length.out=500))
    k <- scale_age(c(13, 20, 30, 60, 90))
  } else if (which=="age (13,30,60,75,90)") {
    grid <- scale_age(seq(13,90, length.out=500))
    k <- scale_age(c(13, 30, 60, 75, 90))
  } else if (which=="age (13,32,52,71,90)") {
    grid <- scale_age(seq(13,90, length.out=500))
    k <- scale_age(seq(13,90, length.out=5))
  } else if (which %in% c("age (13,28,44,60,75)", "age (13,28,44,60,75) +i")) {
    grid <- scale_age(seq(13,75, length.out=500))
    k <- scale_age(round(seq(13,75, length.out=5)))
  } else if (which=="year (00,05,10,15,20)") {
    grid <- scale_year(seq(2000,2022, length.out=500))
    k <- scale_year(seq(2000,2020, length.out=5))
  } else if (which %in% c("year (10,13,17,20,23)", "year (10,13,17,20,23) +i")) {
    grid <- scale_year(seq(2010,2023, length.out=500))
    k <- scale_year(seq(2010,2023, length.out=5))
  } else if (which %in% c("year (10,13,16,19,22)", "year (10,13,16,19,22) +i")) {
    grid <- scale_year(seq(2010,2022, length.out=500))
    k <- scale_year(seq(2010,2022, length.out=5))
  } else if (which=="age (13,20,30,40,60)") {
    grid <- scale_age(seq(13,60, length.out=500))
    k <- scale_age(c(13,20,30,40,60))
  } else if (which=="year (17,...,22)") {
    grid <- scale_year(seq(2017,2022, length.out=500))
    k <- scale_year(seq(2010,2022, length.out=5))
  }
  
  # Construct spline bases
  b2 <- construct_basis("age (13,20,30,60,90)")
  b3 <- construct_basis("age (13,30,60,75,90)")
  b4 <- construct_basis("year (00,05,10,15,20)",
                        window_start=window_start, window_end=window_end)
  b5 <- construct_basis("year (10,13,17,20,23)",
                        window_start=window_start, window_end=window_end)
  b6 <- construct_basis("age (13,28,44,60,75)")
  b7 <- construct_basis("year (10,13,17,20,23) +i",
                        window_start=window_start, window_end=window_end)
  b8 <- construct_basis("age (13,28,44,60,75) +i")
  b9 <- construct_basis("age (13,20,30,40,60)")
  b10 <- construct_basis("year (10,13,16,19,22)", window_start=window_start,
                         window_end=window_end)
  b11 <- construct_basis("year (10,13,16,19,22) +i", window_start=window_start,
                         window_end=window_end)
  b12 <- construct_basis("year (17,...,22)", window_start=window_start,
                         window_end=window_end)
  
  # Apply spline bases to dataframe
  if (model_version %in% c(0:23)) {
    d$dat_i$b2_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,1) }),4)
    d$dat_i$b2_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,2) }),4)
    d$dat_i$b2_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,3) }),4)
    d$dat_i$b2_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,4) }),4)
    d$dat_i$b3_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,1) }),4)
    d$dat_i$b3_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,2) }),4)
    d$dat_i$b3_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,3) }),4)
    d$dat_i$b3_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,4) }),4)
    d$dat_i$b4_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,1) }),4)
    d$dat_i$b4_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,2) }),4)
    d$dat_i$b4_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,3) }),4)
    d$dat_i$b4_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,4) }),4)
  } else if (model_version %in% c(24:25)) {
    d$dat_i$b5_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b5(j,1) }),4)
    d$dat_i$b5_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b5(j,2) }),4)
    d$dat_i$b5_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b5(j,3) }),4)
    d$dat_i$b5_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b5(j,4) }),4)
    d$dat_i$b6_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b6(w_2,1) }),4)
    d$dat_i$b6_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b6(w_2,2) }),4)
    d$dat_i$b6_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b6(w_2,3) }),4)
    d$dat_i$b6_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b6(w_2,4) }),4)
    d$dat_i$b7_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b7(j,1) }),4)
    d$dat_i$b7_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b7(j,2) }),4)
    d$dat_i$b7_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b7(j,3) }),4)
    d$dat_i$b7_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b7(j,4) }),4)
    d$dat_i$b7_5 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b7(j,5) }),4)
    d$dat_i$b8_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b8(w_2,1) }),4)
    d$dat_i$b8_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b8(w_2,2) }),4)
    d$dat_i$b8_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b8(w_2,3) }),4)
    d$dat_i$b8_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b8(w_2,4) }),4)
    d$dat_i$b8_5 <- signif(sapply(d$dat_i$w_2, function(w_2) { b8(w_2,5) }),4)
    d$dat_i$b9_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,1) }),4)
    d$dat_i$b9_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,2) }),4)
    d$dat_i$b9_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,3) }),4)
    d$dat_i$b9_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,4) }),4)
  } else if (model_version %in% c(26:27)) {
    d$dat_i$b9_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,1) }),4)
    d$dat_i$b9_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,2) }),4)
    d$dat_i$b9_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,3) }),4)
    d$dat_i$b9_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,4) }),4)
    d$dat_i$b10_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,1) }),4)
    d$dat_i$b10_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,2) }),4)
    d$dat_i$b10_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,3) }),4)
    d$dat_i$b10_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,4) }),4)
    d$dat_i$b11_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b11(j,1) }),4)
    d$dat_i$b11_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b11(j,2) }),4)
    d$dat_i$b11_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b11(j,3) }),4)
    d$dat_i$b11_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b11(j,4) }),4)
    d$dat_i$b11_5 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b11(j,5) }),4)
  } else if (model_version==28) {
    d$dat_i$b9_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,1) }),4)
    d$dat_i$b9_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,2) }),4)
    d$dat_i$b9_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,3) }),4)
    d$dat_i$b9_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b9(w_2,4) }),4)
    d$dat_i$b10_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,1) }),4)
    d$dat_i$b10_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,2) }),4)
    d$dat_i$b10_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,3) }),4)
    d$dat_i$b10_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,4) }),4)
  } else if (model_version %in% c(29:33)) {
    d$dat_i$b9_1 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,1) }),4)
    d$dat_i$b9_2 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,2) }),4)
    d$dat_i$b9_3 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,3) }),4)
    d$dat_i$b9_4 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,4) }),4)
    d$dat_i$b10_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,1) }),4)
    d$dat_i$b10_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,2) }),4)
    d$dat_i$b10_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,3) }),4)
    d$dat_i$b10_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,4) }),4)
  } else if (model_version==34) {
    d$dat_i$b9_1 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,1) }),4)
    d$dat_i$b9_2 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,2) }),4)
    d$dat_i$b9_3 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,3) }),4)
    d$dat_i$b9_4 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,4) }),4)
    d$dat_i$b12_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,1) }),4)
    d$dat_i$b12_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,2) }),4)
    d$dat_i$b12_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,3) }),4)
    d$dat_i$b12_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,4) }),4)
  }
  
  #' #' Calculate likelihood for individual i
  #' #'
  #' #' @param i Index for an individual
  #' #' @param params TO DO
  #' #' @param inds TO DO
  #' #' @return Numeric likelihood
  #' #' @note dat_objs is accessed globally
  #' lik_fn <- function(i, params, inds) {
  #'   
  #'   print(paste0("i: ", i))
  #'   # Extract data for individual i
  #'   d <- dat_objs[[i]]
  #'   
  #'   # Compute the likelihood for individual i
  #'   f2 <- sum(unlist(lapply(d$X_i_set, function(x) {
  #'     d$dat_i$x <- x
  #'     if (length(x)==1) {
  #'       d$dat_i$x_prev <- 0
  #'     } else {
  #'       d$dat_i$x_prev <- c(0,x[1:(length(x)-1)]) # !!!!! This can be precomputed
  #'     }
  #'     return(prod(apply(X=d$dat_i, MARGIN=1, FUN = function(r) {
  #'       x <- r[["x"]]
  #'       w_ij <- as.numeric(r[inds$w])
  #'       j <- r[["cal_time_sc"]]
  #'       spl_ij <- r[inds$spl]
  #'       return(
  #'         f_x(x=x, x_prev=r[["x_prev"]], w=w_ij, j=j,
  #'             s=r[["init_visit"]], spl=spl_ij, params=params) *
  #'           f_y(y=r[["y"]], x=x, w=w_ij, z=r[["z"]], j=j,
  #'               spl=spl_ij, params=params)
  #'       )
  #'     })))
  #'   })))
  #'   
  #'   browser() # !!!!!
  #' 
  #'   if (f2<=0) {
  #'     f2 <- 1e-10
  #'     # warning("Likelihood of zero")
  #'   }
  #'   
  #'   return(f2)
  #'   
  #' }
  
  #' Calculate likelihood component f_x
  if (cfg$model_version %in% c(0:23)) {
    
    if (cfg$model_version==0) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(params$a_x + params$t_x*j + sum(params$g_x*w))
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==1) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(params$a_x)
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s)
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==2) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(params$a_x + params$g_x1*w[1])
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + params$g_s1*w[1])
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version %in% c(3,4,5)) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(params$a_x + sum(params$g_x*w))
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + sum(params$g_s*w))
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==6) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(params$a_x + params$t_x*j + sum(params$g_x*w))
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + sum(params$g_s*w))
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version %in% c(7:9)) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(params$a_x + params$t_x*j + sum(params$g_x*w))
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==10) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) { 1 }
      
    } else if (cfg$model_version==11) {
      
      b1 <- construct_basis("age (0-100), 4DF")
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$a_x + params$t_x*j + params$g_x[1]*w[1] +
                params$g_x[2]*b1(w[2],1) + params$g_x[3]*b1(w[2],2) +
                params$g_x[4]*b1(w[2],3) + params$g_x[5]*b1(w[2],4)
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==12) {
      
      f_x <- (function() {
        b2 <- construct_basis("age (13,20,30,60,90)")
        f_x <- function(x, x_prev, w, j, s, spl, params) {
          if (s==0) {
            if (x==1 && x_prev==1) {
              return(1)
            } else {
              prob <- icll(
                params$a_x + params$t_x*j + params$g_x[1]*w[1] +
                  params$g_x[2]*b2(w[2],1) + params$g_x[3]*b2(w[2],2) +
                  params$g_x[4]*b2(w[2],3) + params$g_x[5]*b2(w[2],4)
              )
              if (x==1) { return(prob) } else { return(1-prob) }
            }
          } else {
            prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        }
        return(f_x)
      })()
      
    } else if (cfg$model_version %in% c(13,14)) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$a_x + params$t_x[1]*spl[["b4_1"]] +
                params$t_x[2]*spl[["b4_2"]] +
                params$t_x[3]*spl[["b4_3"]] + params$t_x[4]*spl[["b4_4"]] +
                params$g_x[1]*w[1] + params$g_x[2]*spl[["b2_1"]] +
                params$g_x[3]*spl[["b2_2"]] + params$g_x[4]*spl[["b2_3"]] +
                params$g_x[5]*spl[["b2_4"]]
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==15) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$a_x + params$t_x[1]*spl[["b4_1"]] +
                params$t_x[2]*spl[["b4_2"]] + params$t_x[3]*spl[["b4_3"]] +
                params$t_x[4]*spl[["b4_4"]] +
                w[1]*(
                  params$g_x[1]*spl[["b2_1"]] + params$g_x[2]*spl[["b2_2"]] +
                    params$g_x[3]*spl[["b2_3"]] + params$g_x[4]*spl[["b2_4"]]
                ) +
                (1-w[1])*(
                  params$g_x[5]*spl[["b2_1"]] + params$g_x[6]*spl[["b2_2"]] +
                    params$g_x[7]*spl[["b2_3"]] + params$g_x[8]*spl[["b2_4"]]
                )
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==16) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$a_x + params$t_x[1]*spl[["b4_1"]] +
                params$t_x[2]*spl[["b4_2"]] + params$t_x[3]*spl[["b4_3"]] +
                params$t_x[4]*spl[["b4_4"]] +
                w[1]*(
                  params$g_x[1]*spl[["b2_1"]] + params$g_x[2]*spl[["b2_2"]] +
                    params$g_x[3]*spl[["b2_3"]] + params$g_x[4]*spl[["b2_4"]]
                ) +
                (1-w[1])*(
                  params$g_x[5]*spl[["b2_1"]] + params$g_x[6]*spl[["b2_2"]] +
                    params$g_x[7]*spl[["b2_3"]] + params$g_x[8]*spl[["b2_4"]]
                )
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(
            params$a_s + params$t_s[1]*spl[["b4_1"]] + params$t_s[2]*spl[["b4_2"]] +
              params$t_s[3]*spl[["b4_3"]] + params$t_s[4]*spl[["b4_4"]] +
              sum(params$g_s*w)
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version %in% c(17,18)) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$a_x + params$t_x[1]*spl[["b4_1"]] +
                params$t_x[2]*spl[["b4_2"]] + params$t_x[3]*spl[["b4_3"]] +
                params$t_x[4]*spl[["b4_4"]] +
                w[1]*(
                  params$g_x[1]*spl[["b2_1"]] + params$g_x[2]*spl[["b2_2"]] +
                    params$g_x[3]*spl[["b2_3"]] + params$g_x[4]*spl[["b2_4"]]
                ) +
                (1-w[1])*(
                  params$g_x[5]*spl[["b2_1"]] + params$g_x[6]*spl[["b2_2"]] +
                    params$g_x[7]*spl[["b2_3"]] + params$g_x[8]*spl[["b2_4"]]
                )
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(
            params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b3_1"]] +
              params$g_s[3]*spl[["b3_2"]] + params$g_s[4]*spl[["b3_3"]] + 
              params$g_s[5]*spl[["b3_4"]]
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version %in% c(19:21)) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$a_x + params$t_x[1]*spl[["b5_1"]] +
                params$t_x[2]*spl[["b5_2"]] + params$t_x[3]*spl[["b5_3"]] +
                params$t_x[4]*spl[["b5_4"]] +
                w[1]*(
                  params$g_x[1]*spl[["b2_1"]] + params$g_x[2]*spl[["b2_2"]] +
                    params$g_x[3]*spl[["b2_3"]] + params$g_x[4]*spl[["b2_4"]]
                ) +
                (1-w[1])*(
                  params$g_x[5]*spl[["b2_1"]] + params$g_x[6]*spl[["b2_2"]] +
                    params$g_x[7]*spl[["b2_3"]] + params$g_x[8]*spl[["b2_4"]]
                )
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(
            params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b3_1"]] +
              params$g_s[3]*spl[["b3_2"]] + params$g_s[4]*spl[["b3_3"]] + 
              params$g_s[5]*spl[["b3_4"]]
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==22) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$a_x + params$t_x[1]*spl[["b5_1"]] +
                params$t_x[2]*spl[["b5_2"]] + params$t_x[3]*spl[["b5_3"]] +
                params$t_x[4]*spl[["b5_4"]] +
                w[1]*(
                  params$g_x[1]*spl[["b6_1"]] + params$g_x[2]*spl[["b6_2"]] +
                    params$g_x[3]*spl[["b6_3"]] + params$g_x[4]*spl[["b6_4"]]
                ) +
                (1-w[1])*(
                  params$g_x[5]*spl[["b6_1"]] + params$g_x[6]*spl[["b6_2"]] +
                    params$g_x[7]*spl[["b6_3"]] + params$g_x[8]*spl[["b6_4"]]
                )
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(
            params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b6_1"]] +
              params$g_s[3]*spl[["b6_2"]] + params$g_s[4]*spl[["b6_3"]] + 
              params$g_s[5]*spl[["b6_4"]]
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    } else if (cfg$model_version==23) {
      
      f_x <- function(x, x_prev, w, j, s, spl, params) {
        if (s==0) {
          if (x==1 && x_prev==1) {
            return(1)
          } else {
            prob <- icll(
              params$t_x[1]*spl[["b5_1"]] + params$t_x[2]*spl[["b5_2"]] +
                params$t_x[3]*spl[["b5_3"]] + params$t_x[4]*spl[["b5_4"]] +
                w[1]*(
                  params$g_x[1]*spl[["b8_1"]] + params$g_x[2]*spl[["b8_2"]] +
                    params$g_x[3]*spl[["b8_3"]] + params$g_x[4]*spl[["b8_4"]] +
                    params$g_x[5]*spl[["b8_5"]]
                ) +
                (1-w[1])*(
                  params$g_x[6]*spl[["b8_1"]] + params$g_x[7]*spl[["b8_2"]] +
                    params$g_x[8]*spl[["b8_3"]] + params$g_x[9]*spl[["b8_4"]] +
                    params$g_x[10]*spl[["b8_5"]]
                )
            )
            if (x==1) { return(prob) } else { return(1-prob) }
          }
        } else {
          prob <- icll(
            params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b6_1"]] +
              params$g_s[3]*spl[["b6_2"]] + params$g_s[4]*spl[["b6_3"]] + 
              params$g_s[5]*spl[["b6_4"]]
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      }
      
    }
    
  }
  if (cfg$model_version==24) {
    
    f_x <- function(x, x_prev, w, j, s, spl, params) {
      if (s==0) {
        if (x==1 && x_prev==1) {
          return(1)
        } else {
          prob <- icll(
            params$a_x + params$t_x[1]*spl[["b5_1"]] +
              params$t_x[2]*spl[["b5_2"]] + params$t_x[3]*spl[["b5_3"]] +
              params$t_x[4]*spl[["b5_4"]] +
              w[1]*(
                params$g_x[1]*spl[["b6_1"]] + params$g_x[2]*spl[["b6_2"]] +
                  params$g_x[3]*spl[["b6_3"]] + params$g_x[4]*spl[["b6_4"]]
              ) +
              (1-w[1])*(
                params$g_x[5]*spl[["b6_1"]] + params$g_x[6]*spl[["b6_2"]] +
                  params$g_x[7]*spl[["b6_3"]] + params$g_x[8]*spl[["b6_4"]]
              )
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      } else {
        prob <- icll(
          params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b6_1"]] +
            params$g_s[3]*spl[["b6_2"]] + params$g_s[4]*spl[["b6_3"]] + 
            params$g_s[5]*spl[["b6_4"]]
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    }
    
  } else if (cfg$model_version==25) {
    
    f_x <- function(x, x_prev, w, j, s, spl, params) {
      if (s==0) {
        if (x==1 && x_prev==1) {
          return(1)
        } else {
          prob <- icll(
            params$a_x + params$t_x[1]*spl[["b5_1"]] +
              params$t_x[2]*spl[["b5_2"]] + params$t_x[3]*spl[["b5_3"]] +
              params$t_x[4]*spl[["b5_4"]] +
              w[1]*(
                params$g_x[1]*spl[["b9_1"]] + params$g_x[2]*spl[["b9_2"]] +
                  params$g_x[3]*spl[["b9_3"]] + params$g_x[4]*spl[["b9_4"]]
              ) +
              (1-w[1])*(
                params$g_x[5]*spl[["b9_1"]] + params$g_x[6]*spl[["b9_2"]] +
                  params$g_x[7]*spl[["b9_3"]] + params$g_x[8]*spl[["b9_4"]]
              )
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      } else {
        prob <- icll(
          params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b9_1"]] +
            params$g_s[3]*spl[["b9_2"]] + params$g_s[4]*spl[["b9_3"]] + 
            params$g_s[5]*spl[["b9_4"]]
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    }
    
  } else if (cfg$model_version %in% c(26:28)) {
    
    f_x <- function(x, x_prev, w, j, s, spl, params) {
      if (s==0) {
        if (x==1 && x_prev==1) {
          return(1)
        } else {
          prob <- icll(
            params$a_x + params$t_x[1]*spl[["b10_1"]] +
              params$t_x[2]*spl[["b10_2"]] + params$t_x[3]*spl[["b10_3"]] +
              params$t_x[4]*spl[["b10_4"]] +
              w[1]*(
                params$g_x[1]*spl[["b9_1"]] + params$g_x[2]*spl[["b9_2"]] +
                  params$g_x[3]*spl[["b9_3"]] + params$g_x[4]*spl[["b9_4"]]
              ) +
              (1-w[1])*(
                params$g_x[5]*spl[["b9_1"]] + params$g_x[6]*spl[["b9_2"]] +
                  params$g_x[7]*spl[["b9_3"]] + params$g_x[8]*spl[["b9_4"]]
              )
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      } else {
        prob <- icll(
          params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b9_1"]] +
            params$g_s[3]*spl[["b9_2"]] + params$g_s[4]*spl[["b9_3"]] + 
            params$g_s[5]*spl[["b9_4"]]
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    }
    
  }
  
  #' Calculate likelihood component f_y
  if (cfg$model_version %in% c(0:23)) {
    
    if (cfg$model_version==0) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$a_y + params$t_y*j + sum(params$g_y*w) + params$beta_x*x*(1-z) +
            params$beta_z*x*z
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==1) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(params$a_y + params$beta_x*x*(1-z) + params$beta_z*x*z)
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==2) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$a_y + params$g_y1*w[1] + params$beta_x*x*(1-z) + params$beta_z*x*z
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==3) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$a_y + sum(params$g_y*w) + params$beta_x*x*(1-z) + params$beta_z*x*z
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==4) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$a_y + sum(params$g_y*w) + params$beta_x*x*(1-z) + params$beta_z*x*z
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version %in% c(5:7)) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$a_y + params$t_y*j + sum(params$g_y*w) + params$beta_x*x
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==8) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$a_y + params$t_y*j + params$g_y[1]*w[1] + params$g_y[2]*w[2] +
            params$g_y[3]*(w[2]^2) + params$g_y[4]*(w[2]^3) +
            params$beta_x*x*(1-z) + params$beta_z*x*z
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version %in% c(9,11)) {
      
      b1 <- construct_basis("age (0-100), 4DF")
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$a_y + params$t_y*j + params$g_y[1]*w[1] +
            params$g_y[2]*b1(w[2],1) + params$g_y[3]*b1(w[2],2) +
            params$g_y[4]*b1(w[2],3) + params$g_y[5]*b1(w[2],4) +
            params$beta_x*x*(1-z) + params$beta_z*x*z
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==10) {
      
      # Note: there was a bug in this model when it was first run (beta_x excluded)
      b1 <- construct_basis("age (0-100), 4DF")
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$beta_x*x*(1-z) + params$beta_z*x*z + params$a_y +
            params$g_y[1]*w[1] + params$g_y[2]*b1(w[2],1) +
            params$g_y[3]*b1(w[2],2) + params$g_y[4]*b1(w[2],3) +
            params$g_y[5]*b1(w[2],4) + params$t_y*j
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version %in% c(12,13)) {
      
      # Note: there was a bug in this model when it was first run (beta_x excluded)
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$beta_x*x*(1-z) + params$beta_z*x*z + params$a_y +
            params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
            params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
            params$g_y[5]*spl[["b3_4"]] + params$t_y*j
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version %in% c(14:17)) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          params$beta_x*x*(1-z) + params$beta_z*x*z + params$a_y +
            params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
            params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
            params$g_y[5]*spl[["b3_4"]] + params$t_y[1]*spl[["b4_1"]] +
            params$t_y[2]*spl[["b4_2"]] + params$t_y[3]*spl[["b4_3"]] +
            params$t_y[4]*spl[["b4_4"]]
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==18) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          (params$beta_x[1]+params$beta_x[2]*j)*x*(1-z) +
            (params$beta_z[1]+params$beta_z[2]*j)*x*z + params$a_y +
            params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
            params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
            params$g_y[5]*spl[["b3_4"]] + params$t_y[1]*spl[["b4_1"]] +
            params$t_y[2]*spl[["b4_2"]] + params$t_y[3]*spl[["b4_3"]] +
            params$t_y[4]*spl[["b4_4"]]
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==19) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          (params$beta_x[1]+params$beta_x[2]*j)*x*(1-z) +
            (params$beta_z[1]+params$beta_z[2]*j)*x*z + params$a_y +
            params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
            params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
            params$g_y[5]*spl[["b3_4"]] + params$t_y[1]*spl[["b5_1"]] +
            params$t_y[2]*spl[["b5_2"]] + params$t_y[3]*spl[["b5_3"]] +
            params$t_y[4]*spl[["b5_4"]]
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==20) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          x*(1-z)*(
            params$beta_x[1]*spl[["b5_1"]] + params$beta_x[2]*spl[["b5_2"]] +
              params$beta_x[3]*spl[["b5_3"]] + params$beta_x[4]*spl[["b5_4"]]
          ) + x*z*(
            params$beta_z[1]*spl[["b5_1"]] + params$beta_z[2]*spl[["b5_2"]] +
              params$beta_z[3]*spl[["b5_3"]] + params$beta_z[4]*spl[["b5_4"]]
          ) +
            params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
            params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
            params$g_y[5]*spl[["b3_4"]] + params$t_y[1]*spl[["b5_1"]] +
            params$t_y[2]*spl[["b5_2"]] + params$t_y[3]*spl[["b5_3"]] +
            params$t_y[4]*spl[["b5_4"]]
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==21) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          x*(
            params$beta_x[1]*spl[["b5_1"]] + params$beta_x[2]*spl[["b5_2"]] +
              params$beta_x[3]*spl[["b5_3"]] + params$beta_x[4]*spl[["b5_4"]]
          ) +
            params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
            params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
            params$g_y[5]*spl[["b3_4"]] + params$t_y[1]*spl[["b5_1"]] +
            params$t_y[2]*spl[["b5_2"]] + params$t_y[3]*spl[["b5_3"]] +
            params$t_y[4]*spl[["b5_4"]]
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==22) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          x*(
            params$beta_x[1]*spl[["b5_1"]] + params$beta_x[2]*spl[["b5_2"]] +
              params$beta_x[3]*spl[["b5_3"]] + params$beta_x[4]*spl[["b5_4"]]
          ) +
            params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b6_1"]] +
            params$g_y[3]*spl[["b6_2"]] + params$g_y[4]*spl[["b6_3"]] +
            params$g_y[5]*spl[["b6_4"]] + params$t_y[1]*spl[["b5_1"]] +
            params$t_y[2]*spl[["b5_2"]] + params$t_y[3]*spl[["b5_3"]] +
            params$t_y[4]*spl[["b5_4"]]
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    } else if (cfg$model_version==23) {
      
      f_y <- function(y, x, w, z, j, spl, params) {
        prob <- icll(
          x*(
            params$beta_x[1]*spl[["b7_1"]] + params$beta_x[2]*spl[["b7_2"]] +
              params$beta_x[3]*spl[["b7_3"]] + params$beta_x[4]*spl[["b7_4"]] +
              params$beta_x[5]*spl[["b7_5"]]
          ) +
            params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b6_1"]] +
            params$g_y[3]*spl[["b6_2"]] + params$g_y[4]*spl[["b6_3"]] +
            params$g_y[5]*spl[["b6_4"]] + params$t_y[1]*spl[["b5_1"]] +
            params$t_y[2]*spl[["b5_2"]] + params$t_y[3]*spl[["b5_3"]] +
            params$t_y[4]*spl[["b5_4"]]
        )
        if (y==1) { return(prob) } else { return(1-prob) }
      }
      
    }
    
  }
  if (cfg$model_version==24) {
    
    f_y <- function(y, x, w, z, j, spl, params) {
      prob <- icll(
        x*(
          params$beta_x[1]*spl[["b7_1"]] + params$beta_x[2]*spl[["b7_2"]] +
            params$beta_x[3]*spl[["b7_3"]] + params$beta_x[4]*spl[["b7_4"]] +
            params$beta_x[5]*spl[["b7_5"]]
        ) +
          params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b6_1"]] +
          params$g_y[3]*spl[["b6_2"]] + params$g_y[4]*spl[["b6_3"]] +
          params$g_y[5]*spl[["b6_4"]] + params$t_y[1]*spl[["b5_1"]] +
          params$t_y[2]*spl[["b5_2"]] + params$t_y[3]*spl[["b5_3"]] +
          params$t_y[4]*spl[["b5_4"]]
      )
      if (y==1) { return(prob) } else { return(1-prob) }
    }
    
  } else if (cfg$model_version==25) {
    
    f_y <- function(y, x, w, z, j, spl, params) {
      prob <- icll(
        x*(
          params$beta_x[1]*spl[["b7_1"]] + params$beta_x[2]*spl[["b7_2"]] +
            params$beta_x[3]*spl[["b7_3"]] + params$beta_x[4]*spl[["b7_4"]] +
            params$beta_x[5]*spl[["b7_5"]]
        ) +
          params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b9_1"]] +
          params$g_y[3]*spl[["b9_2"]] + params$g_y[4]*spl[["b9_3"]] +
          params$g_y[5]*spl[["b9_4"]] + params$t_y[1]*spl[["b5_1"]] +
          params$t_y[2]*spl[["b5_2"]] + params$t_y[3]*spl[["b5_3"]] +
          params$t_y[4]*spl[["b5_4"]]
      )
      if (y==1) { return(prob) } else { return(1-prob) }
    }
    
  } else if (cfg$model_version==26) {
    
    f_y <- function(y, x, w, z, j, spl, params) {
      prob <- icll(
        x*(
          params$beta_x[1]*spl[["b11_1"]] + params$beta_x[2]*spl[["b11_2"]] +
            params$beta_x[3]*spl[["b11_3"]] + params$beta_x[4]*spl[["b11_4"]] +
            params$beta_x[5]*spl[["b11_5"]]
        ) +
          params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b9_1"]] +
          params$g_y[3]*spl[["b9_2"]] + params$g_y[4]*spl[["b9_3"]] +
          params$g_y[5]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
          params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
          params$t_y[4]*spl[["b10_4"]]
      )
      if (y==1) { return(prob) } else { return(1-prob) }
    }
    
  } else if (cfg$model_version==27) {
    
    f_y <- function(y, x, w, z, j, spl, params) {
      prob <- icll(
        params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b9_1"]] +
          params$g_y[3]*spl[["b9_2"]] + params$g_y[4]*spl[["b9_3"]] +
          params$g_y[5]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
          params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
          params$t_y[4]*spl[["b10_4"]]
      )
      if (y==1) { return(prob) } else { return(1-prob) }
    }
    
  } else if (cfg$model_version==28) {
    
    f_y <- function(y, x, w, z, j, spl, params) {
      prob <- icll(
        x*(
          (params$beta_x[1] + params$beta_x[2]*j) +
            (params$beta_x[3] + params$beta_x[4]*j) * max(w[2]-0.3,0) +
            (params$beta_x[5] + params$beta_x[6]*j) * max(w[2]-0.45,0)
        ) +
          params$a_y + params$g_y[1]*w[1] + params$g_y[2]*spl[["b9_1"]] +
          params$g_y[3]*spl[["b9_2"]] + params$g_y[4]*spl[["b9_3"]] +
          params$g_y[5]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
          params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
          params$t_y[4]*spl[["b10_4"]]
      )
      if (y==1) { return(prob) } else { return(1-prob) }
    }
    
  }
  
  # Convert parameter vector to a named list
  if (model_version==0) {
    params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], t_x=p[9], t_y=p[10], a_s=p[11], t_s=p[12], g_s=p[13:14])
  } else if (model_version==1) {
    params <- list(a_x=p[1], a_y=p[2], beta_x=p[3], beta_z=p[4], a_s=p[5])
  } else if (model_version==2) {
    params <- list(a_x=p[1], a_y=p[2], beta_x=p[3], beta_z=p[4], a_s=p[5], g_x1=p[6], g_y1=p[7], g_s1=p[8])
  } else if (model_version %in% c(3,4)) {
    params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], a_s=p[9], g_s=p[10:11])
  } else if (model_version==5) {
    params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], t_y=p[9], a_s=p[10], g_s=p[11:12])
  } else if (model_version==6) {
    params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], t_x=p[9], t_y=p[10], a_s=p[11], g_s=p[12:13])
  } else if (model_version==7) {
    params <- list(a_x=p[1], g_x=p[2:3], t_x=p[4], a_s=p[5], g_s=p[6:7], t_s=p[8], beta_x=p[9], a_y=p[10], g_y=p[11:12], t_y=p[13])
  } else if (model_version==8) {
    params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:8], beta_x=p[9], beta_z=p[10], t_x=p[11], t_y=p[12], a_s=p[13], t_s=p[14], g_s=p[15:16])
  } else if (model_version==9) {
    params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:9], beta_x=p[10], beta_z=p[11], t_x=p[12], t_y=p[13], a_s=p[14], t_s=p[15], g_s=p[16:17])
  } else if (model_version==10) {
    params <- list(beta_z=p[1], a_y=p[2], g_y=p[3:7], t_y=p[8])
  } else if (model_version %in% c(11,12)) {
    params <- list(a_x=p[1], g_x=p[2:6], t_x=p[7], a_s=p[8], g_s=p[9:10], t_s=p[11], beta_x=p[12], beta_z=p[13], a_y=p[14], g_y=p[15:19], t_y=p[20])
  } else if (model_version==13) {
    params <- list(a_x=p[1], g_x=p[2:6], t_x=p[7:10], a_s=p[11], g_s=p[12:13], t_s=p[14], beta_x=p[15], beta_z=p[16], a_y=p[17], g_y=p[18:22], t_y=p[23])
  } else if (model_version==14) {
    params <- list(a_x=p[1], g_x=p[2:6], t_x=p[7:10], a_s=p[11], g_s=p[12:13], t_s=p[14], beta_x=p[15], beta_z=p[16], a_y=p[17], g_y=p[18:22], t_y=p[23:26])
  } else if (model_version==15) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:16], t_s=p[17], beta_x=p[18], beta_z=p[19], a_y=p[20], g_y=p[21:25], t_y=p[26:29])
  } else if (model_version==16) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:16], t_s=p[17:20], beta_x=p[21], beta_z=p[22], a_y=p[23], g_y=p[24:28], t_y=p[29:32])
  } else if (model_version==17) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20], beta_z=p[21], a_y=p[22], g_y=p[23:27], t_y=p[28:31])
  } else if (model_version %in% c(18,19)) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20:21], beta_z=p[22:23], a_y=p[24], g_y=p[25:29], t_y=p[30:33])
  } else if (model_version==20) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20:23], beta_z=p[24:27], a_y=p[28], g_y=p[29:33], t_y=p[34:37])
  } else if (model_version %in% c(21:22)) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20:23], a_y=p[24], g_y=p[25:29], t_y=p[30:33])
  } else if (model_version==23) {
    params <- list(g_x=p[1:10], t_x=p[11:14], a_s=p[15], g_s=p[16:20], beta_x=p[21:25], a_y=p[26], g_y=p[27:31], t_y=p[32:35])
  } else if (model_version %in% c(24:26)) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20:24], a_y=p[25], g_y=p[26:30], t_y=p[31:34])
  } else if (model_version==27) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], a_y=p[20], g_y=p[21:25], t_y=p[26:29])
  } else if (model_version==28) {
    params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20:25], a_y=p[26], g_y=p[27:31], t_y=p[32:35])
  }
  
  # Set initial parameter estimates
  if (cfg$model_version==0) {
    `par_init` <- c(a_x=-5.603, g_x1=0, g_x2=-0.3655, a_y=-6.020, g_y1=0, g_y2=4.282, beta_x=1.401, beta_z=0.0004, t_x=0.9609, t_y=-3.906, a_s=-1.740, t_s=-2.100, g_s1=0, g_s2=1.271) # Model iteration 0
  } else if (cfg$model_version==1) {
    par_init <- c(a_x=-5.651, a_y=-4.942, beta_x=1.423, beta_z=0.2235, a_s=-2.007)
  } else if (cfg$model_version==2) {
    par_init <- c(a_x=-5.651, a_y=-4.942, beta_x=1.423, beta_z=0.2235, a_s=-2.007, g_x1=0, g_y1=0, g_s1=0)
  } else if (cfg$model_version %in% c(3,4)) {
    par_init <- c(a_x=-5.8039, g_x1=-0.5518, g_x2=0.6733, a_y=-6.4375, g_y1=0.3011, g_y2=4.0686, beta_x=1.6762, beta_z=0.8045, a_s=-2.116, g_s1=-0.3937, g_s2=1.1987)
  } else if (cfg$model_version==5) {
    par_init <- c(a_x=-5.8039, g_x1=-0.5518, g_x2=0.6733, a_y=-6.4375, g_y1=0.3011, g_y2=4.0686, beta_x=1.6762, beta_z=0.8045, t_y=0, a_s=-2.116, g_s1=-0.3937, g_s2=1.1987)
  } else if (cfg$model_version==6) {
    par_init <- c(a_x=-4.6185, g_x1=-1.3117, g_x2=-0.3883, a_y=-6.1581, g_y1=0.3794, g_y2=4.4762, beta_x=1.8497, beta_z=1.708, t_x=0, t_y=-0.45141, a_s=-2.495, g_s1=-0.0177, g_s2=1.684)
  } else if (cfg$model_version==7) {
    par_init <- c(a_x=-6.2967, g_x1=-0.1535, g_x2=0.9796, t_x=0.5343, a_s=-2.3111, g_s1=-0.5649, g_s2=0.6198, t_s=0.4245, beta_x=1.401, a_y=-5.5786, g_y1=0.3278, g_y2=4.2046, t_y=-0.7198)
  } else if (cfg$model_version==8) {
    par_init <- c(a_x=-3.5607, g_x1=-0.3244, g_x2=-0.2809, a_y=-5.7446, g_y1=0.3544, g_y2=4.4057, g_y3=0, g_y4=0, beta_x=1.8096, beta_z=1.8153, t_x=-0.786, t_y=-0.7826, a_s=-2.87, t_s=0.6349, g_s1=-0.3768, g_s2=0.6409)
  } else if (cfg$model_version==9) {
    par_init <- c(a_x=-2.2308, g_x1=-0.4977, g_x2=-0.9101, a_y=-6.3404, g_y1=0.5996, g_y2=2.4098, g_y3=2.8139, g_y4=7.0955, g_y5=6.0127, beta_x=1.295, beta_z=1.2856, t_x=-1.366, t_y=-0.7141, a_s=-2.0978, t_s=0.3321, g_s1=-0.8771, g_s2=0.8316)
  } else if (cfg$model_version==10) {
    par_init <- c(beta_z=0.3, a_y=-9.4955, g_y1=0.3209, g_y2=5.7549, g_y3=5.2759, g_y4=13.7284, g_y5=5.2979, t_y=-0.6637)
  } else if (cfg$model_version==11) {
    par_init <- c(a_x=-2.2308, g_x1=-0.4977, g_x2=0, g_x3=0, g_x4=0, g_x5=0, t_x=-1.366, a_s=-2.0978, g_s1=-0.8771, g_s2=0.8316, t_s=0.3321, beta_x=1.295, beta_z=1.2856, a_y=-6.3404, g_y1=0.5996, g_y2=2.4098, g_y3=2.8139, g_y4=7.0955, g_y5=6.0127, t_y=-0.7141)
  } else if (cfg$model_version==12) {
    par_init <- c(a_x=-4.667, g_x1=-0.7444, g_x2=3.658, g_x3=-1.970, g_x4=-0.8989, g_x5=-6.690, t_x=-1.169, a_s=-3.299, g_s1=-0.6594, g_s2=0.8443, t_s=1.051, beta_x=0.4, beta_z=0.4, a_y=-5.654, g_y1=0.2733, g_y2=1.759, g_y3=2.802, g_y4=6.180, g_y5=2.649, t_y=-0.6332)
  } else if (cfg$model_version==13) {
    par_init <- c(a_x=-6.535, g_x1=-0.6737, g_x2=3.636, g_x3=0.2734, g_x4=0.4366, g_x5=-8.512, t_x1=-1.578, t_x2=-0.2818, t_x3=0.3750, t_x4=-2.160, a_s=-3.369, g_s1=-0.5761, g_s2=0.8899, t_s=1.051, beta_x=1, beta_z=0.5072, a_y=-5.944, g_y1=0.3940, g_y2=1.871, g_y3=2.923, g_y4=6.809, g_y5=3.004, t_y=-0.6077)
  } else if (cfg$model_version==14) {
    par_init <- c(a_x=-6.3567, g_x1=-0.7201, g_x2=3.5877, g_x3=0.8982, g_x4=1.2344, g_x5=-7.8761, t_x1=-2.1962, t_x2=-0.396, t_x3=-0.0074, t_x4=-2.2542, a_s=-3.0662, g_s1=-0.6222, g_s2=0.8341, t_s=0.8966, beta_x=0.9409, beta_z=0.7236, a_y=-6.5048, g_y1=0.4056, g_y2=1.9911, g_y3=3.0964, g_y4=6.9362, g_y5=3.4698, t_y1=-0.2982, t_y2=-0.8746, t_y3=-0.5772, t_y4=-1.0723)
  } else if (cfg$model_version==15) {
    par_init <- c(a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66, g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593, t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0.753, t_s=0.936, beta_x=1.12, beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18, g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12)
  } else if (cfg$model_version==16) {
    par_init <- c(a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66, g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593, t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0.753, t_s1=0, t_s2=0, t_s3=0, t_s4=0, beta_x=1.12, beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18, g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12)
  } else if (cfg$model_version==17) {
    par_init <- c(a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66, g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593, t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0, g_s3=0, g_s4=0, g_s5=0, beta_x=1.12, beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18, g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12)
  } else if (cfg$model_version==18) {
    par_init <- c(a_x=-7.855, g_x1=4.791, g_x2=-0.935, g_x3=0.844, g_x4=-8.116, g_x5=2.804, g_x6=-2.155, g_x7=5.066, g_x8=-5.616, t_x1=-1.084, t_x2=-0.764, t_x3=-1.505, t_x4=-2.598, a_s=-3.043, g_s1=-0.281, g_s2=1.184, g_s3=0.219, g_s4=3.380, g_s5=-3.239, beta_x1=2.636, beta_x2=-1.051, beta_z1=3.983, beta_z2=-1.870, a_y=-7.221, g_y1=0.411, g_y2=1.692, g_y3=3.440, g_y4=6.288, g_y5=3.826, t_y1=0.413, t_y2=0.030, t_y3=0.530, t_y4=-0.235)
  } else if (cfg$model_version==19) {
    par_init <- c(a_x=-8.725, g_x1=3.414, g_x2=-0.561, g_x3=0.040, g_x4=-6.618, g_x5=-0.058, g_x6=0.498, g_x7=6.534, g_x8=-4.481, t_x1=-1.146, t_x2=-1.745, t_x3=-0.895, t_x4=-2.134, a_s=-2.914, g_s1=-0.435, g_s2=1.418, g_s3=-0.779, g_s4=3.834, g_s5=-2.681, beta_x1=1.813, beta_x2=-1.234, beta_z1=3.826, beta_z2=-3.966, a_y=-7.504, g_y1=0.598, g_y2=2.118, g_y3=3.794, g_y4=6.900, g_y5=4.007, t_y1=0.303, t_y2=-0.594, t_y3=-1.262, t_y4=-1.038)
  } else if (cfg$model_version==20) {
    par_init <- c(a_x=-8.182, g_x1=3.290, g_x2=-1.767, g_x3=2.711, g_x4=-2.527, g_x5=-0.778, g_x6=0.247, g_x7=7.216, g_x8=-3.587, t_x1=-0.226, t_x2=-1.011, t_x3=-1.618, t_x4=-1.725, a_s=-3.047, g_s1=-0.430, g_s2=1.355, g_s3=-0.465, g_s4=3.651, g_s5=-3.644, beta_x1=0.693, beta_x2=-0.017, beta_x3=2.102, beta_x4=-1.032, beta_z1=1.920, beta_z2=-0.819, beta_z3=1.866, beta_z4=-5.010, a_y=-7.120, g_y1=0.575, g_y2=2.194, g_y3=3.795, g_y4=7.084, g_y5=3.862, t_y1=-0.804, t_y2=-0.073, t_y3=-0.903, t_y4=-0.350)
  } else if (cfg$model_version==21) {
    par_init <- c(a_x=-8.106, g_x1=2.437, g_x2=-1.533, g_x3=2.540, g_x4=-2.320, g_x5=-1.679, g_x6=1.008, g_x7=7.416, g_x8=-2.923, t_x1=-0.284, t_x2=-1.655, t_x3=-0.782, t_x4=-1.090, a_s=-3.130, g_s1=-0.443, g_s2=1.369, g_s3=-0.218, g_s4=3.286, g_s5=-4.564, beta_x1=0.119, beta_x2=-0.316, beta_x3=2.065, beta_x4=-1.207, a_y=-6.998, g_y1=0.559, g_y2=2.351, g_y3=3.749, g_y4=7.466, g_y5=3.935, t_y1=-0.824, t_y2=0.082, t_y3=-1.909, t_y4=-0.774)
  } else if (cfg$model_version==22) {
    par_init <- c(a_x=-8.106, g_x1=2.437, g_x2=-1.533, g_x3=2.540, g_x4=-2.320, g_x5=-1.679, g_x6=1.008, g_x7=7.416, g_x8=-2.923, t_x1=-0.284, t_x2=-1.655, t_x3=-0.782, t_x4=-1.090, a_s=-3.130, g_s1=-0.443, g_s2=1.369, g_s3=-0.218, g_s4=3.286, g_s5=-4.564, beta_x1=0.119, beta_x2=-0.316, beta_x3=2.065, beta_x4=-1.207, a_y=-6.998, g_y1=0.559, g_y2=2.351, g_y3=3.749, g_y4=7.466, g_y5=3.935, t_y1=-0.824, t_y2=0.082, t_y3=-1.909, t_y4=-0.774)
  } else if (cfg$model_version==23) {
    par_init <- c(g_x1=1.378, g_x2=0.652, g_x3=-0.220, g_x4=-1.242, g_x5=3.564, g_x6=2.113, g_x7=-2.053, g_x8=-0.294, g_x9=-0.560, g_x10=-2.843, t_x1=-2.248, t_x2=-5.615, t_x3=-15.796, t_x4=0.996, a_s=-3.074, g_s1=-0.461, g_s2=3.276, g_s3=0.252, g_s4=4.272, g_s5=-0.583, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, beta_x5=0, a_y=-8.222, g_y1=0.657, g_y2=1.982, g_y3=3.492, g_y4=8.479, g_y5=3.900, t_y1=-0.883, t_y2=-0.201, t_y3=-1.170, t_y4=-1.428)
  } else if (cfg$model_version==24) {
    par_init <- c(a_x=-5.974, g_x1=-0.153, g_x2=-0.892, g_x3=1.491, g_x4=1.431, g_x5=-6.335, g_x6=-0.520, g_x7=4.633, g_x8=-1.664, t_x1=0.254, t_x2=-4.697, t_x3=-4.189, t_x4=0.077, a_s=-2.931, g_s1=-0.461, g_s2=3.276, g_s3=0.252, g_s4=4.272, g_s5=-0.583, beta_x1=1.257, beta_x2=0.545, beta_x3=-0.298, beta_x4=2.249, beta_x5=-1.334, a_y=-8.222, g_y1=0.657, g_y2=1.982, g_y3=3.492, g_y4=8.479, g_y5=3.900, t_y1=-0.883, t_y2=-0.201, t_y3=-1.170, t_y4=-1.428)
  } else if (cfg$model_version==25) {
    par_init <- c(a_x=-6.660, g_x1=2.739, g_x2=-1.388, g_x3=-0.418, g_x4=-2.872, g_x5=1.115, g_x6=-1.687, g_x7=3.922, g_x8=-3.564, t_x1=-0.108, t_x2=-2.928, t_x3=-1.238, t_x4=-0.435, a_s=-3.293, g_s1=-0.518, g_s2=3.702, g_s3=2.628, g_s4=4.248, g_s5=0.814, beta_x1=1.459, beta_x2=0.955, beta_x3=-0.187, beta_x4=2.743, beta_x5=-1.324, a_y=-8.834, g_y1=0.623, g_y2=2.767, g_y3=2.160, g_y4=5.665, g_y5=2.648, t_y1=-0.193, t_y2=0.521, t_y3=-0.161, t_y4=-0.711)
  } else if (cfg$model_version==26) {
    par_init <- c(a_x=-6.586, g_x1=1.622, g_x2=-0.340, g_x3=-1.187, g_x4=-2.215, g_x5=0.786, g_x6=-0.299, g_x7=3.819, g_x8=-2.592, t_x1=0.921, t_x2=-2.219, t_x3=-0.694, t_x4=0.235, a_s=-3.250, g_s1=-0.465, g_s2=3.692, g_s3=2.392, g_s4=4.333, g_s5=0.986, beta_x1=1.030, beta_x2=1.061, beta_x3=-0.070, beta_x4=2.966, beta_x5=-0.965, a_y=-8.645, g_y1=0.595, g_y2=2.588, g_y3=1.931, g_y4=5.047, g_y5=2.706, t_y1=-0.142, t_y2=0.554, t_y3=0.172, t_y4=-0.695)
  } else if (cfg$model_version==27) {
    par_init <- c(a_x=-6.586, g_x1=1.622, g_x2=-0.340, g_x3=-1.187, g_x4=-2.215, g_x5=0.786, g_x6=-0.299, g_x7=3.819, g_x8=-2.592, t_x1=0.921, t_x2=-2.219, t_x3=-0.694, t_x4=0.235, a_s=-3.250, g_s1=-0.465, g_s2=3.692, g_s3=2.392, g_s4=4.333, g_s5=0.986, a_y=-8.645, g_y1=0.595, g_y2=2.588, g_y3=1.931, g_y4=5.047, g_y5=2.706, t_y1=-0.142, t_y2=0.554, t_y3=0.172, t_y4=-0.695)
  } else if (cfg$model_version==28) {
    par_init <- c(a_x=-6.586, g_x1=1.622, g_x2=-0.340, g_x3=-1.187, g_x4=-2.215, g_x5=0.786, g_x6=-0.299, g_x7=3.819, g_x8=-2.592, t_x1=0.921, t_x2=-2.219, t_x3=-0.694, t_x4=0.235, a_s=-3.250, g_s1=-0.465, g_s2=3.692, g_s3=2.392, g_s4=4.333, g_s5=0.986, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, beta_x5=0, beta_x6=0, a_y=-8.645, g_y1=0.595, g_y2=2.588, g_y3=1.931, g_y4=5.047, g_y5=2.706, t_y1=-0.142, t_y2=0.554, t_y3=0.172, t_y4=-0.695)
  }
  
  # Construct spline bases (process.R)
  b1 <- construct_basis("age (0-100), 4DF")
  b2 <- construct_basis("age (13,20,30,60,90)")
  b3 <- construct_basis("age (13,30,60,75,90)")
  b4 <- construct_basis("year (00,05,10,15,20)", window_start=cfg2$w_start)
  b5 <- construct_basis("year (10,13,17,20,23)", window_start=cfg2$w_start)
  b6 <- construct_basis("age (13,28,44,60,75)")
  b7 <- construct_basis("year (10,13,17,20,23) +i", window_start=cfg2$w_start)
  b8 <- construct_basis("age (13,28,44,60,75) +i")
  b9 <- construct_basis("age (13,20,30,40,60)")
  b10 <- construct_basis("year (10,13,16,19,22)", window_start=cfg2$w_start,
                         window_end=cfg2$w_end)
  b11 <- construct_basis("year (10,13,16,19,22) +i", window_start=cfg2$w_start,
                         window_end=cfg2$w_end)
  b12 <- construct_basis("year (17,...,22)", window_start=cfg2$w_start,
                         window_end=cfg2$w_end)
  
  # Functions to return spline basis function as a matrix
  A_b5 <- function(j) {
    t(matrix(c(b5(j,1),b5(j,2),b5(j,3),b5(j,4))))
  }
  A_b6 <- function(w_2) {
    t(matrix(c(b6(w_2,1), b6(w_2,2), b6(w_2,3), b6(w_2,4))))
  }
  A_b7 <- function(j) {
    t(matrix(c(b7(j,1),b7(j,2),b7(j,3),b7(j,4),b7(j,5))))
  }
  A_b8 <- function(w_2) {
    t(matrix(c(b8(w_2,1), b8(w_2,2), b8(w_2,3), b8(w_2,4), b8(w_2,5))))
  }
  A_b9 <- function(w_2) {
    t(matrix(c(b9(w_2,1), b9(w_2,2), b9(w_2,3), b9(w_2,4))))
  }
  A_b10 <- function(j) {
    t(matrix(c(b10(j,1),b10(j,2),b10(j,3),b10(j,4))))
  }
  A_b11 <- function(j) {
    t(matrix(c(b11(j,1),b11(j,2),b11(j,3),b11(j,4),b11(j,5))))
  }
  A_b12 <- function(j) {
    t(matrix(c(b12(j,1),b12(j,2),b12(j,3),b12(j,4))))
  }
  
  #
  
}

#################################.
##### Modified exp function #####
#################################.

#' #' Modified exp function (see scratch for derivation)
#' #' 
#' #' @param x Numeric input
#' #' @return Numeric output
#' exp2 <- (function() {
#'   
#'   expit <- function(x) {1/(1+exp(-x))}
#'   logit <- function(x) { log(x/(1-x)) }
#'   e <- -0.1 # This value configurable but hard-coded
#'   ell <- logit(exp(e))
#'   x_0 <- e - (ell*exp(ell))/(exp(e)*(1+exp(ell))^2)
#'   k_0 <- exp(e-ell)*(1+exp(ell))^2
#'   exp2 <- function(x) {
#'     if (x<=e) {
#'       return(exp(x))
#'     } else {
#'       return(1/(1+exp(k_0*(x_0-x))))
#'     }
#'   }
#'   return(exp2)
#' })()



######################################################.
##### Simpson's sum approximation (did not work) #####
######################################################.

if (F) {
  
  # len_set <- length(d$X_i_set)
  # if (len_set>=7) {
  # 
  #   f2_first <- f2_fnc(d$X_i_set[[1]])
  #   f2_last <- f2_fnc(d$X_i_set[[len_set]])
  #   f2_mid_1 <- f2_fnc(d$X_i_set[[2]])
  #   f2_mid_2 <- f2_fnc(d$X_i_set[[round((len_set+1)/2)]])
  #   f2_mid_3 <- f2_fnc(d$X_i_set[[len_set-1]])
  # 
  #   simpson_sum <- (1/6)*(f2_mid_1+4*f2_mid_2+f2_mid_3)
  #   f2 <- sum(f2_first+f2_last+(len_set-2)*simpson_sum)
  # 
  # } else {
  
}



##########################.
##### Code profiling #####
##########################.

if (F) {
  
  cl <- parallel::makeCluster(cfg$sim_n_cores)
  parallel::clusterEvalQ(cl, sink(paste0("C:/Users/ak811/Desktop/Avi/Research/HIVMI/output", Sys.getpid(), ".txt"))) # !!!!!
  parallel::parLapply(cl, c(1:5), function(i) {
    print("Check 1")
    print(pryr::mem_used())
  })
  parallel::clusterExport(cl, c("f_x", "f_y", "icll"), envir=.GlobalEnv)
  parallel::parLapply(cl, c(1:5), function(i) {
    print("Check 2")
    print(pryr::mem_used())
  })
  parallel::clusterExport(cl, c("lik_fn2"), envir=.GlobalEnv)
  parallel::parLapply(cl, c(1:5), function(i) {
    print("Check 3")
    print(pryr::mem_used())
  })
  parallel::clusterExport(cl, c("inds", "batches"), envir=.GlobalEnv)
  parallel::parLapply(cl, c(1:5), function(i) {
    print("Check 4")
    print(pryr::mem_used())
  })
  print("Check 4.2")
  print("pryr::object_size(dat_objs)")
  print(pryr::object_size(dat_objs))
  print("pryr::object_size(get('dat_objs', envir=.GlobalEnv))")
  print(pryr::object_size(get("dat_objs", envir=.GlobalEnv)))
  parallel::clusterExport(cl, c("dat_objs"), envir=.GlobalEnv)
  parallel::parLapply(cl, c(1:5), function(i) {
    print("Check 5")
    print(pryr::mem_used())
    print("pryr::object_size(dat_objs)")
    print(pryr::object_size(dat_objs))
    # dat_objs[[i]]
  })
  
  parallel::clusterExport(cl, c("dat"), envir=.GlobalEnv)
  parallel::parLapply(cl, c(1:5), function(i) {
    pryr::object_size(dat)
    # pryr::mem_used()
  })
  
  par_init <- c(a_x=-6.535, g_x1=-0.6737, g_x2=3.636, g_x3=0.2734, g_x4=0.4366, g_x5=-8.512, t_x1=-1.578, t_x2=-0.2818, t_x3=0.3750, t_x4=-2.160, a_s=-3.369, g_s1=-0.5761, g_s2=0.8899, t_s1=0, t_s2=0, t_s3=0, t_s4=0, beta_x=1, beta_z=0.5072, a_y=-5.944, g_y1=0.3940, g_y2=1.871, g_y3=2.923, g_y4=6.809, g_y5=3.004, t_y=-0.6077)
  
  negloglik <- construct_negloglik(parallelize=T, cfg$model_version)
  system.time({ qqq1 <- negloglik(par_init) })
  print(qqq1)
  
  negloglik <- construct_negloglik(parallelize=F, cfg$model_version)
  system.time({ qqq2 <- negloglik(par_init) })
  print(qqq2)
  
}

#############################.
##### Old plotting code #####
#############################.

if (F) {
  
  #' Return plot of modeled mortality probabilities (HIV- vs. HIV+) with CIs
  #' @param sex Boolean integer; male=1, female=0
  #' @param age An integer; calendar year
  #' @param m An integer representing the model version number
  #' @param w_start An integer representing the window start calendar year
  #' @param y_max Maximum Y value for the plot
  #' @return ggplot2 object
  plot_mort2 <- function(sex, age, m, w_start, y_max) {
    
    if (w_start==2000) {
      grid <- seq(2000,2023,0.01) %>% (function(x) { x-(w_start-1) })
    } else if (w_start==2010) {
      grid <- seq(2010,2023,0.01) %>% (function(x) { x-(w_start-1) })
    }
    
    plot_data <- data.frame(
      x = rep(grid,2) + (w_start-1),
      Probability = 1000 * c(
        sapply(grid, function(j) { prob(type="mort (HIV-)", m=m, j=j, w_1=sex, w_2=age, which="est") }),
        sapply(grid, function(j) { prob(type="mort (HIV+)", m=m, j=j, w_1=sex, w_2=age, which="est") })
      ),
      ci_lo = 1000 * c(
        sapply(grid, function(j) { prob(type="mort (HIV-)", m=m, j=j, w_1=sex, w_2=age, which="ci_lo") }),
        sapply(grid, function(j) { prob(type="mort (HIV+)", m=m, j=j, w_1=sex, w_2=age, which="ci_lo") })
      ),
      ci_up = 1000 * c(
        sapply(grid, function(j) { prob(type="mort (HIV-)", m=m, j=j, w_1=sex, w_2=age, which="ci_up") }),
        sapply(grid, function(j) { prob(type="mort (HIV+)", m=m, j=j, w_1=sex, w_2=age, which="ci_up") })
      ),
      color = rep(c("HIV-","HIV+"), each=length(grid))
    )
    
    if (sex==0) {
      title <- paste0("Mortality rate among females aged ", age, "; model")
    } else if (sex==1) {
      title <- paste0("Mortality rate among males aged ", age, "; model")
    }
    title <- paste0(title, " v", m)
    
    plot <- ggplot(
      plot_data,
      aes(x=x, y=Probability, color=color)
    ) +
      geom_line() +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_up, fill=color),
        alpha = 0.2,
        linetype = "dotted"
      ) +
      coord_cartesian(ylim=c(0,y_max)) +
      labs(title=title, x="Year", color="HIV Status", fill="HIV Status",
           y="Deaths per 1,000 person-years") +
      theme(
        plot.background = element_rect(color="black"),
        legend.position = "bottom"
      ) +
      scale_color_manual(values=c("forestgreen", "#56B4E9")) +
      scale_fill_manual(values=c("forestgreen", "#56B4E9"))
    
    return(plot)
    
  }
  
  p11 <- plot_mort2(sex=0, age=20, m=m, w_start=w_start, y_max=y_max)
  p12 <- plot_mort2(sex=0, age=35, m=m, w_start=w_start, y_max=y_max)
  p13 <- plot_mort2(sex=0, age=50, m=m, w_start=w_start, y_max=y_max)
  p14 <- plot_mort2(sex=1, age=20, m=m, w_start=w_start, y_max=y_max)
  p15 <- plot_mort2(sex=1, age=35, m=m, w_start=w_start, y_max=y_max)
  p16 <- plot_mort2(sex=1, age=50, m=m, w_start=w_start, y_max=y_max)
  plot_05 <- ggpubr::ggarrange(
    p11, p12, p13, p14, p15, p16,
    ncol = 3,
    nrow = 2,
    legend = "bottom",
    common.legend = T
  )
  
}



##################################.
##### Old negloglik function #####
##################################.

if (F) {
  
  #' Negative log-likelihood across individuals and time
  #' @param dat A dataset returned by generate_dataset()
  #' @param par Vector of parameters governing the distribution.
  #' @return Numeric likelihood
  #' @notes This corresponds to the missing data structure
  negloglik_miss_old <- function(dat, par) {
    
    # Convert parameter vector to a named list
    p <- as.numeric(par)
    params <- list(a_x=p[1], g_x=c(p[2],p[3]), a_y=p[4], g_y=c(p[5],p[6]),
                   beta_x=p[7], beta_z=p[8])
    
    # Compute the negative likelihood across individuals
    n <- attr(dat, "n")
    -1 * sum(log(unlist(lapply(c(1:n), function(i) {
      
      dat_i <- filter(dat, id==i)
      J <- nrow(dat_i)
      
      # Calculate vectors for patient i
      w <- subset(dat_i, select=c(w_1,w_2))
      y <- dat_i$y
      z <- dat_i$z
      v <- dat_i$v
      u <- dat_i$u
      d <- dat_i$d
      xs <- dat_i$xs
      u_prev <- c(0, dat_i$u[c(1:(J-1))])
      
      # Calculate the set X_i to sum over
      X_i_set <- list()
      for (j in c(1:(J+1))) {
        x_ <- c(rep(0,J-j+1), rep(1,j-1))
        case_i_ <- case(x_,v)
        T_pm <- T_plusminus(case_i_, J, x_, v)
        if (case_i_!=9 &&
            all(xs==x_*d) &&
            all(d==g_delta(case_i_, J, T_pm$T_minus, T_pm$T_plus)) &&
            prod(u[1:J]==u_prev[1:J]+(1-u_prev[1:J])*v[1:J]*x_[1:J])
        ) {
          X_i_set <- c(X_i_set, list(x_))
        }
      }
      
      # Compute the likelihood for individual i
      w <- t(w)
      f2 <- sum(unlist(lapply(X_i_set, function(x) {
        x_prev <- c(0,x[1:(length(x)-1)])
        prod(unlist(lapply(c(1:J), function(j) {
          w_ij <- as.numeric(w[,j])
          f_x(x=x[j], x_prev=x_prev[j], w=w_ij, params=params) *
            f_y(y=y[j], x=x[j], w=w_ij, z=z[j], params=params)
        })))
      })))
      if (f2<=0) {
        f2 <- 1e-10
        # warning("Likelihood of zero")
      }
      
      return(f2)
      
    }))))
    
  }
  
}



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
