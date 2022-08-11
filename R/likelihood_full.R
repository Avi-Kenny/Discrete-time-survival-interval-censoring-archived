#' Negative log-likelihood across individuals and time
#' @param dat A dataset returned by generate_dataset()
#' @param par Vector of parameters governing the distribution.
#' @return Numeric likelihood
#' @notes This corresponds to the full data (unknown) structure and is for
#'        testing purposes only
negloglik_full <- function(dat, par) {
  
  # Convert parameter vector to a named list
  p <- par
  params <- list(a_x=p[1], g_x=c(p[2],p[3]), a_y=p[4], g_y=c(p[5],p[6]),
                 beta=p[7])
  
  # Calculate negative log likelihood
  -1 * sum(apply(
    X = dat,
    MARGIN = 1,
    FUN = function(r) {
      loglik_full_ij(
        x_j = r[["x"]],
        x_j1 = r[["x_prev"]],
        y_j = r[["y"]],
        w_j = c(r[["w_sex"]], r[["w_age"]]),
        params = params
      )
    }
  ))
  
}



#' Log-likelihood for a single individual
#' @param x_j X_j Serostatus at time j
#' @param x_j1 X_{j-1} serostatus at time j-1
#' @param y_j Y_j outcome indicator at time j
#' @param w_j W_j covariate vector at time j
#' @param params Named list of parameters
#' @return Numeric likelihood
#' @notes This corresponds to the full data (unknown) structure and is for
#'        testing purposes only
loglik_full_ij <- function(x_j, x_j1, y_j, w_j, params) {
  
  p <- params
  
  if (x_j1==1) {
    piece_x <- 1
  } else {
    exp_lin_x <- exp(p$a_x + sum(p$g_x*w_j))
    piece_x <- ifelse(x_j==1, exp_lin_x, 1-exp_lin_x)
  }
  
  exp_lin_y <- exp(p$a_y + sum(p$g_y*w_j) + p$b*x_j)
  piece_y <- ifelse(y_j==1, exp_lin_y, 1-exp_lin_y)
  
  return(log(max(piece_x,1e-8))+log(max(piece_y,1e-8)))
  
}
