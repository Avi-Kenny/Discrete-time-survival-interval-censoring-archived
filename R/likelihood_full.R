# Log-likelihood for a single individual
#' @param x_j X_j Serostatus at time j
#' @param x_j1 X_{j-1} serostatus at time j-1
#' @param y_j Y_j outcome indicator at time j
#' @param w_j W_j covariate vector at time j
#' @param prm Parameter vector: (alpha_j, g_x, beta_j, g_y, beta)
loglik_ij <- function(x_j, x_j1, y_j, w_j, prm) {
  
  if (x_j1==1) {
    piece_x <- 1
  } else {
    # !!!!! Add a "t" argument for linear baseline hazard
    exp_lin_x <- exp(prm[1] + prm[2]*w_j[1] + prm[3]*w_j[2])
    piece_x <- ifelse(x_j==1, exp_lin_x, 1-exp_lin_x)
  }
  
  exp_lin_y <- exp(prm[4] + prm[5]*w_j[1] + prm[6]*w_j[2] + prm[7]*x_j)
  piece_y <- ifelse(y_j==1, exp_lin_y, 1-exp_lin_y)
  
  return(log(max(piece_x,1e-8))+log(max(piece_y,1e-8)))
  
}


# Log-likelihood across individuals and time
negloglik <- function(prm) {
  
  -1 * sum(apply(
    X = dat,
    MARGIN = 1,
    FUN = function(r) {
      loglik_ij(
        x_j = r[["x"]],
        x_j1 = r[["x_prev"]],
        y_j = r[["y"]],
        w_j = c(r[["w_sex"]], r[["w_age"]]),
        prm = prm
      )
    }
  ))
  
}
