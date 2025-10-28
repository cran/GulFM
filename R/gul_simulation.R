#' Single-replication GUL simulation
#'
#' Generates one synthetic data set, estimates loadings with the GUL, and evaluates estimation accuracy.
#'
#' @param n      Integer: sample size.
#' @param p      Integer: number of observed variables.
#' @param m      Integer: number of latent factors (both layers).
#' @param g_fun  Function: element-wise, smooth transformation applied to the
#'               latent factors (e.g. `tanh`, `sin`).
#'
#' @return Named numeric vector with components
#'   error_F      : Frobenius norm ||hat(Ag) - Ag||_F
#' @examples
#' gul_simulation(200, 50, 5, g_fun = tanh)
gul_simulation <- function(n, p, m, g_fun) {
  ## ---------- 1.  generate synthetic data ---------------------------------
  data      <- generate_gfm_data(n = n, p = p, m = m, g_fun = g_fun)
  X         <- data$X          # n x p centred, scaled data matrix
  Ag_true   <- data$Ag         # p x m true overall loading matrix

  ## ---------- 2.  estimate loadings via GUL --------------------------------
  gul_res   <- estimate_gul_loadings(X = X, m = m)
  Ag_hat    <- gul_res$hat_Ag  # p x m estimated overall loading matrix

  ## ---------- 3.  compute Frobenius error ----------------------------------
  error_F <- norm(Ag_hat - Ag_true, type = "F")  # base-R function

  ## ---------- 5.  return results -------------------------------------------
   return(list(error_F = error_F))
}
