#' Simulation wrapper for Theorems 9 & 10
#'
#' One Monte-Carlo replicate; returns empirical error, exceedance indicator,
#' theoretical bounds, and assumption-check flags.
#'
#' @param n       sample size
#' @param p       number of observed variables
#' @param m       number of latent factors
#' @param g_type  character: "linear", "weak_nonlinear", "strong_nonlinear"
#' @param epsilon error threshold
#' @param zero_tol zero-mean tolerance (default 0.02)
#' @return        one-row data-frame
#' @examples
#'   df <- g_theorem(500, 200, 5, "linear", 0.6)
g_theorem <- function(n, p, m, g_type, epsilon, zero_tol = 0.02) {
  ## 1. pick smooth link and its Lipschitz constant
  g_info <- g_fun(g_type)
  g_fun  <- g_info$g_fun
  Lg     <- g_info$L_g

  ## 2. verify Theorem-10 assumptions (quiet)
  zero_ok <- verify_mean(g_fun, m = m, tol = zero_tol)
  subG_ok <- verify_subgaussian(g_fun, m = m)

  ## 3. generate GFM data
  dat      <- generate_gfm_data(n = n, p = p, m = m, g_fun = g_fun)
  Ag_true  <- dat$Ag

  ## 4. GUL estimate
  Ag_hat   <- estimate_gul_loadings(X = dat$X, m = m)$hat_Ag

  ## 5. error and empirical probability (1 replicate)
  error_F   <- norm(Ag_hat - Ag_true, "F")
  prob_exceed <- as.integer(error_F > epsilon)   # 0/1

  ## 6. theoretical bounds (empirical C = 15)
  C <- 15
  bound9  <- C * Lg^2 * m * (m + log(p)) / (n * epsilon^2)
  bound10 <- C * Lg^2 * m * log(p)           / (n * epsilon^2)

  ## 7. return everything
  data.frame(n = n, p = p, m = m, g_type = g_type, L_g = Lg,
             epsilon = epsilon, error_F = error_F,
             prob_exceed = prob_exceed,
             bound_thm9  = bound9,
             bound_thm10 = bound10,
             zero_ok     = zero_ok,
             subG_ok     = subG_ok)
}
