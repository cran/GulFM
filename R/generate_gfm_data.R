#' Generate general factor model with smooth latent transformation
#'
#' @param n      Integer: sample size.
#' @param p      Integer: number of observed variables.
#' @param m      Integer: number of latent factors (both layers).
#' @param seed   1.
#' @param g_fun  Function: smooth, element-wise transformation applied to latent
#'               factors.  Must be vectorised, e.g. `sin`, `tanh`, `scale`.
#' @param sigma_V Numeric: standard deviation of the idiosyncratic noise
#'               (default 0.1 => Var = 0.01).
#' @return List with components
#'   X  : n * p matrix of standardised observations.
#'   A1 : p * m first-layer loading matrix.
#'   A2 : m * m second-layer loading matrix.
#'   Ag : p * m overall loading matrix (Ag = A1 %*% t(A2)).
#'   F1 : n * m latent factors (before transformation).
#'   gF1: n * m latent factors (after transformation).
#'   V1 : n * p noise matrix (for diagnostics).
#'
#' @examples
#' dat <- generate_gfm_data(200, 50, 5, g_fun = tanh)
generate_gfm_data <- function(n, p, m, g_fun,seed=1,
                              sigma_V = 0.1) {

  ## 1.  Loadings -----------------------------------------------------------
  ## Draw rows i.i.d.  N_m(0, (1/m)I) so that A1A1^T  I_p when p is large.
  A1 <- MASS::mvrnorm(n = p, mu = rep(0, m), Sigma = diag(1/m, m))
  A2 <- MASS::mvrnorm(n = m, mu = rep(0, m), Sigma = diag(1/m, m))
  Ag <- A1 %*% t(A2)
  Ag <- Ag / norm(Ag, "2")          # scale so that spectral norm = 1

  ## 2.  Latent factors ------------------------------------------------------
  F1 <- MASS::mvrnorm(n = n, mu = rep(0, m), Sigma = diag(m))
  gF1 <- matrix(g_fun(F1), nrow = n, ncol = m)  # vectorised  fast

  ## 3.  Idiosyncratic noise -------------------------------------------------
  V1 <- MASS::mvrnorm(n = n, mu = rep(0, p),
                      Sigma = diag(sigma_V^2, p))

  ## 4.  Observed data -------------------------------------------------------
  X <- gF1 %*% t(A1) + V1

  ## 5.  Standardisation (eq. 2.1) ------------------------------------------
  X <- scale(X, center = TRUE, scale = TRUE)       # centre columns
  X <- X %*% diag(1 / sqrt(matrixStats::colVars(X))) # unit column variance

  return(list(X = X, A1 = A1, A2 = A2, Ag = Ag,
       F1 = F1, gF1 = gF1, V1 = V1))
}

