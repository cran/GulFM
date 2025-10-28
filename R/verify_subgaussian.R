#' Verify sub-Gaussian preservation
#'
#' Draws n i.i.d. N(0, I_m) latent factors, applies g component-wise,
#' and checks whether E[exp(g(x))] remains below an empirical cut-off.
#' This is a quick proxy for finite sub-Gaussian norm.
#'
#' @param g_fun  vectorised map g: R -> R
#' @param m      latent dimension
#' @param n      Monte-Carlo sample size
#' @param cut    empirical threshold (default exp(2) & 7.389)
#' @return       logical TRUE if E[exp(g)] < cut on all coords
#' @examples
#'   tmp <- g_fun("strong_nonlinear")
#'   verify_subgaussian(tmp$g_fun, m = 5)
verify_subgaussian <- function(g_fun, m = 5, n = 1000, cut = exp(2)) {
  F1  <- MASS::mvrnorm(n = n, mu = rep(0, m), Sigma = diag(1, m))
  gF1 <- apply(F1, 2, g_fun)
  Eg  <- colMeans(exp(gF1))
  ok  <- all(Eg < cut)
  invisible(ok)
}
