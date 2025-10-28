#' Verify zero-mean preservation (Theorem 10 assumption 2a)
#'
#' Draws n i.i.d. N(0, I_m) latent factors, applies g component-wise,
#' and checks whether |E[g(x)]| < tol on every coordinate.
#'
#' @param g_fun  vectorised map g: R -> R
#' @param m      latent dimension
#' @param n      Monte-Carlo sample size
#' @param tol    numerical tolerance (default 1e-3)
#' @return       logical TRUE if |mean| < tol on all coords
#' @examples
#'   tmp <- g_fun("weak_nonlinear")
#'   verify_mean(tmp$g_fun, m = 5)
verify_mean <- function(g_fun, m = 5, n = 10000, tol = 1e-3) {
  F1  <- MASS::mvrnorm(n = n, mu = rep(0, m), Sigma = diag(1, m))
  gF1 <- apply(F1, 2, g_fun)
  mg  <- colMeans(gF1)
  invisible(all(abs(mg) < tol))
}
