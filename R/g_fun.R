#' Smooth link functions compliant with Theorems 9&10
#'
#' Returns a vectorised map \eqn{g(\cdot)} and its exact Lipschitz constant
#' \eqn{L_g} for three increasingly nonlinear choices.
#'
#' @param type Character string selecting the map:
#'   \code{"linear"}, \code{"weak_nonlinear"}, or \code{"strong_nonlinear"}.
#'
#' @return Named list with components
#'   \item{g_fun}{vectorised function \eqn{g(\cdot)}}
#'   \item{L_g}{scalar Lipschitz constant of \eqn{g}}
#'
#' @examples
#' ## pick a link with L_g = 1
#' tmp  <- g_fun("linear")
#' dat  <- generate_gfm_data(n = 500, p = 200, m = 5, g_fun = tmp$g_fun)
#' est  <- estimate_gul_loadings(dat$X, m = 5)
#' err  <- norm(est$hat_Ag - dat$Ag, "F")
#' sprintf("F-error (L_g = %d) = %.3f", tmp$L_g, err)
g_fun <- function(type = c("linear",
                           "weak_nonlinear",
                           "strong_nonlinear")) {
  type <- match.arg(type)

  switch(type,
         "linear" = {
           g_fun <- function(x) x
           L_g   <- 1
         },
         "weak_nonlinear" = {
           g_fun <- function(x) 0.5 * tanh(2 * x)   # L_g=1
           L_g   <- 1
         },
         "strong_nonlinear" = {
           g_fun <- function(x) 2 * tanh(x)         # L_g=2
           L_g   <- 2
         })

  g_fun <- Vectorize(g_fun)
  list(g_fun = g_fun, L_g = L_g)
}
