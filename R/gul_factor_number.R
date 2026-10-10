#' Two-stage eigenvalue-ratio estimator for the number of factors
#'
#' This function performs a two-stage eigenvalue-ratio (ER) estimation
#' of the number of factors in a General Factor Model (GFM), tailored
#' to the GUL's two-layer structure.
#'
#' @param X Observation matrix with dimension \eqn{n \times p},
#'   where each row represents an observation.
#' @param k_max Upper bound for the first-stage search.
#'   Defaults to \eqn{\min(n-1, p-1)}.
#' @param g_fun Known smooth transformation applied element-wise to
#'   the first-layer factor scores. Defaults to \code{identity}.
#' @param p_norm Logical. If \code{TRUE}, normalize first-layer scores
#'   by \eqn{\sqrt{p}} before forming the second-layer covariance.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{m_hat_1}: First-stage ER estimate, an integer bounded
#'     above by \eqn{\min(n-1, p-1)}.
#'   \item \code{m_hat}: Second-stage refined ER estimate.
#'   \item \code{lambda_hat}: Top eigenvalues of the first-layer sample
#'     covariance \eqn{\hat\Sigma_X}.
#'   \item \code{mu_hat}: Eigenvalues of the second-layer sample
#'     covariance \eqn{\hat\Sigma^*_{g(F_1)}}.
#'   \item \code{ratio_1}: First-stage eigenvalue ratios
#'     \eqn{\hat\lambda_k / \hat\lambda_{k+1}}.
#'   \item \code{ratio_2}: Second-stage eigenvalue ratios
#'     \eqn{\hat\mu_k / \hat\mu_{k+1}}.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500; p <- 200; m_true <- 5
#' F  <- matrix(stats::rnorm(n * m_true), n, m_true)
#' gF <- sin(F)
#' A1 <- matrix(stats::rnorm(p * m_true, 0, 1/sqrt(p)), p, m_true)
#' V1 <- matrix(stats::rnorm(n * p, sd = 2), n, p)
#' X  <- gF %*% t(A1) + V1
#' res <- gul_factor_number(X, g_fun = sin)
#' res$m_hat
gul_factor_number <- function(X, k_max = NULL,
                              g_fun = identity,
                              p_norm = TRUE) {

  n <- nrow(X); p <- ncol(X)

  if (is.null(k_max)) k_max <- min(n - 1, p - 1)
  k_max <- min(k_max, n - 1, p - 1)

  ## Stage 1: ER on sample covariance of X
  svdX <- svd(X, nu = k_max, nv = k_max)
  lam  <- svdX$d^2 / (n - 1)
  ratio_1 <- lam[-length(lam)] / lam[-1]
  ratio_1[!is.finite(ratio_1)] <- 0
  m_hat_1 <- max(1L, which.max(ratio_1))

  ## Stage 2: if m_hat_1 <= 1, no second-stage ER possible
  if (m_hat_1 <= 1L) {
    return(list(m_hat_1 = m_hat_1, m_hat = m_hat_1,
                lambda_hat = lam, mu_hat = numeric(0),
                ratio_1 = ratio_1, ratio_2 = numeric(0)))
  }

  ## First-layer factor scores (revised normalization: Q1 only)
  Q1 <- svdX$v[, 1:m_hat_1, drop = FALSE]
  F1_hat <- X %*% Q1

  ## Apply known transformation g(.)
  gF1_hat <- g_fun(F1_hat)

  ## Normalize by sqrt(p) if requested
  if (p_norm) gF1_hat <- gF1_hat / sqrt(p)

  ## Second-layer sample covariance
  gF1_hat_c <- scale(gF1_hat, center = TRUE, scale = FALSE)
  Sigma2 <- crossprod(gF1_hat_c) / (n - 1)

  ## Eigen-decomposition
  eig2 <- eigen(Sigma2, symmetric = TRUE)
  mu   <- eig2$values

  ## Second-stage ER with convention mu_{m1+1} = 0
  ratio_2 <- mu[-length(mu)] / mu[-1]
  ratio_2[!is.finite(ratio_2)] <- 0
  m_hat <- max(1L, which.max(ratio_2))

  return(list(
    m_hat_1    = m_hat_1,
    m_hat      = m_hat,
    lambda_hat = lam,
    mu_hat     = mu,
    ratio_1    = ratio_1,
    ratio_2    = ratio_2
  ))
}
