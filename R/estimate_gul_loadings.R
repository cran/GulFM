#'   General unilateral load  Estimator
#'
#' @param X   n *p data matrix (already centred and scaled if desired).
#' @param m   number of latent factors (both layers).
#'
#' @return A list with
#'   hat_A1  : p * m  1st-layer loadings
#'   hat_A2  : m * m  2nd-layer loadings
#'   hat_Ag  : p * m  overall loadings
#'   Sigma1  : p * p  sample cov(X)  (for diagnostics)
#'   Sigma2  : m * m  sample cov(hat_gF1)
#'   hat_gF1 : n * m  estimated transformed latent factors
#'   eig1    : eigen-values of Sigma1
#'   eig2    : eigen-values of Sigma2
#'
#' @details
#' Step 1: PCA on X  to get  hat_A1
#' Step 2: Regress X on hat_A1  to get  hat_gF1
#' Step 3: PCA on hat_gF1  to get  hat_A2
#' Step 4: hat_Ag = hat_A1 %*% t(hat_A2)
#'
#' @examples
#' dat  <- generate_gfm_data(500, 50, 5, tanh, seed = 1)
#' est  <- estimate_gul_loadings(dat$X, m = 5)
#' err  <- sqrt(mean((est$hat_Ag - dat$Ag)^2))  # overall RMSE
estimate_gul_loadings <- function(X, m) {

  n <- nrow(X)
  p <- ncol(X)

  ## ---------- Step 1: PCA on X to get hat_A1 ------------------------------
  Sigma1 <- cov(X)
  eig1 <- eigen(Sigma1, symmetric = TRUE)
  idx1 <- order(eig1$values, decreasing = TRUE)[1:m]
  Lambda1 <- diag(eig1$values[idx1])
  Q1 <- eig1$vectors[, idx1]
  hat_A1 <- Q1 %*% sqrt(Lambda1)

  ## ---------- Step 2: estimate transformed latent factors -----------------
  # Least-squares:  X %*% hat_gF1 %*% t(hat_A1)
  # hat_gF1 = X %*% hat_A1 %*% solve(t(hat_A1)%*%hat_A1)
  # Use chol2inv for numerical stability
  hat_gF1 <- X %*% hat_A1 %*% chol2inv(chol(crossprod(hat_A1)))

  ## ---------- Step 3: PCA on hat_gF1 to get hat_A2 ------------------------
  Sigma2 <- cov(hat_gF1)
  eig2 <- eigen(Sigma2, symmetric = TRUE)
  idx2 <- order(eig2$values, decreasing = TRUE)[1:m]
  Lambda2 <- diag(eig2$values[idx2])
  Q2 <- eig2$vectors[, idx2]
  hat_A2 <- Q2 %*% sqrt(Lambda2)

  ## ---------- Step 4: overall loading -------------------------------------
  hat_Ag <- hat_A1 %*% t(hat_A2)

  return(list(hat_A1 = hat_A1,
                 hat_A2 = hat_A2,
                 hat_Ag = hat_Ag,
                 Sigma1 = Sigma1,
                 Sigma2 = Sigma2,
                 hat_gF1 = hat_gF1,
                 eig1 = eig1$values,
                 eig2 = eig2$values))
}
