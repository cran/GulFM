#' Multi-metric evaluation of factor loading matrix estimation error
#' @param A_true True loading matrix (p x m)
#' @param A_hat  Estimated loading matrix (p x m)
#' @return data.frame with MSE, RMSE, MAE, MaxDev, and Cosine similarity
#' @examples
#' ## simulated example
#' p <- 100; m <- 5
#' Ag_true  <- matrix(rnorm(p*m), p, m)
#' Ag_hat   <- Ag_true + matrix(rnorm(p*m, 0, 0.1), p, m)
#' metrics  <- loading_metrics(Ag_true, Ag_hat)
#' print(metrics)

loading_metrics <- function(A_true, A_hat) {
  if (!all(dim(A_true) == dim(A_hat)))
    stop("Dimensions of A_true and A_hat must match.")
  
  D <- A_true - A_hat
  p <- nrow(A_true)
  m <- ncol(A_true)
  
  MSE    <- mean(D^2)
  RMSE   <- sqrt(MSE)
  MAE    <- mean(abs(D))
  MaxDev <- max(abs(D))
  
  # Average column-wise cosine similarity
  cos_col <- function(u, v) sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
  Cosine  <- mean(mapply(cos_col,
                         asplit(A_true, 2),   # split into column vectors
                         asplit(A_hat,  2)))
  
  data.frame(MSE = MSE,
             RMSE = RMSE,
             MAE = MAE,
             MaxDev = MaxDev,
             Cosine = Cosine)
}