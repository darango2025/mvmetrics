# ============================================================
# mvmetrics — Internal utilities and normalisation helpers
# ============================================================

#' Normalise a matrix of observations or predictions
#'
#' Subtracts the column mean and divides by the column standard deviation,
#' placing all continuous variables on a comparable dimensionless scale before
#' computing joint metrics.
#'
#' @param mat   Numeric matrix \[N x K\]. Rows = observations, columns = variables.
#' @param mu    Numeric vector of length K. Column means (from training partition).
#' @param sigma Numeric vector of length K. Column standard deviations (training).
#'
#' @return Matrix of the same dimensions as `mat` with normalised values.
#'
#' @examples
#' Y <- matrix(rnorm(200), 100, 2)
#' Y_norm <- normalise_matrix(Y, mu = colMeans(Y), sigma = apply(Y, 2, sd))
#'
#' @export
normalise_matrix <- function(mat, mu, sigma) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)
  stopifnot(length(mu) == ncol(mat), length(sigma) == ncol(mat))
  sweep(sweep(mat, 2, mu, "-"), 2, pmax(sigma, 1e-10), "/")
}

# ── Internal helpers (not exported) ──────────────────────────────────────────

.clean_cor <- function(mat, method = "spearman") {
  ok <- apply(mat, 1, function(r) all(is.finite(r)))
  K  <- ncol(mat)
  if (sum(ok) < max(5, K + 1)) return(diag(K))
  C  <- stats::cor(mat[ok, ], method = method, use = "complete.obs")
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  C
}

.check_dims <- function(y_obs, y_pred, y_draws = NULL) {
  y_obs  <- as.matrix(y_obs)
  y_pred <- as.matrix(y_pred)
  if (!identical(dim(y_obs), dim(y_pred)))
    stop("y_obs and y_pred must have the same dimensions [N x K].")
  if (!is.null(y_draws)) {
    if (length(dim(y_draws)) != 3)
      stop("y_draws must be a 3-dimensional array [N x K x B].")
    if (dim(y_draws)[1] != nrow(y_obs) || dim(y_draws)[2] != ncol(y_obs))
      stop("y_draws dimensions [N x K x B] must match y_obs [N x K].")
  }
  list(y_obs = y_obs, y_pred = y_pred)
}
