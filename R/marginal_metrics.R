# ============================================================
# mvmetrics — Marginal (univariate) metric functions
# ============================================================

#' Root Mean Squared Error per variable
#'
#' Computes RMSE for each column of a prediction matrix independently.
#'
#' @param y_obs  Numeric matrix \[N x K\]. Observed values.
#' @param y_pred Numeric matrix \[N x K\]. Predicted values.
#'
#' @return Named numeric vector of length K with RMSE per variable.
#'
#' @examples
#' Y <- matrix(rnorm(300), 100, 3)
#' Yhat <- Y + matrix(rnorm(300, 0, 0.5), 100, 3)
#' marginal_rmse(Y, Yhat)
#'
#' @export
marginal_rmse <- function(y_obs, y_pred) {
  chk <- .check_dims(y_obs, y_pred)
  y_obs <- chk$y_obs; y_pred <- chk$y_pred
  out <- sapply(seq_len(ncol(y_obs)), function(k)
    sqrt(mean((y_obs[, k] - y_pred[, k])^2, na.rm = TRUE)))
  names(out) <- colnames(y_obs)
  out
}

#' Mean Absolute Error per variable
#'
#' @param y_obs  Numeric matrix \[N x K\]. Observed values.
#' @param y_pred Numeric matrix \[N x K\]. Predicted values.
#'
#' @return Named numeric vector of length K with MAE per variable.
#'
#' @examples
#' Y <- matrix(rnorm(300), 100, 3)
#' Yhat <- Y + matrix(rnorm(300, 0, 0.5), 100, 3)
#' marginal_mae(Y, Yhat)
#'
#' @export
marginal_mae <- function(y_obs, y_pred) {
  chk <- .check_dims(y_obs, y_pred)
  y_obs <- chk$y_obs; y_pred <- chk$y_pred
  out <- sapply(seq_len(ncol(y_obs)), function(k)
    mean(abs(y_obs[, k] - y_pred[, k]), na.rm = TRUE))
  names(out) <- colnames(y_obs)
  out
}

#' Log-Loss for binary variables
#'
#' Computes the binary cross-entropy (Log-Loss) for predicted probabilities.
#' Predicted values are clipped to \[1e-7, 1 - 1e-7\] to avoid log(0).
#'
#' @param y_obs  Numeric vector or matrix column of 0/1 observations.
#' @param p_pred Numeric vector or matrix column of predicted probabilities in (0,1).
#'
#' @return Scalar Log-Loss value (lower = better).
#'
#' @examples
#' y    <- rbinom(100, 1, 0.4)
#' phat <- runif(100, 0.2, 0.8)
#' marginal_logloss(y, phat)
#'
#' @export
marginal_logloss <- function(y_obs, p_pred) {
  y_obs  <- as.vector(y_obs)
  p_pred <- as.vector(p_pred)
  if (length(y_obs) != length(p_pred))
    stop("y_obs and p_pred must have the same length.")
  p_pred <- pmin(pmax(p_pred, 1e-7), 1 - 1e-7)
  -mean(y_obs * log(p_pred) + (1 - y_obs) * log(1 - p_pred), na.rm = TRUE)
}

#' Classification Accuracy for binary variables
#'
#' Threshold at 0.5: predicted probability >= 0.5 is classified as 1.
#'
#' @param y_obs  Numeric vector of 0/1 observations.
#' @param p_pred Numeric vector of predicted probabilities.
#'
#' @return Scalar accuracy in \[0, 1\] (higher = better).
#'
#' @examples
#' y    <- rbinom(100, 1, 0.4)
#' phat <- runif(100, 0.2, 0.8)
#' marginal_accuracy(y, phat)
#'
#' @export
marginal_accuracy <- function(y_obs, p_pred) {
  y_obs  <- as.vector(y_obs)
  p_pred <- as.vector(p_pred)
  mean(as.integer(p_pred >= 0.5) == y_obs, na.rm = TRUE)
}

#' Metric Family A — Weighted Sum of Marginal Scores
#'
#' Aggregates RMSE for continuous variables and Log-Loss for the binary variable
#' using uniform weights (1/K_c per continuous variable, w_bin for binary).
#' This metric is **blind to cross-variable dependence** by construction and
#' serves as the marginal baseline.
#'
#' @param y_obs  Numeric matrix \[N x K\]. Observations (continuous + binary).
#'   Continuous variables in columns 1:K_c; binary variable in column K.
#' @param y_pred Numeric matrix \[N x K\]. Point predictions. Binary column
#'   must contain probabilities in (0, 1).
#' @param K_c    Integer. Number of continuous variables (default = ncol(y_obs) - 1).
#' @param w_bin  Numeric. Weight for the binary Log-Loss term (default = 1.0).
#'
#' @return A list with components:
#'   \describe{
#'     \item{total}{Weighted aggregate score (lower = better).}
#'     \item{rmse_per_var}{Named vector of RMSE per continuous variable.}
#'     \item{rmse_mean}{Unweighted mean RMSE across continuous variables.}
#'     \item{logloss}{Log-Loss for the binary variable.}
#'   }
#'
#' @references
#' Koochali, A. et al. (2022). Random noise vs state-of-the-art probabilistic
#' forecasting methods. \emph{Applied Sciences}, 12(10), 5104.
#'
#' @examples
#' set.seed(1)
#' Y    <- cbind(matrix(rnorm(400), 100, 4), rbinom(100, 1, 0.4))
#' Yhat <- cbind(Y[, 1:4] + rnorm(400, 0, 0.3),
#'               runif(100, 0.3, 0.7))
#' metric_A(Y, Yhat, K_c = 4)
#'
#' @export
metric_A <- function(y_obs, y_pred, K_c = NULL, w_bin = 1.0) {
  chk  <- .check_dims(y_obs, y_pred)
  y_obs <- chk$y_obs; y_pred <- chk$y_pred
  K    <- ncol(y_obs)
  if (is.null(K_c)) K_c <- K - 1L
  if (K_c >= K) stop("K_c must be < ncol(y_obs); last column assumed binary.")

  rmse_k <- marginal_rmse(y_obs[, seq_len(K_c), drop = FALSE],
                           y_pred[, seq_len(K_c), drop = FALSE])
  ll     <- marginal_logloss(y_obs[, K], y_pred[, K])
  score  <- mean(rmse_k) + w_bin * ll

  list(total = score, rmse_per_var = rmse_k,
       rmse_mean = mean(rmse_k), logloss = ll)
}
