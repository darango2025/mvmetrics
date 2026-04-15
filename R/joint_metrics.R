# ============================================================
# mvmetrics — Joint (multivariate) metric functions
# ============================================================

#' Metric Family B — Multivariate Energy Score
#'
#' The Energy Score generalises the CRPS to R^d and is **strictly proper**.
#' It penalises both marginal miscalibration and errors in the joint
#' dependence structure.
#'
#' The binary component is embedded in R via the posterior predictive
#' probability (each draw is Bernoulli(p_hat)), giving a hybrid Euclidean
#' distance that remains computable and strictly proper under mild conditions.
#'
#' **Input scale:** `y_obs` and `y_draws` must be **normalised** (zero mean,
#' unit variance per variable) before calling this function. Use
#' [normalise_matrix()] with training-set parameters.
#'
#' @param y_obs   Numeric matrix \[N x K\]. Normalised observations.
#' @param y_draws Numeric array \[N x K x B\]. Normalised predictive draws.
#'   B = number of posterior predictive samples per observation.
#'
#' @return Scalar Energy Score (lower = better).
#'
#' @references
#' Gneiting, T. & Raftery, A.E. (2007). Strictly proper scoring rules,
#' prediction, and estimation. \emph{JASA}, 102(477), 359-378.
#'
#' @examples
#' set.seed(1)
#' N <- 50; K <- 3; B <- 20
#' Y      <- matrix(rnorm(N * K), N, K)
#' Ydraws <- array(rnorm(N * K * B), dim = c(N, K, B))
#' metric_B_energy_score(Y, Ydraws)
#'
#' @export
metric_B_energy_score <- function(y_obs, y_draws) {
  y_obs <- as.matrix(y_obs)
  if (length(dim(y_draws)) != 3)
    stop("y_draws must be a 3-D array [N x K x B].")
  N <- nrow(y_obs); B <- dim(y_draws)[3]

  es_vals <- vapply(seq_len(N), function(i) {
    yi   <- y_obs[i, ]
    Di   <- y_draws[i, , ]                        # K x B

    # E[||Y - y||]
    diffs1 <- Di - yi
    term1  <- mean(sqrt(colSums(diffs1^2)))

    # (1/2) E[||Y - Y'||]
    idx    <- utils::combn(B, 2)
    diffs2 <- Di[, idx[1, ], drop = FALSE] - Di[, idx[2, ], drop = FALSE]
    term2  <- mean(sqrt(colSums(diffs2^2)))

    term1 - 0.5 * term2
  }, numeric(1))

  mean(es_vals)
}

#' Metric Family C — Variogram Score
#'
#' The Variogram Score penalises errors in the cross-variable covariance
#' structure. It is **proper but not strictly proper**. Its practical advantage
#' is that it is more directly sensitive to pairwise covariance errors than the
#' Energy Score.
#'
#' **Input scale:** same normalisation requirement as [metric_B_energy_score()].
#'
#' @param y_obs   Numeric matrix \[N x K\]. Normalised observations.
#' @param y_draws Numeric array \[N x K x B\]. Normalised predictive draws.
#' @param p       Numeric. Order parameter (default = 0.5, recommended by
#'   Scheuerer & Hamill 2015).
#'
#' @return Scalar Variogram Score (lower = better).
#'
#' @references
#' Scheuerer, M. & Hamill, T.M. (2015). Variogram-based proper scoring rules
#' for probabilistic forecasts of multivariate quantities.
#' \emph{Monthly Weather Review}, 143(4), 1321-1334.
#'
#' @examples
#' set.seed(1)
#' N <- 30; K <- 3; B <- 15
#' Y      <- matrix(rnorm(N * K), N, K)
#' Ydraws <- array(rnorm(N * K * B), dim = c(N, K, B))
#' metric_C_variogram_score(Y, Ydraws)
#'
#' @export
metric_C_variogram_score <- function(y_obs, y_draws, p = 0.5) {
  y_obs <- as.matrix(y_obs)
  K     <- ncol(y_obs)
  N     <- nrow(y_obs)
  pairs <- utils::combn(K, 2)

  vs_vals <- vapply(seq_len(N), function(i) {
    yi  <- y_obs[i, ]
    Di  <- y_draws[i, , ]
    vs_i <- 0
    for (j in seq_len(ncol(pairs))) {
      k1 <- pairs[1, j]; k2 <- pairs[2, j]
      obs_diff  <- abs(yi[k1] - yi[k2])^p
      pred_diff <- mean(abs(Di[k1, ] - Di[k2, ])^p)
      vs_i <- vs_i + (obs_diff - pred_diff)^2
    }
    vs_i
  }, numeric(1))

  mean(vs_vals)
}

#' Metric Family D — Marginal-Dependence Decomposition (MDD)
#'
#' Separates predictive performance into a **marginal component** (S_marg,
#' from Metric A) and a **dependence component** (S_dep, from the Frobenius
#' distance between observed and model PIT correlation matrices). This is the
#' most scientifically informative metric: if a model beats the independent
#' baseline primarily through S_dep, that is direct evidence that its shared
#' structure captures genuine cross-variable coupling.
#'
#' **Input scale:** `y_obs`, `y_pred`, and `y_draws` in **original** (physical)
#' scale. The PIT transform is scale-invariant.
#'
#' @param y_obs   Numeric matrix \[N x K\]. Observations (original scale).
#' @param y_pred  Numeric matrix \[N x K\]. Point predictions (original scale).
#'   Binary column must contain probabilities.
#' @param y_draws Numeric array \[N x K x B\]. Predictive draws (original scale).
#' @param K_c     Integer. Number of continuous variables (default = K - 1).
#' @param alpha   Numeric in \[0, 1\]. Weight on the marginal component (default 0.5).
#' @param w_bin   Numeric. Weight for binary Log-Loss in S_marg (default 1.0).
#'
#' @return A list with components:
#'   \describe{
#'     \item{total}{Composite MDD score = alpha * S_marg + (1-alpha) * S_dep.}
#'     \item{S_marg}{Marginal score component (Metric A total).}
#'     \item{S_dep}{Dependence score (Frobenius distance of PIT correlations).}
#'     \item{C_obs}{K x K Spearman correlation matrix of observed PITs.}
#'     \item{C_model}{K x K Spearman correlation matrix of model PITs.}
#'   }
#'
#' @references
#' Ziel, F. & Berk, K. (2019). Multivariate forecasting evaluation: On
#' sensitive and strictly proper scoring rules. \emph{arXiv:1910.07325}.
#'
#' @examples
#' set.seed(1)
#' N <- 80; K <- 3; B <- 20
#' Y      <- cbind(matrix(rnorm(N * 2), N, 2), rbinom(N, 1, 0.4))
#' Yhat   <- cbind(Y[, 1:2] + rnorm(N * 2, 0, 0.3), runif(N, 0.3, 0.7))
#' Ydraws <- array(rnorm(N * K * B), dim = c(N, K, B))
#' Ydraws[, 3, ] <- rbinom(N * B, 1, 0.4)
#' metric_D_mdd(Y, Yhat, Ydraws, K_c = 2)
#'
#' @export
metric_D_mdd <- function(y_obs, y_pred, y_draws,
                          K_c = NULL, alpha = 0.5, w_bin = 1.0) {
  chk   <- .check_dims(y_obs, y_pred, y_draws)
  y_obs <- chk$y_obs; y_pred <- chk$y_pred
  K     <- ncol(y_obs); N <- nrow(y_obs); B <- dim(y_draws)[3]
  if (is.null(K_c)) K_c <- K - 1L

  # Marginal part (Metric A)
  S_marg <- metric_A(y_obs, y_pred, K_c = K_c, w_bin = w_bin)$total

  # PIT matrices
  PIT_obs   <- matrix(NA_real_, N, K)
  PIT_model <- matrix(NA_real_, N, K)

  for (k in seq_len(K_c)) {
    draws_k        <- y_draws[, k, ]               # N x B
    PIT_obs[, k]   <- rowMeans(draws_k <= y_obs[, k], na.rm = TRUE)
    rk_mat         <- apply(draws_k, 1, function(d)
                        rank(d, ties.method = "average"))
    diag_idx       <- ((seq_len(N) - 1L) %% B) + 1L
    PIT_model[, k] <- rk_mat[cbind(diag_idx, seq_len(N))] / (B + 1)
  }
  p_hat              <- pmin(pmax(y_pred[, K], 1e-7), 1 - 1e-7)
  PIT_obs[, K]       <- p_hat
  PIT_model[, K]     <- p_hat

  C_obs   <- .clean_cor(PIT_obs)
  C_model <- .clean_cor(PIT_model)
  S_dep   <- norm(C_obs - C_model, type = "F")^2

  list(
    total   = alpha * S_marg + (1 - alpha) * S_dep,
    S_marg  = S_marg,
    S_dep   = S_dep,
    C_obs   = C_obs,
    C_model = C_model
  )
}

#' Metric Family E — Enhanced CV-Weighted Marginal-Dependence (E-CVWMD)
#'
#' A novel joint scoring metric proposed by Arango Londoño (2026). Addresses
#' three limitations of prior metrics simultaneously:
#' \enumerate{
#'   \item Equal variable weighting regardless of predictive difficulty
#'         (fixed by CV-derived weights).
#'   \item Dependence score requiring full predictive CDF specification
#'         (replaced by a residual-correlation Frobenius distance).
#'   \item Fixed marginal/dependence trade-off parameter
#'         (replaced by a data-adaptive alpha* from a multivariate
#'         independence test).
#' }
#'
#' **Scale contract:** `y_obs`, `y_pred`, and `y_draws` must be in the
#' **original physical scale**. The `mu_train` and `sd_train` arguments
#' supply training-set parameters used for CV computation and residual
#' standardisation.
#'
#' @param y_obs    Numeric matrix \[N x K\]. Observations (original scale).
#' @param y_pred   Numeric matrix \[N x K\]. Predictions (original scale).
#' @param y_draws  Numeric array \[N x K x B\]. Predictive draws (original scale).
#' @param mu_train Numeric vector \[K_c\]. Training-set column means (continuous vars).
#' @param sd_train Numeric vector \[K_c\]. Training-set column SDs (continuous vars).
#' @param K_c      Integer. Number of continuous variables (default = K - 1).
#' @param alpha    Numeric or NULL. Fixed weight override. If NULL (default),
#'   an adaptive alpha* is computed via a distance-correlation test.
#' @param run_dcor Logical. Whether to run the dcov.test from the \pkg{energy}
#'   package (default TRUE). Set FALSE for speed in bootstrap loops.
#'
#' @return A list with components:
#'   \describe{
#'     \item{total}{E-CVWMD composite score (lower = better).}
#'     \item{S_E_marg}{CV-weighted marginal score.}
#'     \item{S_E_dep}{Residual-correlation dependence score.}
#'     \item{alpha_star}{Effective alpha used (adaptive or fixed).}
#'     \item{p_mv}{p-value from distance-correlation test (NA if not run).}
#'     \item{rmse_k}{Named RMSE vector for continuous variables.}
#'     \item{acc_bin}{Binary classification accuracy.}
#'     \item{w_cont}{CV-derived weights for continuous variables.}
#'     \item{w_bin}{Weight for the binary variable.}
#'     \item{R_obs}{K x K Spearman correlation matrix of standardised residuals.}
#'     \item{R_model}{K x K Spearman correlation matrix of model residuals.}
#'   }
#'
#' @references
#' Arango Londoño, D. (2026). Joint Modeling of Precipitation and Temperature
#' in Colombia. Doctoral dissertation, Universidad Nacional de Colombia.
#'
#' @examples
#' set.seed(1)
#' N <- 80; K <- 3; B <- 20
#' Y      <- cbind(matrix(rnorm(N * 2, mean = c(20, 30), sd = c(2, 4)),
#'                        N, 2, byrow = FALSE),
#'                 rbinom(N, 1, 0.4))
#' Yhat   <- cbind(Y[, 1:2] + rnorm(N * 2, 0, 0.5), runif(N, 0.3, 0.7))
#' Ydraws <- array(rnorm(N * K * B), dim = c(N, K, B))
#' Ydraws[, 3, ] <- rbinom(N * B, 1, 0.4)
#' metric_E_cvwmd(Y, Yhat, Ydraws,
#'                mu_train = c(20, 30),
#'                sd_train = c(2, 4),
#'                K_c = 2, run_dcor = FALSE)
#'
#' @export
metric_E_cvwmd <- function(y_obs, y_pred, y_draws,
                            mu_train, sd_train,
                            K_c = NULL, alpha = NULL,
                            run_dcor = TRUE) {
  chk   <- .check_dims(y_obs, y_pred, y_draws)
  y_obs <- chk$y_obs; y_pred <- chk$y_pred
  K     <- ncol(y_obs); N <- nrow(y_obs); B <- dim(y_draws)[3]
  if (is.null(K_c)) K_c <- K - 1L

  if (length(mu_train) != K_c || length(sd_train) != K_c)
    stop("mu_train and sd_train must have length K_c.")

  # ── 1. CV-derived weights ─────────────────────────────────────────────────
  cv_k    <- sd_train / pmax(abs(mu_train), 1e-6)
  w_k_raw <- cv_k / sum(cv_k)
  w_bin   <- 1 / (K_c + 1)
  w_cont  <- w_k_raw * (K_c / (K_c + 1))

  # ── 2. RMSE per continuous variable ──────────────────────────────────────
  rmse_k <- marginal_rmse(y_obs[, seq_len(K_c), drop = FALSE],
                           y_pred[, seq_len(K_c), drop = FALSE])

  # ── 3. Accuracy for binary variable ──────────────────────────────────────
  acc_bin     <- marginal_accuracy(y_obs[, K], y_pred[, K])
  penalty_bin <- 1 - acc_bin

  # ── 4. CV-weighted marginal score ─────────────────────────────────────────
  S_E_marg <- sum(w_cont * rmse_k) + w_bin * penalty_bin

  # ── 5. Residual-correlation dependence score ─────────────────────────────
  resid_obs   <- matrix(NA_real_, N, K)
  resid_model <- matrix(NA_real_, N, K)

  for (k in seq_len(K_c)) {
    resid_obs[, k]   <- (y_obs[, k]  - y_pred[, k])       / pmax(sd_train[k], 1e-6)
    resid_model[, k] <- (rowMeans(y_draws[, k, ]) - y_pred[, k]) / pmax(sd_train[k], 1e-6)
  }
  resid_obs[, K]   <- y_obs[, K]           - y_pred[, K]
  resid_model[, K] <- rowMeans(y_draws[, K, ]) - y_pred[, K]

  R_obs   <- .clean_cor(resid_obs)
  R_model <- .clean_cor(resid_model)
  S_E_dep <- norm(R_obs - R_model, type = "F")^2

  # ── 6. Adaptive alpha* ────────────────────────────────────────────────────
  p_mv <- NA_real_

  if (is.null(alpha)) {
    if (run_dcor && requireNamespace("energy", quietly = TRUE)) {
      n_dcor   <- min(500L, N)
      idx_dcor <- sample.int(N, n_dcor)
      sub      <- resid_obs[idx_dcor, seq_len(K_c), drop = FALSE]
      ok       <- apply(sub, 1, function(r) all(is.finite(r)))
      if (sum(ok) >= 20L) {
        dcor_res <- tryCatch(
          energy::dcov.test(sub[ok, ], R = 199L, seed = 1L),
          error = function(e) list(p.value = NA_real_)
        )
        p_mv <- dcor_res$p.value
      }
    }
    alpha_star <- if (!is.na(p_mv) && is.finite(p_mv) && p_mv < 0.05)
                    0.5 + 0.3 * (1 - p_mv) else 0.5
  } else {
    alpha_star <- alpha
  }

  # ── 7. Composite score ────────────────────────────────────────────────────
  S_E <- alpha_star * S_E_marg + (1 - alpha_star) * S_E_dep

  list(
    total      = S_E,
    S_E_marg   = S_E_marg,
    S_E_dep    = S_E_dep,
    alpha_star = alpha_star,
    p_mv       = p_mv,
    rmse_k     = rmse_k,
    acc_bin    = acc_bin,
    w_cont     = w_cont,
    w_bin      = w_bin,
    R_obs      = R_obs,
    R_model    = R_model
  )
}
