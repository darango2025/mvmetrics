# ============================================================
# mvmetrics — Main evaluation wrapper
# ============================================================

#' Evaluate multiple models with all joint metric families
#'
#' The primary user-facing function of \pkg{mvmetrics}. Accepts a list of
#' observed/predicted matrix pairs (one entry per model), computes all five
#' metric families (A–E), produces a ranked summary table, and recommends the
#' metric family best suited to the data structure detected.
#'
#' @param obs      Numeric matrix \[N x K\]. Observed values in **original**
#'   physical scale. The last column is treated as binary (0/1) unless
#'   `K_c` is specified. Columns should be named for readable output.
#' @param pred_list Named list of predicted matrices. Each element is a numeric
#'   matrix \[N x K\] with the same column structure as `obs`. The list names
#'   become model labels (e.g. `list(M1 = ..., M2 = ...)`).
#' @param draws_list Named list of predictive draw arrays. Each element is a
#'   numeric array \[N x K x B\] in the same scale as `obs`. Must have the
#'   same names as `pred_list`. Required for Metrics B, C, D, E.
#'   Pass NULL to skip probabilistic metrics (only Metric A will be computed).
#' @param mu_train Numeric vector \[K_c\]. Training-set column means for
#'   continuous variables. Required for Metric E. If NULL, Metric E is skipped.
#' @param sd_train Numeric vector \[K_c\]. Training-set column SDs.
#'   Required for Metric E.
#' @param K_c     Integer. Number of continuous variables (default = K - 1).
#' @param alpha_D Numeric in \[0, 1\]. Fixed alpha for Metric D (default 0.5).
#' @param alpha_E Numeric or NULL. Fixed alpha for Metric E. NULL = adaptive.
#' @param run_dcor Logical. Run distance-correlation test for Metric E alpha*?
#'   (default FALSE for speed; set TRUE for final analyses).
#' @param bootstrap Logical. Compute bootstrap 95% CIs? (default FALSE).
#' @param B_boot Integer. Number of bootstrap resamples (default 200).
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return An object of class `"mvmetrics_report"` — a list with:
#'   \describe{
#'     \item{scores}{data.frame with one row per model and columns for each
#'       metric (A, B, C, D, D_marg, D_dep, E, E_marg, E_dep, E_alpha_star).}
#'     \item{rankings}{data.frame with rank per model per metric (1 = best).}
#'     \item{recommendation}{Character. Recommended metric family with rationale.}
#'     \item{bootstrap_ci}{data.frame of bootstrap CIs (NULL if bootstrap=FALSE).}
#'     \item{marginal_detail}{List of per-variable marginal metrics.}
#'     \item{call}{The matched call.}
#'   }
#'   Use [print.mvmetrics_report()], [summary.mvmetrics_report()], and
#'   [plot_metrics()] for formatted output.
#'
#' @examples
#' set.seed(42)
#' N <- 200; K <- 5; B <- 30
#'
#' # Observed: 4 continuous + 1 binary
#' mu  <- c(18.9, 30.6, 80.7, 3.8, NA)
#' sds <- c(2.5,  4.0,  8.0,  0.8, NA)
#' obs <- cbind(
#'   matrix(rnorm(N * 4, rep(mu[1:4], each = N), rep(sds[1:4], each = N)), N, 4),
#'   rbinom(N, 1, 0.35)
#' )
#' colnames(obs) <- c("Tmin","Tmax","HR","Rad","Pbin")
#'
#' # M1: independent predictions
#' m1_pred   <- cbind(matrix(rnorm(N * 4, rep(mu[1:4], each = N),
#'                                 rep(sds[1:4] * 1.1, each = N)), N, 4),
#'                    runif(N, 0.3, 0.5))
#' m1_draws  <- array(rnorm(N * K * B), dim = c(N, K, B))
#'
#' # M3: joint predictions (slightly better)
#' m3_pred   <- cbind(matrix(rnorm(N * 4, rep(mu[1:4], each = N),
#'                                 rep(sds[1:4] * 0.9, each = N)), N, 4),
#'                    runif(N, 0.35, 0.55))
#' m3_draws  <- array(rnorm(N * K * B), dim = c(N, K, B))
#' colnames(m1_pred) <- colnames(m3_pred) <- colnames(obs)
#'
#' report <- evaluate_models(
#'   obs        = obs,
#'   pred_list  = list(M1 = m1_pred, M3 = m3_pred),
#'   draws_list = list(M1 = m1_draws, M3 = m3_draws),
#'   mu_train   = mu[1:4],
#'   sd_train   = sds[1:4],
#'   K_c        = 4
#' )
#' print(report)
#'
#' @export
evaluate_models <- function(obs,
                             pred_list,
                             draws_list = NULL,
                             mu_train   = NULL,
                             sd_train   = NULL,
                             K_c        = NULL,
                             alpha_D    = 0.5,
                             alpha_E    = NULL,
                             run_dcor   = FALSE,
                             bootstrap  = FALSE,
                             B_boot     = 200L,
                             verbose    = TRUE) {

  mc <- match.call()

  # ── Input validation ───────────────────────────────────────────────────────
  obs <- as.matrix(obs)
  K   <- ncol(obs)
  N   <- nrow(obs)
  if (is.null(K_c)) K_c <- K - 1L
  if (!is.list(pred_list) || length(pred_list) == 0)
    stop("pred_list must be a non-empty named list.")
  if (is.null(names(pred_list)))
    names(pred_list) <- paste0("M", seq_along(pred_list))

  n_models <- length(pred_list)
  model_names <- names(pred_list)

  has_draws <- !is.null(draws_list)
  has_E     <- has_draws && !is.null(mu_train) && !is.null(sd_train)

  if (has_draws && !identical(names(draws_list), model_names))
    stop("draws_list must have the same names as pred_list.")

  # Training parameters for normalisation (needed for B and C)
  mu_norm  <- if (!is.null(mu_train)) mu_train else
                colMeans(obs[, seq_len(K_c), drop = FALSE], na.rm = TRUE)
  sd_norm  <- if (!is.null(sd_train)) sd_train else
                apply(obs[, seq_len(K_c), drop = FALSE], 2, stats::sd, na.rm = TRUE)

  # Normalised observations
  obs_norm        <- obs
  obs_norm[, seq_len(K_c)] <- normalise_matrix(
    obs[, seq_len(K_c), drop = FALSE], mu_norm, sd_norm)

  # ── Compute metrics for each model ────────────────────────────────────────
  scores_list  <- vector("list", n_models)
  marginal_det <- vector("list", n_models)
  names(scores_list) <- names(marginal_det) <- model_names

  for (m in model_names) {
    if (verbose) message("  Computing metrics for model: ", m)

    y_pred <- as.matrix(pred_list[[m]])
    y_pred_norm        <- y_pred
    y_pred_norm[, seq_len(K_c)] <- normalise_matrix(
      y_pred[, seq_len(K_c), drop = FALSE], mu_norm, sd_norm)

    sc <- list()

    # Metric A
    res_A <- metric_A(obs, y_pred, K_c = K_c)
    sc$A  <- res_A$total
    marginal_det[[m]] <- list(
      rmse_per_var = res_A$rmse_per_var,
      logloss      = res_A$logloss,
      accuracy     = marginal_accuracy(obs[, K], y_pred[, K])
    )

    if (has_draws) {
      y_drw      <- draws_list[[m]]
      y_drw_norm <- y_drw
      for (k in seq_len(K_c))
        y_drw_norm[, k, ] <- (y_drw[, k, ] - mu_norm[k]) / pmax(sd_norm[k], 1e-10)

      sc$B <- metric_B_energy_score(obs_norm, y_drw_norm)
      sc$C <- metric_C_variogram_score(obs_norm, y_drw_norm)

      res_D    <- metric_D_mdd(obs, y_pred, y_drw, K_c = K_c, alpha = alpha_D)
      sc$D     <- res_D$total
      sc$D_marg <- res_D$S_marg
      sc$D_dep  <- res_D$S_dep
    }

    if (has_E) {
      res_E         <- metric_E_cvwmd(obs, y_pred, y_drw,
                                       mu_train = mu_train, sd_train = sd_train,
                                       K_c = K_c, alpha = alpha_E,
                                       run_dcor = run_dcor)
      sc$E           <- res_E$total
      sc$E_marg      <- res_E$S_E_marg
      sc$E_dep       <- res_E$S_E_dep
      sc$E_alpha_star <- res_E$alpha_star
    }

    scores_list[[m]] <- sc
  }

  # ── Assemble scores data.frame ─────────────────────────────────────────────
  all_cols <- c("A","B","C","D","D_marg","D_dep","E","E_marg","E_dep","E_alpha_star")
  scores_df <- do.call(rbind, lapply(model_names, function(m) {
    sc  <- scores_list[[m]]
    row <- setNames(rep(NA_real_, length(all_cols)), all_cols)
    for (col in names(sc)) if (col %in% all_cols) row[col] <- sc[[col]]
    as.data.frame(t(row))
  }))
  rownames(scores_df) <- model_names

  # Remove all-NA columns
  keep <- apply(scores_df, 2, function(x) any(!is.na(x)))
  scores_df <- scores_df[, keep, drop = FALSE]

  # ── Rankings (rank by column, lower score = rank 1 = better) ──────────────
  rank_cols <- setdiff(names(scores_df),
                       c("D_marg","D_dep","E_marg","E_dep","E_alpha_star"))
  rank_df <- scores_df[, rank_cols, drop = FALSE]
  for (col in names(rank_df))
    rank_df[[col]] <- rank(rank_df[[col]], ties.method = "min")

  # ── Metric recommendation ──────────────────────────────────────────────────
  rec <- recommend_metric(scores_df, n_models = n_models,
                           has_draws = has_draws, has_E = has_E)

  # ── Bootstrap CIs ─────────────────────────────────────────────────────────
  boot_ci <- NULL
  if (bootstrap && has_draws) {
    if (verbose) message("  Running bootstrap (B = ", B_boot, ")...")
    boot_ci <- bootstrap_metrics(
      obs = obs, pred_list = pred_list, draws_list = draws_list,
      mu_train = mu_train, sd_train = sd_train,
      K_c = K_c, alpha_D = alpha_D, B_boot = B_boot, verbose = verbose
    )
  }

  structure(
    list(
      scores         = scores_df,
      rankings       = rank_df,
      recommendation = rec,
      bootstrap_ci   = boot_ci,
      marginal_detail = marginal_det,
      call           = mc
    ),
    class = "mvmetrics_report"
  )
}

# ── Print method ──────────────────────────────────────────────────────────────

#' Print an mvmetrics report
#' @param x   Object of class `mvmetrics_report`.
#' @param ... Ignored.
#' @export
print.mvmetrics_report <- function(x, ...) {
  cat("\n=== mvmetrics: Joint Performance Metric Report ===\n\n")

  cat("--- Point Estimates (lower = better) ---\n")
  print(round(x$scores, 4))

  cat("\n--- Rankings (1 = best per metric) ---\n")
  print(x$rankings)

  cat("\n--- Recommendation ---\n")
  cat(x$recommendation, "\n")

  if (!is.null(x$bootstrap_ci)) {
    cat("\n--- Bootstrap 95% CIs ---\n")
    print(x$bootstrap_ci[, c("Model","Metric","Estimate","CI_lo","CI_hi",
                               "CIs_overlap"), drop = FALSE])
  }
  invisible(x)
}

#' Summary method for mvmetrics report
#' @param object Object of class `mvmetrics_report`.
#' @param ...    Ignored.
#' @export
summary.mvmetrics_report <- function(object, ...) {
  cat("\n=== mvmetrics Summary ===\n")
  cat("Models evaluated:", paste(rownames(object$scores), collapse = ", "), "\n")
  cat("Metrics computed:", paste(names(object$scores), collapse = ", "), "\n")
  cat("\nBest model per metric:\n")
  for (col in names(object$rankings)) {
    best <- rownames(object$rankings)[which.min(object$rankings[[col]])]
    cat(sprintf("  %-12s -> %s\n", col, best))
  }
  cat("\nRecommendation:\n", object$recommendation, "\n")
  invisible(object)
}
