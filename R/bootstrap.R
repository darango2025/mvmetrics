# ============================================================
# mvmetrics — Bootstrap uncertainty quantification
# ============================================================

#' Non-parametric bootstrap for joint metric confidence intervals
#'
#' Resamples test observations with replacement and recomputes all
#' available metrics (A, B, D, E) for each resample. Returns 95%
#' percentile CIs and a pairwise overlap test between models.
#'
#' @param obs        Numeric matrix \[N x K\]. Observations.
#' @param pred_list  Named list of prediction matrices.
#' @param draws_list Named list of draw arrays.
#' @param mu_train   Numeric vector \[K_c\] (for Metric E normalisation).
#' @param sd_train   Numeric vector \[K_c\].
#' @param K_c        Integer. Number of continuous variables.
#' @param alpha_D    Numeric. Alpha for Metric D (default 0.5).
#' @param B_boot     Integer. Number of bootstrap resamples (default 200).
#' @param verbose    Logical. Progress messages (default TRUE).
#'
#' @return data.frame with columns:
#'   Model, Metric, Estimate, SE, CI_lo, CI_hi, CIs_overlap (pairwise M1 vs M2,
#'   only when exactly two models are provided).
#'
#' @export
bootstrap_metrics <- function(obs, pred_list, draws_list,
                               mu_train = NULL, sd_train = NULL,
                               K_c = NULL, alpha_D = 0.5,
                               B_boot = 200L, verbose = TRUE) {
  obs  <- as.matrix(obs)
  K    <- ncol(obs); N <- nrow(obs)
  if (is.null(K_c)) K_c <- K - 1L
  model_names <- names(pred_list)

  mu_norm <- if (!is.null(mu_train)) mu_train else
               colMeans(obs[, seq_len(K_c), drop = FALSE], na.rm = TRUE)
  sd_norm <- if (!is.null(sd_train)) sd_train else
               apply(obs[, seq_len(K_c), drop = FALSE], 2, stats::sd, na.rm = TRUE)

  has_E <- !is.null(mu_train) && !is.null(sd_train)

  # Metric keys to collect
  metric_keys <- c("A","B","D","D_dep")
  if (has_E) metric_keys <- c(metric_keys, "E","E_dep")

  # Storage: list of matrices [B_boot x n_metrics] per model
  boot_store <- lapply(model_names, function(m)
    matrix(NA_real_, B_boot, length(metric_keys),
           dimnames = list(NULL, metric_keys)))
  names(boot_store) <- model_names

  for (b in seq_len(B_boot)) {
    idx <- sample.int(N, N, replace = TRUE)
    obs_b <- obs[idx, , drop = FALSE]

    for (m in model_names) {
      y_pred <- as.matrix(pred_list[[m]])[idx, , drop = FALSE]
      y_drw  <- draws_list[[m]][idx, , , drop = FALSE]

      y_pred_n <- y_pred
      y_pred_n[, seq_len(K_c)] <- normalise_matrix(
        y_pred[, seq_len(K_c), drop = FALSE], mu_norm, sd_norm)
      obs_n <- obs_b
      obs_n[, seq_len(K_c)] <- normalise_matrix(
        obs_b[, seq_len(K_c), drop = FALSE], mu_norm, sd_norm)
      y_drw_n <- y_drw
      for (k in seq_len(K_c))
        y_drw_n[, k, ] <- (y_drw[, k, ] - mu_norm[k]) / pmax(sd_norm[k], 1e-10)

      sc <- tryCatch({
        rA <- metric_A(obs_b, y_pred, K_c = K_c)$total
        rB <- metric_B_energy_score(obs_n, y_drw_n)
        rD <- metric_D_mdd(obs_b, y_pred, y_drw, K_c = K_c, alpha = alpha_D)
        row <- c(A = rA, B = rB, D = rD$total, D_dep = rD$S_dep)
        if (has_E) {
          rE  <- metric_E_cvwmd(obs_b, y_pred, y_drw,
                                 mu_train = mu_train, sd_train = sd_train,
                                 K_c = K_c, alpha = 0.5, run_dcor = FALSE)
          row <- c(row, E = rE$total, E_dep = rE$S_E_dep)
        }
        row
      }, error = function(e) rep(NA_real_, length(metric_keys)))

      boot_store[[m]][b, names(sc)] <- sc
    }

    if (verbose && b %% 50 == 0)
      message("    Bootstrap resample ", b, " / ", B_boot)
  }

  # Summarise: point estimates + percentile CIs
  point_est <- lapply(model_names, function(m) {
    obs_n <- obs
    obs_n[, seq_len(K_c)] <- normalise_matrix(
      obs[, seq_len(K_c), drop = FALSE], mu_norm, sd_norm)
    y_pred <- as.matrix(pred_list[[m]])
    y_drw  <- draws_list[[m]]
    y_pred_n <- y_pred
    y_pred_n[, seq_len(K_c)] <- normalise_matrix(
      y_pred[, seq_len(K_c), drop = FALSE], mu_norm, sd_norm)
    y_drw_n <- y_drw
    for (k in seq_len(K_c))
      y_drw_n[, k, ] <- (y_drw[, k, ] - mu_norm[k]) / pmax(sd_norm[k], 1e-10)

    rA <- metric_A(obs, y_pred, K_c = K_c)$total
    rB <- metric_B_energy_score(obs_n, y_drw_n)
    rD <- metric_D_mdd(obs, y_pred, y_drw, K_c = K_c, alpha = alpha_D)
    pt <- c(A = rA, B = rB, D = rD$total, D_dep = rD$S_dep)
    if (has_E) {
      rE  <- metric_E_cvwmd(obs, y_pred, y_drw,
                             mu_train = mu_train, sd_train = sd_train,
                             K_c = K_c, alpha = 0.5, run_dcor = FALSE)
      pt  <- c(pt, E = rE$total, E_dep = rE$S_E_dep)
    }
    pt
  })
  names(point_est) <- model_names

  ci_rows <- do.call(rbind, lapply(model_names, function(m) {
    mat <- boot_store[[m]]
    pt  <- point_est[[m]]
    do.call(rbind, lapply(metric_keys, function(nm) {
      v  <- mat[, nm]; v <- v[is.finite(v)]
      data.frame(
        Model    = m,
        Metric   = nm,
        Estimate = round(pt[nm], 4),
        SE       = round(stats::sd(v), 4),
        CI_lo    = round(stats::quantile(v, 0.025), 4),
        CI_hi    = round(stats::quantile(v, 0.975), 4),
        stringsAsFactors = FALSE
      )
    }))
  }))

  # Pairwise overlap test (only for exactly 2 models)
  ci_rows$CIs_overlap <- NA
  if (length(model_names) == 2) {
    wide <- stats::reshape(ci_rows[, c("Model","Metric","CI_lo","CI_hi")],
                           idvar = "Metric", timevar = "Model",
                           direction = "wide")
    m1  <- model_names[1]; m2 <- model_names[2]
    lo1 <- wide[[paste0("CI_lo.", m1)]]; hi1 <- wide[[paste0("CI_hi.", m1)]]
    lo2 <- wide[[paste0("CI_lo.", m2)]]; hi2 <- wide[[paste0("CI_hi.", m2)]]
    overlap <- !(hi1 < lo2 | hi2 < lo1)
    for (i in seq_len(nrow(ci_rows))) {
      met <- ci_rows$Metric[i]
      idx <- which(wide$Metric == met)
      if (length(idx) == 1) ci_rows$CIs_overlap[i] <- overlap[idx]
    }
  }

  rownames(ci_rows) <- NULL
  ci_rows
}
