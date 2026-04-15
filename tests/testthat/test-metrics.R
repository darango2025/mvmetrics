library(testthat)
library(mvmetrics)

# ── Reproducible test data ─────────────────────────────────────────────────
set.seed(42)
N <- 100; K <- 3; K_c <- 2; B <- 15

make_test_data <- function() {
  obs   <- cbind(matrix(rnorm(N * K_c, mean = c(20, 30), sd = c(2, 4)),
                         N, K_c),
                 rbinom(N, 1, 0.4))
  colnames(obs) <- c("X1","X2","Pbin")
  pred  <- cbind(obs[, 1:K_c] + rnorm(N * K_c, 0, 0.5), runif(N, 0.3, 0.6))
  colnames(pred) <- colnames(obs)
  draws <- array(rnorm(N * K * B), dim = c(N, K, B))
  list(obs = obs, pred = pred, draws = draws)
}

d <- make_test_data()

# ── normalise_matrix ─────────────────────────────────────────────────────────
test_that("normalise_matrix produces zero-mean unit-variance columns", {
  mu  <- colMeans(d$obs[, 1:K_c])
  sig <- apply(d$obs[, 1:K_c], 2, sd)
  Y_n <- normalise_matrix(d$obs[, 1:K_c], mu, sig)
  expect_equal(round(colMeans(Y_n), 8), c(X1 = 0, X2 = 0))
  expect_equal(round(apply(Y_n, 2, sd), 6), c(X1 = 1, X2 = 1))
})

test_that("normalise_matrix errors on dimension mismatch", {
  expect_error(normalise_matrix(d$obs[, 1:K_c], c(1), c(1, 2)))
})

# ── marginal_rmse ─────────────────────────────────────────────────────────────
test_that("marginal_rmse returns non-negative values", {
  r <- marginal_rmse(d$obs[, 1:K_c], d$pred[, 1:K_c])
  expect_true(all(r >= 0))
  expect_equal(length(r), K_c)
})

test_that("marginal_rmse is zero for perfect predictions", {
  r <- marginal_rmse(d$obs[, 1:K_c], d$obs[, 1:K_c])
  expect_equal(as.numeric(r), c(0, 0))
})

# ── marginal_logloss ──────────────────────────────────────────────────────────
test_that("marginal_logloss is positive", {
  ll <- marginal_logloss(d$obs[, K], d$pred[, K])
  expect_true(is.finite(ll) && ll > 0)
})

# ── marginal_accuracy ─────────────────────────────────────────────────────────
test_that("marginal_accuracy is in [0, 1]", {
  acc <- marginal_accuracy(d$obs[, K], d$pred[, K])
  expect_true(acc >= 0 && acc <= 1)
})

# ── metric_A ─────────────────────────────────────────────────────────────────
test_that("metric_A returns a list with finite total", {
  res <- metric_A(d$obs, d$pred, K_c = K_c)
  expect_true(is.list(res))
  expect_true(is.finite(res$total))
  expect_true(res$total >= 0)
})

# ── metric_B_energy_score ─────────────────────────────────────────────────────
test_that("metric_B_energy_score returns finite positive value", {
  obs_n <- normalise_matrix(d$obs[, 1:K_c],
                             colMeans(d$obs[, 1:K_c]),
                             apply(d$obs[, 1:K_c], 2, sd))
  # Pad binary column
  obs_full <- cbind(obs_n, d$obs[, K])
  drw_full <- d$draws
  es <- metric_B_energy_score(obs_full, drw_full)
  expect_true(is.finite(es) && es > 0)
})

# ── metric_D_mdd ─────────────────────────────────────────────────────────────
test_that("metric_D_mdd returns components with correct signs", {
  res <- metric_D_mdd(d$obs, d$pred, d$draws, K_c = K_c)
  expect_true(all(c("total","S_marg","S_dep","C_obs","C_model") %in% names(res)))
  expect_true(res$S_dep >= 0)
  expect_true(res$S_marg >= 0)
  expect_equal(dim(res$C_obs), c(K, K))
})

# ── metric_E_cvwmd ────────────────────────────────────────────────────────────
test_that("metric_E_cvwmd returns valid weights summing to 1", {
  res <- metric_E_cvwmd(d$obs, d$pred, d$draws,
                         mu_train = c(20, 30),
                         sd_train = c(2, 4),
                         K_c = K_c, run_dcor = FALSE)
  expect_true(is.finite(res$total))
  expect_equal(round(sum(res$w_cont) + res$w_bin, 6), 1)
  expect_true(res$S_E_dep >= 0)
})

test_that("metric_E_cvwmd fixed alpha overrides adaptive", {
  res <- metric_E_cvwmd(d$obs, d$pred, d$draws,
                         mu_train = c(20, 30),
                         sd_train = c(2, 4),
                         K_c = K_c, alpha = 0.3, run_dcor = FALSE)
  expect_equal(res$alpha_star, 0.3)
})

# ── evaluate_models ───────────────────────────────────────────────────────────
test_that("evaluate_models returns mvmetrics_report with correct structure", {
  report <- evaluate_models(
    obs        = d$obs,
    pred_list  = list(M1 = d$pred, M2 = d$pred),
    draws_list = list(M1 = d$draws, M2 = d$draws),
    mu_train   = c(20, 30),
    sd_train   = c(2, 4),
    K_c        = K_c,
    verbose    = FALSE
  )
  expect_s3_class(report, "mvmetrics_report")
  expect_true("A" %in% names(report$scores))
  expect_true("B" %in% names(report$scores))
  expect_equal(nrow(report$scores), 2)
  expect_true(is.character(report$recommendation))
  expect_equal(nrow(report$rankings), 2)
})

test_that("evaluate_models works without draws (Metric A only)", {
  report <- evaluate_models(
    obs       = d$obs,
    pred_list = list(M1 = d$pred),
    K_c       = K_c,
    verbose   = FALSE
  )
  expect_true("A" %in% names(report$scores))
  expect_false("B" %in% names(report$scores))
})

# ── recommend_metric ──────────────────────────────────────────────────────────
test_that("recommend_metric returns character string", {
  scores_df <- data.frame(A = c(1.5, 1.6), B = c(1.2, 0.9),
                           D = c(2.0, 1.0), D_marg = c(1.5, 1.5),
                           D_dep = c(0.5, 0.5), row.names = c("M1","M2"))
  rec <- recommend_metric(scores_df, n_models = 2,
                           has_draws = TRUE, has_E = FALSE)
  expect_true(is.character(rec))
  expect_true(nchar(rec) > 10)
})
