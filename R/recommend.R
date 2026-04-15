# ============================================================
# mvmetrics — Metric recommendation engine
# ============================================================

#' Recommend the most appropriate metric family for the data
#'
#' Analyses the pattern of scores and rankings across metric families and
#' returns a character string explaining which metric family is most
#' informative for the specific models and data provided.
#'
#' The recommendation logic follows the decision tree documented in
#' Arango Londoño (2026), Chapter 9:
#' \enumerate{
#'   \item If all models are tied on Metric A but separated on B/C/D/E
#'         → dependence metrics are informative.
#'   \item If S_dep >> S_marg for the best model → use Metric D or E.
#'   \item If scale heterogeneity is detected (Metric A explodes for some models)
#'         → use Metric B (Energy Score) as primary.
#'   \item If only point predictions are available (no draws) → Metric A only.
#' }
#'
#' @param scores_df  data.frame of metric scores as returned by
#'   [evaluate_models()].
#' @param n_models   Integer. Number of models evaluated.
#' @param has_draws  Logical. Are predictive draw arrays available?
#' @param has_E      Logical. Is Metric E (E-CVWMD) available?
#'
#' @return Character string with the recommendation and its rationale.
#'
#' @export
recommend_metric <- function(scores_df, n_models, has_draws, has_E) {

  if (!has_draws) {
    return(paste(
      "RECOMMENDED: Metric A (Weighted Marginal Sum).\n",
      "Reason: No predictive draws provided. Metrics B, C, D, E require",
      "a draw array [N x K x B]. Provide draws_list to unlock joint metrics."
    ))
  }

  cols <- names(scores_df)

  # --- Detect Metric A range ---
  if ("A" %in% cols) {
    A_vals  <- scores_df[["A"]]
    A_range <- diff(range(A_vals, na.rm = TRUE))
    A_max   <- max(A_vals, na.rm = TRUE)
    A_min   <- min(A_vals, na.rm = TRUE)
    # Scale explosion: max > 5 * min (rough heuristic)
    scale_explosion <- is.finite(A_max) && is.finite(A_min) &&
                       A_min > 0 && (A_max / A_min) > 5
  } else {
    A_range <- 0; scale_explosion <- FALSE
  }

  # --- Metric A flatness: all within 5% of each other ---
  A_flat <- ("A" %in% cols) && (A_range / (mean(scores_df[["A"]], na.rm=TRUE) + 1e-10)) < 0.05

  # --- Dependence signal: is S_dep large relative to S_marg? ---
  dep_signal <- FALSE
  if (all(c("D_dep","D_marg") %in% cols)) {
    best_dep  <- min(scores_df[["D_dep"]],  na.rm = TRUE)
    best_marg <- min(scores_df[["D_marg"]], na.rm = TRUE)
    dep_signal <- is.finite(best_dep) && is.finite(best_marg) &&
                  best_dep > 0.1 * best_marg
  }

  # --- Build recommendation ---
  if (scale_explosion) {
    msg <- paste0(
      "RECOMMENDED: Metric B — Multivariate Energy Score.\n",
      "Reason: Scale heterogeneity detected — Metric A varies by a factor of ",
      round(A_max / A_min, 1), "x across models, indicating that the marginal ",
      "sum is dominated by the highest-variance variable and its rankings are ",
      "an artefact of scale rather than predictive quality. Metric B normalises ",
      "internally and is robust to this pathology."
    )
    if (has_E)
      msg <- paste0(msg, "\nNote: Metric E (E-CVWMD) also explodes under scale ",
                    "heterogeneity because CV weights amplify the dominant variable. ",
                    "Apply normalise_matrix() before computing Metric A or E.")
    return(msg)
  }

  if (A_flat && dep_signal) {
    msg <- paste0(
      "RECOMMENDED: Metric D — Marginal-Dependence Decomposition (MDD).\n",
      "Reason: Metric A is essentially flat across models (range = ",
      round(A_range, 4), "), indicating that marginal scoring rules cannot ",
      "distinguish them. The S_dep component of Metric D shows meaningful ",
      "separation, suggesting that the performance difference lies in how well ",
      "each model recovers the cross-variable dependence structure.\n",
      "Also report: Metric B (Energy Score) as a strictly proper cross-check."
    )
    if (has_E) {
      e_dep_vals <- if ("E_dep" %in% cols) scores_df[["E_dep"]] else NULL
      if (!is.null(e_dep_vals) && any(!is.na(e_dep_vals))) {
        msg <- paste0(msg, "\n",
          "Metric E (E-CVWMD) also shows dependence separation and provides ",
          "an independent corroboration via residual correlations (no CDF ",
          "specification required). Include E alongside D for robustness.")
      }
    }
    return(msg)
  }

  if (!A_flat && !dep_signal) {
    return(paste0(
      "RECOMMENDED: Metric A + Metric B.\n",
      "Reason: Models differ on both marginal accuracy (Metric A range = ",
      round(A_range, 4), ") and joint calibration (Metric B). Report both to ",
      "characterise whether the performance advantage is marginal, joint, or both.\n",
      "Use Metric D's S_marg / S_dep decomposition to attribute the source of ",
      "improvement before selecting a single primary metric."
    ))
  }

  # Default
  primary <- if (has_E) "Metric E (E-CVWMD)" else "Metric D (MDD)"
  paste0(
    "RECOMMENDED: ", primary, " + Metric B (Energy Score).\n",
    "Reason: Models show both marginal and dependence differences. ",
    primary, " provides a comprehensive summary combining CV-weighted ",
    "marginal accuracy and residual-correlation dependence recovery. ",
    "Metric B serves as a strictly proper cross-check. Report Metric A ",
    "as the marginal baseline for comparison with existing literature."
  )
}
