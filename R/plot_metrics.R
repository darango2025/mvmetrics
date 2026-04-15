# ============================================================
# mvmetrics — Visualisation
# ============================================================

#' Plot joint metric scores and bootstrap CIs
#'
#' Generates diagnostic plots from an `mvmetrics_report` object. Requires
#' the \pkg{ggplot2} package (suggested dependency).
#'
#' @param report  Object of class `mvmetrics_report` from [evaluate_models()].
#' @param type    Character. One of:
#'   \describe{
#'     \item{"bars"}{Side-by-side bar chart of scores per metric family.}
#'     \item{"forest"}{Forest plot of bootstrap 95% CIs (requires
#'       `report$bootstrap_ci`).}
#'     \item{"heatmap"}{Heatmap of metric rankings.}
#'     \item{"decomp"}{Stacked bar chart of S_marg vs S_dep for Metrics D and E.}
#'   }
#' @param metrics Character vector. Which metrics to plot.
#'   Default: all available primary metrics (A, B, C, D, E).
#' @param ...     Additional arguments passed to ggplot2 theme functions.
#'
#' @return A `ggplot` object (invisible if ggplot2 is not installed).
#'
#' @examples
#' \dontrun{
#' # After running evaluate_models():
#' plot_metrics(report, type = "bars")
#' plot_metrics(report, type = "forest")   # requires bootstrap = TRUE
#' plot_metrics(report, type = "heatmap")
#' plot_metrics(report, type = "decomp")
#' }
#'
#' @export
plot_metrics <- function(report, type = "bars",
                          metrics = c("A","B","C","D","E"), ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_metrics(). ",
         "Install it with: install.packages('ggplot2')")

  scores   <- report$scores
  boot_ci  <- report$bootstrap_ci
  models   <- rownames(scores)

  # Colour palette consistent with thesis figures
  model_cols <- c("#c62828","#2e7d32","#ba7517","#1a3a5c",
                  "#534AB7","#0F6E56","#D84315","#888780")
  names(model_cols) <- models

  # ── BARS ──────────────────────────────────────────────────────────────────
  if (type == "bars") {
    avail_m <- intersect(metrics, names(scores))
    long <- do.call(rbind, lapply(models, function(m) {
      data.frame(Model  = m,
                 Metric = avail_m,
                 Score  = unlist(scores[m, avail_m]),
                 stringsAsFactors = FALSE)
    }))
    long$Metric <- factor(long$Metric, levels = avail_m)

    p <- ggplot2::ggplot(long,
           ggplot2::aes(x = Model, y = Score,
                        fill = Model, alpha = Metric)) +
      ggplot2::geom_col(position = "dodge", width = 0.7) +
      ggplot2::facet_wrap(~ Metric, scales = "free_y") +
      ggplot2::scale_fill_manual(values = model_cols) +
      ggplot2::scale_alpha_manual(values = rep(0.88, length(avail_m)),
                                   guide = "none") +
      ggplot2::labs(title = "Joint Metric Scores by Model",
                    subtitle = "Lower = better for all metrics",
                    x = "Model", y = "Score", fill = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom",
                     strip.text = ggplot2::element_text(face = "bold"))
    return(p)
  }

  # ── FOREST ────────────────────────────────────────────────────────────────
  if (type == "forest") {
    if (is.null(boot_ci))
      stop("Forest plot requires bootstrap CIs. Re-run evaluate_models() with bootstrap = TRUE.")

    plot_data <- boot_ci[boot_ci$Metric %in% metrics, ]
    plot_data$Metric <- factor(plot_data$Metric, levels = rev(metrics))

    p <- ggplot2::ggplot(plot_data,
           ggplot2::aes(x = Estimate, y = Metric, colour = Model,
                        xmin = CI_lo, xmax = CI_hi)) +
      ggplot2::geom_pointrange(
        position = ggplot2::position_dodge(width = 0.5),
        size = 0.7, linewidth = 0.9) +
      ggplot2::scale_colour_manual(values = model_cols) +
      ggplot2::labs(
        title    = "Bootstrap 95% Confidence Intervals",
        subtitle = "Non-overlapping CIs indicate statistically distinguishable performance",
        x = "Score (lower = better)", y = NULL, colour = "Model") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "top")

    if ("CIs_overlap" %in% names(plot_data)) {
      ann <- plot_data[plot_data$Model == models[1], ]
      p <- p + ggplot2::geom_text(
        data = ann,
        ggplot2::aes(x = CI_hi, y = Metric,
                     label = ifelse(!is.na(CIs_overlap) & !CIs_overlap,
                                    "no overlap", "")),
        hjust = -0.1, size = 2.8, colour = "#c62828",
        inherit.aes = FALSE)
    }
    return(p)
  }

  # ── HEATMAP ───────────────────────────────────────────────────────────────
  if (type == "heatmap") {
    rank_df <- report$rankings
    avail_r <- intersect(metrics, names(rank_df))
    long_r  <- do.call(rbind, lapply(models, function(m) {
      data.frame(Model  = m,
                 Metric = avail_r,
                 Rank   = unlist(rank_df[m, avail_r]),
                 stringsAsFactors = FALSE)
    }))
    long_r$Metric <- factor(long_r$Metric, levels = avail_r)
    long_r$Model  <- factor(long_r$Model,  levels = rev(models))

    p <- ggplot2::ggplot(long_r,
           ggplot2::aes(x = Metric, y = Model, fill = Rank)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = Rank), size = 4) +
      ggplot2::scale_fill_gradient(low = "#2e7d32", high = "#c62828",
                                    name = "Rank\n(1=best)") +
      ggplot2::labs(title = "Model Rankings by Metric Family",
                    x = "Metric", y = "Model") +
      ggplot2::theme_minimal(base_size = 11)
    return(p)
  }

  # ── DECOMP ────────────────────────────────────────────────────────────────
  if (type == "decomp") {
    d_cols <- intersect(c("D_marg","D_dep","E_marg","E_dep"), names(scores))
    if (length(d_cols) == 0)
      stop("No decomposition columns (D_marg, D_dep, E_marg, E_dep) found.")

    long_d <- do.call(rbind, lapply(models, function(m) {
      do.call(rbind, lapply(d_cols, function(col) {
        family    <- ifelse(grepl("^D", col), "D: MDD", "E: E-CVWMD")
        component <- ifelse(grepl("marg", col), "S_marg (marginal)",
                                                "S_dep (dependence)")
        data.frame(Model = m, Family = family, Component = component,
                   Score = scores[m, col], stringsAsFactors = FALSE)
      }))
    }))

    p <- ggplot2::ggplot(long_d,
           ggplot2::aes(x = Model, y = Score, fill = Component)) +
      ggplot2::geom_col(position = "dodge", width = 0.6, alpha = 0.88) +
      ggplot2::facet_wrap(~ Family, scales = "free_y") +
      ggplot2::scale_fill_manual(
        values = c("S_marg (marginal)" = "#1565C0",
                   "S_dep (dependence)" = "#D84315")) +
      ggplot2::labs(
        title    = "Marginal vs. Dependence Decomposition",
        subtitle = "Large S_dep for the joint model = genuine dependence recovery",
        x = "Model", y = "Score", fill = "Component") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    return(p)
  }

  stop("type must be one of: 'bars', 'forest', 'heatmap', 'decomp'.")
}
