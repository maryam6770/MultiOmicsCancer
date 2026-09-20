#' Visualization Functions for Pan-Cancer Multi-Omics Analysis
#'
#' A collection of publication-quality plotting functions for visualizing
#' results from differential methylation, differential expression, multi-omics
#' integration, and survival analysis. All plots are based on \code{ggplot2}
#' and return ggplot objects that can be further customized.
#'
#' @param data A data frame returned by \code{\link{analyze_methylation}},
#'   \code{\link{analyze_expression}}, or \code{\link{integrate_meth_expr}}.
#' @param plot_type Character. Type of plot to generate. One of
#'   \code{"volcano"}, \code{"heatmap"}, \code{"correlation"}, or
#'   \code{"survival"}.
#' @param title Optional. Plot title. If \code{NULL}, a default title is used.
#' @param p_cutoff Numeric. Adjusted p-value threshold for highlighting
#'   significant features. Default \code{0.05}.
#' @param fc_cutoff Numeric. Log fold change or delta-beta threshold for
#'   highlighting significant features. Default \code{1} for expression and
#'   \code{0.1} for methylation. If \code{NULL}, uses the \code{significant}
#'   column in the data.
#' @param color_up Character. Color for upregulated/hypermethylated features.
#'   Default \code{"firebrick"}.
#' @param color_down Character. Color for downregulated/hypomethylated features.
#'   Default \code{"steelblue"}.
#' @param color_ns Character. Color for non-significant features.
#'   Default \code{"grey70"}.
#' @param max_labels Integer. Maximum number of labels to display on volcano
#'   plots. Default \code{10}.
#' @param surv_result Optional. A list returned by \code{\link{analyze_survival}}
#'   for survival plots.
#' @param verbose Logical. Whether to print progress messages. Default
#'   \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # Volcano plot of differential expression
#' deg <- analyze_expression(brca)
#' plot_cancer(deg, plot_type = "volcano")
#'
#' # Volcano plot of differential methylation
#' dmr <- analyze_methylation(brca)
#' plot_cancer(dmr, plot_type = "volcano")
#'
#' # Survival plot
#' surv <- analyze_survival(brca, feature_name = "TP53")
#' plot_cancer(plot_type = "survival", surv_result = surv)
#' }
#'
#' @importFrom rlang .data
#' @importFrom utils head
#' @export
plot_cancer <- function(data = NULL,
                        plot_type = c("volcano", "heatmap", "correlation", "survival"),
                        title = NULL,
                        p_cutoff = 0.05,
                        fc_cutoff = NULL,
                        color_up = "firebrick",
                        color_down = "steelblue",
                        color_ns = "grey70",
                        max_labels = 10,
                        surv_result = NULL,
                        verbose = TRUE) {

  plot_type <- match.arg(plot_type)

  if (plot_type == "survival") {
    return(plot_survival(surv_result, title = title))
  }

  if (is.null(data)) {
    stop("'data' must be provided for plot_type = '", plot_type, "'.")
  }

  if (plot_type == "volcano") {
    return(plot_volcano(data, title = title, p_cutoff = p_cutoff,
                        fc_cutoff = fc_cutoff, color_up = color_up,
                        color_down = color_down, color_ns = color_ns,
                        max_labels = max_labels, verbose = verbose))
  }

  if (plot_type == "heatmap") {
    stop("Heatmap plotting requires a MultiAssayExperiment. ",
         "Use plot_heatmap() directly (coming in a future version).")
  }

  if (plot_type == "correlation") {
    return(plot_correlation(data, title = title, p_cutoff = p_cutoff,
                            color_up = color_up, color_down = color_down,
                            color_ns = color_ns, verbose = verbose))
  }
}


#' @keywords internal
#' @noRd
plot_volcano <- function(data, title, p_cutoff, fc_cutoff,
                         color_up, color_down, color_ns, max_labels, verbose) {

  # Detect whether this is expression or methylation data
  is_methylation <- "delta_beta" %in% colnames(data)

  # Determine effect-size column and threshold
  if (is_methylation) {
    effect_col <- "delta_beta"
    if (is.null(fc_cutoff)) fc_cutoff <- 0.1
    x_label <- "Delta Beta"
  } else {
    effect_col <- "logFC"
    if (is.null(fc_cutoff)) fc_cutoff <- 1
    x_label <- "Log2 Fold Change"
  }

  # Check required columns
  if (!all(c(effect_col, "adj.P.Val") %in% colnames(data))) {
    stop("'data' must contain columns '", effect_col, "' and 'adj.P.Val'.")
  }

  # Determine label column
  label_col <- if ("gene" %in% colnames(data)) "gene" else
    if ("probe" %in% colnames(data)) "probe" else NULL

  # Classify features
  data$category <- "Not significant"
  data$category[data$adj.P.Val < p_cutoff & data[[effect_col]] >= fc_cutoff] <-
    if (is_methylation) "Hypermethylated" else "Upregulated"
  data$category[data$adj.P.Val < p_cutoff & data[[effect_col]] <= -fc_cutoff] <-
    if (is_methylation) "Hypomethylated" else "Downregulated"
  data$category <- factor(data$category,
                          levels = c("Upregulated", "Hypermethylated",
                                     "Downregulated", "Hypomethylated",
                                     "Not significant"))

  # Colors
  color_map <- c(
    "Upregulated"     = color_up,
    "Hypermethylated" = color_up,
    "Downregulated"   = color_down,
    "Hypomethylated"  = color_down,
    "Not significant" = color_ns
  )

  # Compute -log10(p)
  data$neg_log10_p <- -log10(data$adj.P.Val)

  # Build plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[effect_col]],
                                          y = .data$neg_log10_p,
                                          color = .data$category)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = color_map, name = "") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff),
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                        linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      title = if (is.null(title)) paste("Volcano Plot -", if (is_methylation) "Methylation" else "Expression") else title,
      x = x_label,
      y = expression(-log[10]("Adjusted P-value"))
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  # Add labels for top features
  if (!is.null(label_col) && max_labels > 0) {
    top_features <- data[data$category != "Not significant", ]
    top_features <- top_features[order(top_features$adj.P.Val), ]
    top_features <- utils::head(top_features, max_labels)

    if (nrow(top_features) > 0) {
      p <- p + ggplot2::geom_text(
        data = top_features,
        ggplot2::aes(label = .data[[label_col]]),
        vjust = -0.8, size = 3, check_overlap = TRUE
      )
    }
  }

  if (verbose) message("Volcano plot generated.")
  return(p)
}


#' @keywords internal
#' @noRd
plot_correlation <- function(data, title, p_cutoff,
                             color_up, color_down, color_ns, verbose) {

  if (!all(c("correlation", "adj_p_value") %in% colnames(data))) {
    stop("'data' must contain columns 'correlation' and 'adj_p_value'. ",
         "Use output from integrate_meth_expr().")
  }

  # Classify
  data$category <- "Not significant"
  data$category[data$adj_p_value < p_cutoff & data$correlation < 0] <- "Negative"
  data$category[data$adj_p_value < p_cutoff & data$correlation > 0] <- "Positive"
  data$category <- factor(data$category,
                          levels = c("Negative", "Positive", "Not significant"))

  color_map <- c(
    "Negative"        = color_down,
    "Positive"        = color_up,
    "Not significant" = color_ns
  )

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$correlation,
                                          fill = .data$category)) +
    ggplot2::geom_histogram(bins = 40, alpha = 0.8, color = "white") +
    ggplot2::scale_fill_manual(values = color_map, name = "") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      title = if (is.null(title)) "Methylation-Expression Correlation" else title,
      x = "Correlation Coefficient",
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  if (verbose) message("Correlation plot generated.")
  return(p)
}


#' @keywords internal
#' @noRd
plot_survival <- function(surv_result, title) {

  if (is.null(surv_result)) {
    stop("'surv_result' must be provided for survival plots. ",
         "Use output from analyze_survival().")
  }

  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is required for survival plots. ",
         "Install with: install.packages('survminer')")
  }

  if (is.null(surv_result$km_fit)) {
    stop("'surv_result' does not contain a Kaplan-Meier fit. ",
         "Make sure you specified a feature_name in analyze_survival().")
  }

  p <- survminer::ggsurvplot(
    surv_result$km_fit,
    data = surv_result$data,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    palette = c("steelblue", "firebrick"),
    legend.title = "Group",
    legend.labs = c("Low", "High"),
    title = if (is.null(title)) paste("Survival by", surv_result$feature) else title,
    xlab = "Time",
    ylab = "Survival Probability",
    ggtheme = ggplot2::theme_minimal(base_size = 12)
  )

  return(p)
}
