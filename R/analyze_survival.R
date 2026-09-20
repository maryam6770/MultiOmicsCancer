#' Survival Analysis for Pan-Cancer Studies
#'
#' Performs survival analysis using clinical data and a molecular feature
#' (gene expression, methylation, or a risk score). Supports Kaplan-Meier
#' curves, log-rank tests, and univariate Cox proportional hazards models.
#'
#' @param mae A \code{MultiAssayExperiment} object created by
#'   \code{\link{load_cancer_data}}, containing clinical survival data.
#' @param time_column Character. Name of the clinical column containing
#'   follow-up time (in days or months). Default \code{"OS.time"}.
#' @param event_column Character. Name of the clinical column containing
#'   event status (1 = event/death, 0 = censored). Default \code{"OS"}.
#' @param feature_name Character. Name of the molecular feature to test
#'   (gene symbol or probe ID). If \code{NULL}, only clinical summary is
#'   returned.
#' @param assay_name Character. Which assay to use for the feature. One of
#'   \code{"expression"} (default), \code{"methylation"}, or \code{"mirna"}.
#' @param split_method Character. How to split samples into high/low groups.
#'   One of \code{"median"} (default), \code{"tertile"}, or \code{"quartile"}.
#' @param verbose Logical. Whether to print progress messages. Default
#'   \code{TRUE}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{km_fit}{A \code{survfit} object for Kaplan-Meier curves}
#'   \item{logrank_p}{P-value from the log-rank test}
#'   \item{cox_model}{A \code{coxph} object for univariate Cox regression}
#'   \item{cox_summary}{Summary table of the Cox model}
#'   \item{data}{Data frame used for the analysis}
#'   \item{group}{Factor of high/low groups}
#' }
#'
#' @examples
#' \dontrun{
#' # Load data
#' brca <- load_cancer_data("brca_tcga")
#'
#' # Survival analysis for TP53 expression
#' surv <- analyze_survival(brca, feature_name = "TP53")
#'
#' # Survival analysis for a methylation probe
#' surv_meth <- analyze_survival(brca,
#'                               feature_name = "cg00000029",
#'                               assay_name = "methylation")
#' }
#'
#' @export
analyze_survival <- function(mae,
                             time_column = "OS.time",
                             event_column = "OS",
                             feature_name = NULL,
                             assay_name = "expression",
                             split_method = "median",
                             verbose = TRUE) {

  # --- Validate inputs ---
  if (!methods::is(mae, "MultiAssayExperiment")) {
    stop("'mae' must be a MultiAssayExperiment object.")
  }

  if (!split_method %in% c("median", "tertile", "quartile")) {
    stop("'split_method' must be one of: 'median', 'tertile', 'quartile'.")
  }

  # --- Extract clinical data ---
  clinical <- as.data.frame(MultiAssayExperiment::colData(mae))

  if (!time_column %in% colnames(clinical)) {
    stop("Time column '", time_column, "' not found in clinical data. ",
         "Available columns: ", paste(colnames(clinical), collapse = ", "))
  }

  if (!event_column %in% colnames(clinical)) {
    stop("Event column '", event_column, "' not found in clinical data. ",
         "Available columns: ", paste(colnames(clinical), collapse = ", "))
  }

  # --- Clean survival data ---
  surv_data <- data.frame(
    sample = rownames(clinical),
    time = as.numeric(clinical[[time_column]]),
    event = as.numeric(clinical[[event_column]]),
    stringsAsFactors = FALSE
  )

  # Remove samples with missing or invalid survival data
  surv_data <- surv_data[!is.na(surv_data$time) & !is.na(surv_data$event) &
                           surv_data$time > 0, ]

  if (nrow(surv_data) < 10) {
    stop("Not enough samples with valid survival data. Found ",
         nrow(surv_data), " samples.")
  }

  if (verbose) message("Using ", nrow(surv_data), " samples with valid survival data.")

  # --- If no feature specified, return clinical summary only ---
  if (is.null(feature_name)) {
    fit <- survival::survfit(survival::Surv(time, event) ~ 1, data = surv_data)
    return(list(
      km_fit = fit,
      data = surv_data,
      message = "No feature specified. Returning overall survival summary."
    ))
  }

  # --- Extract feature data ---
  if (!assay_name %in% names(MultiAssayExperiment::experiments(mae))) {
    stop("Assay '", assay_name, "' not found. ",
         "Available assays: ",
         paste(names(MultiAssayExperiment::experiments(mae)), collapse = ", "))
  }

  feature_se <- MultiAssayExperiment::experiments(mae)[[assay_name]]
  feature_data <- SummarizedExperiment::assay(feature_se)

  if (!feature_name %in% rownames(feature_data)) {
    stop("Feature '", feature_name, "' not found in assay '", assay_name, "'. ",
         "Check the feature name or use a different assay.")
  }

  # --- Match samples ---
  common_samples <- intersect(surv_data$sample, colnames(feature_data))

  if (length(common_samples) < 10) {
    stop("Not enough common samples between survival and feature data. ",
         "Found ", length(common_samples), " samples.")
  }

  surv_data <- surv_data[surv_data$sample %in% common_samples, ]
  feature_vec <- as.numeric(feature_data[feature_name, surv_data$sample])

  surv_data$feature_value <- feature_vec

  # --- Split into high/low groups ---
  if (split_method == "median") {
    cutoff <- stats::median(surv_data$feature_value, na.rm = TRUE)
    surv_data$group <- ifelse(surv_data$feature_value > cutoff, "High", "Low")
  } else if (split_method == "tertile") {
    cutoffs <- stats::quantile(surv_data$feature_value,
                               probs = c(1/3, 2/3), na.rm = TRUE)
    surv_data$group <- ifelse(surv_data$feature_value > cutoffs[2], "High",
                              ifelse(surv_data$feature_value < cutoffs[1], "Low", "Middle"))
    surv_data <- surv_data[surv_data$group != "Middle", ]
  } else if (split_method == "quartile") {
    cutoffs <- stats::quantile(surv_data$feature_value,
                               probs = c(0.25, 0.75), na.rm = TRUE)
    surv_data$group <- ifelse(surv_data$feature_value > cutoffs[2], "High",
                              ifelse(surv_data$feature_value < cutoffs[1], "Low", "Middle"))
    surv_data <- surv_data[surv_data$group != "Middle", ]
  }

  surv_data$group <- factor(surv_data$group, levels = c("Low", "High"))

  if (verbose) message("Split into ", sum(surv_data$group == "High"), " High and ",
                       sum(surv_data$group == "Low"), " Low samples.")

  # --- Kaplan-Meier fit ---
  km_fit <- survival::survfit(survival::Surv(time, event) ~ group, data = surv_data)

  # --- Log-rank test ---
  logrank <- survival::survdiff(survival::Surv(time, event) ~ group, data = surv_data)
  logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1,
                             lower.tail = FALSE)

  # --- Cox proportional hazards model ---
  cox_model <- survival::coxph(survival::Surv(time, event) ~ group, data = surv_data)
  cox_summary <- summary(cox_model)

  if (verbose) {
    hr <- cox_summary$coefficients[1, "exp(coef)"]
    message("Done. Hazard Ratio (High vs Low) = ", round(hr, 3),
            ", Log-rank p = ", format.pval(logrank_p, digits = 3))
  }

  return(list(
    km_fit = km_fit,
    logrank_p = logrank_p,
    cox_model = cox_model,
    cox_summary = cox_summary,
    data = surv_data,
    group = surv_data$group,
    feature = feature_name,
    assay = assay_name
  ))
}
