#' Differential Methylation Analysis for Pan-Cancer Studies
#'
#' Performs differential methylation analysis between tumor and normal samples
#' (or any two groups defined in the clinical data) for any cancer type. Uses
#' the \code{limma} framework for linear modeling of methylation beta values.
#'
#' @param mae A \code{MultiAssayExperiment} object created by
#'   \code{\link{load_cancer_data}}.
#' @param assay_name Character. Name of the methylation assay in the
#'   \code{MultiAssayExperiment}. Default \code{"methylation"}.
#' @param group_column Character. Name of the clinical column used to define
#'   groups. Default \code{"sample_type"}.
#' @param group_levels Character vector of length 2. The two groups to compare.
#'   The first is the reference (e.g., normal), the second is the comparison
#'   (e.g., tumor). If \code{NULL}, the function attempts to auto-detect
#'   common tumor/normal labels.
#' @param method Character. Differential methylation method. Currently supports
#'   \code{"limma"} (default).
#' @param adj_method Character. Multiple testing correction method passed to
#'   \code{limma::topTable}. Default \code{"BH"} (Benjamini-Hochberg).
#' @param p_cutoff Numeric. Adjusted p-value threshold for significance.
#'   Default \code{0.05}. Lower values increase specificity; higher values
#'   increase sensitivity.
#' @param delta_beta Numeric. Minimum absolute difference in mean beta value
#'   between groups to call a probe differentially methylated. Default
#'   \code{0.1} (10\%). Lower values (e.g., \code{0.05}) return more probes;
#'   higher values (e.g., \code{0.2}) return fewer, more robust probes.
#' @param verbose Logical. Whether to print progress messages. Default
#'   \code{TRUE}.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{probe}{Probe or CpG identifier}
#'   \item{logFC}{Difference in mean methylation (group2 - group1)}
#'   \item{AveExpr}{Average methylation across all samples}
#'   \item{t}{Moderated t-statistic}
#'   \item{P.Value}{Raw p-value}
#'   \item{adj.P.Val}{Adjusted p-value}
#'   \item{delta_beta}{Absolute difference in mean beta between groups}
#'   \item{direction}{"Hyper" or "Hypo" methylation in group2 vs group1}
#'   \item{significant}{Logical, whether the probe passes both p-value and
#'     delta-beta cutoffs}
#' }
#'
#' @examples
#' \dontrun{
#' # Load data
#' brca <- load_cancer_data("brca_tcga")
#'
#' # Differential methylation: normal vs tumor
#' dmr <- analyze_methylation(
#'   brca,
#'   group_column = "sample_type",
#'   group_levels = c("Solid Tissue Normal", "Primary Tumor")
#' )
#'
#' # Filter significant probes
#' sig_dmr <- dmr[dmr$significant, ]
#'
#' # More sensitive analysis with lower thresholds
#' dmr_sensitive <- analyze_methylation(brca, delta_beta = 0.05,
#'                                       p_cutoff = 0.01)
#' }
#'
#' @export
analyze_methylation <- function(mae,
                                assay_name = "methylation",
                                group_column = "sample_type",
                                group_levels = NULL,
                                method = "limma",
                                adj_method = "BH",
                                p_cutoff = 0.05,
                                delta_beta = 0.1,
                                verbose = TRUE) {

  # --- Validate inputs ---
  if (!methods::is(mae, "MultiAssayExperiment")) {
    stop("'mae' must be a MultiAssayExperiment object.")
  }

  if (!assay_name %in% names(MultiAssayExperiment::experiments(mae))) {
    stop("Assay '", assay_name, "' not found in the MultiAssayExperiment. ",
         "Available assays: ",
         paste(names(MultiAssayExperiment::experiments(mae)), collapse = ", "))
  }

  if (!method %in% c("limma")) {
    stop("'method' must be 'limma'. Other methods will be added in future versions.")
  }

  # --- Extract data ---
  meth_se <- MultiAssayExperiment::experiments(mae)[[assay_name]]
  meth_data <- SummarizedExperiment::assay(meth_se)
  clinical <- as.data.frame(MultiAssayExperiment::colData(mae))

  if (!group_column %in% colnames(clinical)) {
    stop("Column '", group_column, "' not found in clinical data. ",
         "Available columns: ", paste(colnames(clinical), collapse = ", "))
  }

  # --- Determine groups ---
  if (is.null(group_levels)) {
    # Attempt automatic detection of tumor and normal labels
    unique_vals <- unique(as.character(clinical[[group_column]]))
    tumor_patterns <- c("tumor", "tumour", "primary", "cancer", "malignant")
    normal_patterns <- c("normal", "solid tissue", "healthy", "control", "adjacent")

    tumor_val <- unique_vals[grepl(paste(tumor_patterns, collapse = "|"),
                                   unique_vals, ignore.case = TRUE)][1]
    normal_val <- unique_vals[grepl(paste(normal_patterns, collapse = "|"),
                                    unique_vals, ignore.case = TRUE)][1]

    if (is.na(tumor_val) || is.na(normal_val)) {
      stop("Could not auto-detect group levels. Please specify 'group_levels' ",
           "manually. Unique values: ", paste(unique_vals, collapse = ", "))
    }
    group_levels <- c(normal_val, tumor_val)
    if (verbose) message("Auto-detected groups: reference = '", normal_val,
                         "', comparison = '", tumor_val, "'")
  }

  if (length(group_levels) != 2) {
    stop("'group_levels' must be a character vector of length 2.")
  }

  # --- Filter samples ---
  keep_samples <- clinical[[group_column]] %in% group_levels
  if (sum(keep_samples) < 4) {
    stop("Not enough samples in the specified groups. Found ",
         sum(keep_samples), " samples.")
  }

  meth_data <- meth_data[, keep_samples, drop = FALSE]
  clinical <- clinical[keep_samples, , drop = FALSE]

  # --- Ensure sample order matches ---
  if (!identical(colnames(meth_data), rownames(clinical))) {
    common <- intersect(colnames(meth_data), rownames(clinical))
    if (length(common) == 0) {
      stop("Sample names in methylation data and clinical data do not match.")
    }
    meth_data <- meth_data[, common, drop = FALSE]
    clinical <- clinical[common, , drop = FALSE]
  }

  # --- Remove zero-variance probes ---
  probe_var <- apply(meth_data, 1, stats::var, na.rm = TRUE)
  keep_probes <- !is.na(probe_var) & probe_var > 0
  meth_data <- meth_data[keep_probes, , drop = FALSE]

  if (nrow(meth_data) < 10) {
    stop("Too few informative probes after filtering.")
  }

  # --- Build design matrix ---
  group_factor <- factor(clinical[[group_column]], levels = group_levels)
  design <- stats::model.matrix(~ group_factor)
  colnames(design) <- c("Intercept", "GroupComparison")

  # --- Differential methylation with limma ---
  if (verbose) message("Running limma on ", nrow(meth_data), " probes and ",
                       ncol(meth_data), " samples...")

  fit <- limma::lmFit(meth_data, design)
  fit <- limma::eBayes(fit)

  results <- limma::topTable(fit, coef = "GroupComparison",
                             number = Inf, adjust.method = adj_method,
                             sort.by = "P")

  # --- Compute delta beta ---
  group1_samples <- clinical[[group_column]] == group_levels[1]
  group2_samples <- clinical[[group_column]] == group_levels[2]

  mean_group1 <- rowMeans(meth_data[, group1_samples, drop = FALSE], na.rm = TRUE)
  mean_group2 <- rowMeans(meth_data[, group2_samples, drop = FALSE], na.rm = TRUE)
  delta_beta_vec <- mean_group2 - mean_group1

  # --- Apply significance thresholds ---
  # Lower delta_beta (e.g., 0.05) -> more probes, higher sensitivity
  # Higher delta_beta (e.g., 0.2) -> fewer probes, higher specificity
  # Adjust p_cutoff similarly for statistical stringency
  results$probe <- rownames(results)
  results$delta_beta <- delta_beta_vec[rownames(results)]
  results$direction <- ifelse(results$delta_beta > 0, "Hyper", "Hypo")
  results$significant <- results$adj.P.Val < p_cutoff &
    abs(results$delta_beta) >= delta_beta

  # --- Reorder columns ---
  results <- results[, c("probe", "logFC", "AveExpr", "t", "P.Value",
                         "adj.P.Val", "delta_beta", "direction", "significant")]

  if (verbose) {
    n_sig <- sum(results$significant, na.rm = TRUE)
    message("Done. Found ", n_sig, " significant probes ",
            "(adj.P < ", p_cutoff, " & |delta_beta| >= ", delta_beta, ").")
  }

  return(results)
}
