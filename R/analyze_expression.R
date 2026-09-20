#' Differential Gene Expression Analysis for Pan-Cancer Studies
#'
#' Performs differential gene expression analysis between tumor and normal
#' samples (or any two groups defined in the clinical data) for any cancer
#' type. Supports both \code{limma} (for log-transformed or normalized data)
#' and \code{DESeq2} (for raw count data).
#'
#' @param mae A \code{MultiAssayExperiment} object created by
#'   \code{\link{load_cancer_data}}.
#' @param assay_name Character. Name of the expression assay in the
#'   \code{MultiAssayExperiment}. Default \code{"expression"}.
#' @param group_column Character. Name of the clinical column used to define
#'   groups. Default \code{"sample_type"}.
#' @param group_levels Character vector of length 2. The two groups to compare.
#'   The first is the reference (e.g., normal), the second is the comparison
#'   (e.g., tumor). If \code{NULL}, the function attempts to auto-detect
#'   common tumor/normal labels.
#' @param method Character. Method for differential expression. One of
#'   \code{"limma"} (default) or \code{"DESeq2"}. Use \code{"DESeq2"} only for
#'   raw count data.
#' @param log_transform Logical. If \code{TRUE} and \code{method = "limma"},
#'   applies \code{log2(x + 1)} transformation to the data. Default \code{TRUE}.
#' @param adj_method Character. Multiple testing correction method. Default
#'   \code{"BH"} (Benjamini-Hochberg).
#' @param p_cutoff Numeric. Adjusted p-value threshold for significance.
#'   Default \code{0.05}. Lower values increase specificity; higher values
#'   increase sensitivity.
#' @param logFC_cutoff Numeric. Minimum absolute log2 fold change to call a
#'   gene differentially expressed. Default \code{1} (i.e., 2-fold). Lower
#'   values (e.g., \code{0.5}) increase sensitivity and return more genes,
#'   but may include false positives. Higher values (e.g., \code{1.5}) increase
#'   specificity and return fewer genes.
#' @param verbose Logical. Whether to print progress messages. Default
#'   \code{TRUE}.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{gene}{Gene identifier}
#'   \item{logFC}{Log2 fold change (group2 - group1)}
#'   \item{AveExpr}{Average expression across all samples}
#'   \item{t}{Moderated t-statistic (limma) or Wald statistic (DESeq2)}
#'   \item{P.Value}{Raw p-value}
#'   \item{adj.P.Val}{Adjusted p-value}
#'   \item{direction}{"Up" or "Down" in group2 vs group1}
#'   \item{significant}{Logical, whether the gene passes both p-value and
#'     logFC cutoffs}
#' }
#'
#' @examples
#' \dontrun{
#' # Load data
#' brca <- load_cancer_data("brca_tcga")
#'
#' # Differential expression: normal vs tumor with limma
#' deg <- analyze_expression(brca, method = "limma")
#'
#' # With DESeq2 on raw counts
#' deg_deseq <- analyze_expression(brca, method = "DESeq2")
#'
#' # Filter significant genes
#' sig_deg <- deg[deg$significant, ]
#'
#' # Lower the logFC cutoff to get more genes (higher sensitivity)
#' deg_sensitive <- analyze_expression(brca, logFC_cutoff = 0.5)
#'
#' # Raise the logFC cutoff for stricter results (higher specificity)
#' deg_strict <- analyze_expression(brca, logFC_cutoff = 1.5)
#' }
#'
#' @export
analyze_expression <- function(mae,
                               assay_name = "expression",
                               group_column = "sample_type",
                               group_levels = NULL,
                               method = "limma",
                               log_transform = TRUE,
                               adj_method = "BH",
                               p_cutoff = 0.05,
                               logFC_cutoff = 1,
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

  if (!method %in% c("limma", "DESeq2")) {
    stop("'method' must be one of: 'limma', 'DESeq2'.")
  }

  # --- Extract data ---
  expr_se <- MultiAssayExperiment::experiments(mae)[[assay_name]]
  expr_data <- SummarizedExperiment::assay(expr_se)
  clinical <- as.data.frame(MultiAssayExperiment::colData(mae))

  if (!group_column %in% colnames(clinical)) {
    stop("Column '", group_column, "' not found in clinical data. ",
         "Available columns: ", paste(colnames(clinical), collapse = ", "))
  }

  # --- Determine groups ---
  if (is.null(group_levels)) {
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

  expr_data <- expr_data[, keep_samples, drop = FALSE]
  clinical <- clinical[keep_samples, , drop = FALSE]

  # --- Ensure sample order matches ---
  if (!identical(colnames(expr_data), rownames(clinical))) {
    common <- intersect(colnames(expr_data), rownames(clinical))
    if (length(common) == 0) {
      stop("Sample names in expression data and clinical data do not match.")
    }
    expr_data <- expr_data[, common, drop = FALSE]
    clinical <- clinical[common, , drop = FALSE]
  }

  # --- Remove zero-variance genes ---
  gene_var <- apply(expr_data, 1, stats::var, na.rm = TRUE)
  keep_genes <- !is.na(gene_var) & gene_var > 0
  expr_data <- expr_data[keep_genes, , drop = FALSE]

  if (nrow(expr_data) < 10) {
    stop("Too few informative genes after filtering.")
  }

  # --- Group factor ---
  group_factor <- factor(clinical[[group_column]], levels = group_levels)

  # --- Analysis ---
  if (method == "limma") {

    if (log_transform) {
      if (verbose) message("Applying log2(x + 1) transformation...")
      expr_data <- log2(expr_data + 1)
    }

    design <- stats::model.matrix(~ group_factor)
    colnames(design) <- c("Intercept", "GroupComparison")

    if (verbose) message("Running limma on ", nrow(expr_data), " genes and ",
                         ncol(expr_data), " samples...")

    fit <- limma::lmFit(expr_data, design)
    fit <- limma::eBayes(fit)

    results <- limma::topTable(fit, coef = "GroupComparison",
                               number = Inf, adjust.method = adj_method,
                               sort.by = "P")
  }

  if (method == "DESeq2") {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Package 'DESeq2' is required. Install with: ",
           "BiocManager::install('DESeq2')")
    }

    if (verbose) message("Running DESeq2 on ", nrow(expr_data), " genes and ",
                         ncol(expr_data), " samples...")

    # DESeq2 requires integer counts
    expr_data <- round(expr_data)

    col_data <- data.frame(
      condition = group_factor,
      row.names = colnames(expr_data)
    )

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = expr_data,
      colData = col_data,
      design = ~ condition
    )

    dds <- DESeq2::DESeq(dds, quiet = !verbose)
    res <- DESeq2::results(dds, contrast = c("condition", group_levels[2], group_levels[1]))

    results <- as.data.frame(res)
    results <- results[!is.na(results$padj), ]
    results <- results[order(results$padj), ]

    # Rename columns to match limma output
    colnames(results)[colnames(results) == "log2FoldChange"] <- "logFC"
    colnames(results)[colnames(results) == "pvalue"] <- "P.Value"
    colnames(results)[colnames(results) == "padj"] <- "adj.P.Val"

    # Compute average expression
    results$AveExpr <- rowMeans(log2(expr_data[rownames(results), , drop = FALSE] + 1))
    results$t <- results$stat
  }

  # --- Apply significance thresholds ---
  # Lower logFC_cutoff (e.g., 0.5) -> more genes, higher sensitivity
  # Higher logFC_cutoff (e.g., 1.5) -> fewer genes, higher specificity
  # Adjust p_cutoff similarly for statistical stringency
  results$gene <- rownames(results)
  results$direction <- ifelse(results$logFC > 0, "Up", "Down")
  results$significant <- results$adj.P.Val < p_cutoff &
    abs(results$logFC) >= logFC_cutoff

  # --- Reorder columns ---
  keep_cols <- c("gene", "logFC", "AveExpr", "t", "P.Value",
                 "adj.P.Val", "direction", "significant")
  results <- results[, keep_cols]

  if (verbose) {
    n_sig <- sum(results$significant, na.rm = TRUE)
    message("Done. Found ", n_sig, " significant genes ",
            "(adj.P < ", p_cutoff, " & |logFC| >= ", logFC_cutoff, ").")
  }

  return(results)
}
