#' Integrate Methylation and Gene Expression Data
#'
#' Integrates DNA methylation and gene expression data to identify genes whose
#' expression is significantly correlated with their promoter methylation
#' status. This is a key step in multi-omics analysis for discovering
#' epigenetically regulated genes in any cancer type.
#'
#' @param mae A \code{MultiAssayExperiment} object created by
#'   \code{\link{load_cancer_data}}, containing both methylation and
#'   expression assays.
#' @param meth_assay Character. Name of the methylation assay. Default
#'   \code{"methylation"}.
#' @param expr_assay Character. Name of the expression assay. Default
#'   \code{"expression"}.
#' @param gene_map Optional. A data frame mapping methylation probes to genes.
#'   Must have columns \code{probe} and \code{gene}. If \code{NULL}, the
#'   function assumes probe and gene identifiers match (i.e., row names of
#'   both assays are gene symbols).
#' @param correlation_method Character. Correlation method, one of
#'   \code{"pearson"} (default) or \code{"spearman"}.
#' @param min_samples Integer. Minimum number of samples required for
#'   correlation analysis. Default \code{5}.
#' @param p_cutoff Numeric. Adjusted p-value threshold for significance.
#'   Default \code{0.05}.
#' @param cor_cutoff Numeric. Minimum absolute correlation coefficient to
#'   consider a methylation-expression pair significant. Default \code{0.3}.
#' @param verbose Logical. Whether to print progress messages. Default
#'   \code{TRUE}.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{gene}{Gene identifier}
#'   \item{probe}{Methylation probe identifier}
#'   \item{correlation}{Correlation coefficient between methylation and expression}
#'   \item{p_value}{Raw p-value of the correlation test}
#'   \item{adj_p_value}{Adjusted p-value (Benjamini-Hochberg)}
#'   \item{n_samples}{Number of samples used in the correlation}
#'   \item{regulation}{"Negative" if methylation is inversely correlated with
#'     expression (expected for promoter methylation), "Positive" otherwise}
#'   \item{significant}{Logical, whether the pair passes both p-value and
#'     correlation cutoffs}
#' }
#'
#' @examples
#' \dontrun{
#' # Load multi-omics data
#' brca <- load_cancer_data("brca_tcga",
#'                          data_types = c("methylation", "expression"))
#'
#' # Integrate methylation and expression
#' integrated <- integrate_meth_expr(brca)
#'
#' # Filter significant negative correlations
#' epi_genes <- integrated[integrated$significant &
#'                          integrated$regulation == "Negative", ]
#' }
#'
#' @export
integrate_meth_expr <- function(mae,
                                meth_assay = "methylation",
                                expr_assay = "expression",
                                gene_map = NULL,
                                correlation_method = "pearson",
                                min_samples = 5,
                                p_cutoff = 0.05,
                                cor_cutoff = 0.3,
                                verbose = TRUE) {

  # --- Validate inputs ---
  if (!methods::is(mae, "MultiAssayExperiment")) {
    stop("'mae' must be a MultiAssayExperiment object.")
  }

  available_assays <- names(MultiAssayExperiment::experiments(mae))

  if (!meth_assay %in% available_assays) {
    stop("Methylation assay '", meth_assay, "' not found. ",
         "Available assays: ", paste(available_assays, collapse = ", "))
  }

  if (!expr_assay %in% available_assays) {
    stop("Expression assay '", expr_assay, "' not found. ",
         "Available assays: ", paste(available_assays, collapse = ", "))
  }

  if (!correlation_method %in% c("pearson", "spearman")) {
    stop("'correlation_method' must be one of: 'pearson', 'spearman'.")
  }

  # --- Extract data ---
  if (verbose) message("Extracting methylation and expression data...")

  meth_se <- MultiAssayExperiment::experiments(mae)[[meth_assay]]
  expr_se <- MultiAssayExperiment::experiments(mae)[[expr_assay]]

  meth_data <- SummarizedExperiment::assay(meth_se)
  expr_data <- SummarizedExperiment::assay(expr_se)

  # --- Find common samples ---
  common_samples <- intersect(colnames(meth_data), colnames(expr_data))

  if (length(common_samples) < min_samples) {
    stop("Not enough common samples between methylation and expression data. ",
         "Found ", length(common_samples), " samples, need at least ",
         min_samples, ".")
  }

  if (verbose) message("Using ", length(common_samples), " common samples.")

  meth_data <- meth_data[, common_samples, drop = FALSE]
  expr_data <- expr_data[, common_samples, drop = FALSE]

  # --- Build probe-to-gene mapping ---
  if (is.null(gene_map)) {
    # Assume probe and gene identifiers match
    common_genes <- intersect(rownames(meth_data), rownames(expr_data))

    if (length(common_genes) == 0) {
      stop("No common identifiers between methylation and expression data. ",
           "Please provide a 'gene_map' data frame with columns 'probe' and 'gene'.")
    }

    if (verbose) message("No gene map provided. Using ", length(common_genes),
                         " matching identifiers directly.")

    gene_map <- data.frame(
      probe = common_genes,
      gene = common_genes,
      stringsAsFactors = FALSE
    )
  } else {
    if (!all(c("probe", "gene") %in% colnames(gene_map))) {
      stop("'gene_map' must have columns 'probe' and 'gene'.")
    }
    gene_map <- as.data.frame(gene_map, stringsAsFactors = FALSE)

    # Keep only probes and genes present in both datasets
    gene_map <- gene_map[gene_map$probe %in% rownames(meth_data) &
                           gene_map$gene %in% rownames(expr_data), ]

    if (nrow(gene_map) == 0) {
      stop("No valid probe-gene pairs after filtering. ",
           "Check that probe and gene identifiers match your data.")
    }
  }

  if (verbose) message("Testing ", nrow(gene_map), " probe-gene pairs...")

  # --- Correlation analysis ---
  results <- lapply(seq_len(nrow(gene_map)), function(i) {
    probe <- gene_map$probe[i]
    gene <- gene_map$gene[i]

    meth_vec <- as.numeric(meth_data[probe, ])
    expr_vec <- as.numeric(expr_data[gene, ])

    # Skip pairs with zero variance
    if (stats::var(meth_vec, na.rm = TRUE) == 0 ||
        stats::var(expr_vec, na.rm = TRUE) == 0) {
      return(NULL)
    }

    # Skip pairs with too many NAs
    valid <- !is.na(meth_vec) & !is.na(expr_vec)
    if (sum(valid) < min_samples) {
      return(NULL)
    }

    test <- tryCatch(
      stats::cor.test(meth_vec[valid], expr_vec[valid],
                      method = correlation_method),
      error = function(e) NULL
    )

    if (is.null(test)) return(NULL)

    data.frame(
      gene = gene,
      probe = probe,
      correlation = unname(test$estimate),
      p_value = test$p.value,
      n_samples = sum(valid),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results)

  if (is.null(results) || nrow(results) == 0) {
    stop("No valid correlation results. Check your data.")
  }

  # --- Multiple testing correction ---
  results$adj_p_value <- stats::p.adjust(results$p_value, method = "BH")

  # --- Annotate regulation direction ---
  # Negative correlation: higher methylation -> lower expression (typical for promoter methylation)
  results$regulation <- ifelse(results$correlation < 0, "Negative", "Positive")

  # --- Apply significance thresholds ---
  results$significant <- results$adj_p_value < p_cutoff &
    abs(results$correlation) >= cor_cutoff

  # --- Sort by adjusted p-value ---
  results <- results[order(results$adj_p_value), ]
  rownames(results) <- NULL

  if (verbose) {
    n_sig <- sum(results$significant, na.rm = TRUE)
    n_neg <- sum(results$significant & results$regulation == "Negative", na.rm = TRUE)
    message("Done. Found ", n_sig, " significant pairs (",
            n_neg, " with negative correlation, i.e., epigenetic silencing).")
  }

  return(results)
}
