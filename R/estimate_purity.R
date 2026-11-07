#' Estimate Tumor Purity and Immune Score
#'
#' @description
#' Computes ESTIMATE scores (Stromal, Immune, ESTIMATE) and Tumor Purity
#' from gene expression data using the ESTIMATE algorithm.
#' Requires package 'tidyestimate'.
#'
#' @param expr_data A matrix of gene expression (genes x samples).
#' @param is_affymetrix Logical. TRUE if data is from Affymetrix arrays.
#'   Default FALSE (RNA-seq).
#' @return A data frame with columns: sample, stromal, immune, estimate, purity.
#' @examples
#' # Source: Mock data simulated to demonstrate the ESTIMATE algorithm.
#' set.seed(123)
#' expr_mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(expr_mat) <- paste0("GENE", 1:100)
#' colnames(expr_mat) <- paste0("Sample", 1:10)
#' tryCatch({
#'   res <- estimate_purity(expr_mat)
#'   head(res)
#' }, error = function(e) {
#'   message("estimate_purity requires 'tidyestimate' package.")
#' })
#' @export
estimate_purity <- function(expr_data, is_affymetrix = FALSE) {
  if (!requireNamespace('tidyestimate', quietly = TRUE)) {
    stop('Please install tidyestimate: install.packages("tidyestimate")')
  }
  expr_mat <- as.matrix(expr_data)
  filtered <- tidyestimate::filter_common_genes(
    expr_mat, id = 'hgnc_symbol',
    tell_missing = FALSE, find_alias = TRUE
  )
  scores <- tidyestimate::estimate_score(filtered, is_affymetrix = is_affymetrix)
  scores_df <- as.data.frame(scores)
  if ('sample' %in% colnames(scores_df)) {
    rownames(scores_df) <- scores_df$sample
    scores_df <- scores_df[, -1]
  }
  scores_df$purity <- cos(0.6049872018 + 0.0001467884 * scores_df$estimate)
  scores_df
}

