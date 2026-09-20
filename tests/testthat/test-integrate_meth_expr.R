# Helper to build a test MultiAssayExperiment with matched methylation + expression
.make_test_mae_integrated <- function(n_samples = 20, n_genes = 100) {
  set.seed(42)
  meth <- matrix(runif(n_genes * n_samples, 0, 1),
                 nrow = n_genes, ncol = n_samples)
  expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
                 nrow = n_genes, ncol = n_samples)
  rownames(meth) <- paste0("Gene", seq_len(n_genes))
  rownames(expr) <- paste0("Gene", seq_len(n_genes))
  colnames(meth) <- paste0("Sample", seq_len(n_samples))
  colnames(expr) <- paste0("Sample", seq_len(n_samples))

  # Induce strong negative correlation for first 20 genes
  for (i in 1:20) {
    expr[i, ] <- 15 - 5 * meth[i, ] + rnorm(n_samples, 0, 0.3)
  }

  se_meth <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = meth))
  se_expr <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = expr))

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(methylation = se_meth, expression = se_expr)
  )
}


test_that("integrate_meth_expr returns expected columns", {
  mae <- .make_test_mae_integrated()
  result <- integrate_meth_expr(mae, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expected_cols <- c("gene", "probe", "correlation", "p_value",
                     "n_samples", "adj_p_value", "regulation", "significant")
  expect_true(all(expected_cols %in% colnames(result)))
})


test_that("integrate_meth_expr detects negative correlations", {
  mae <- .make_test_mae_integrated()
  result <- integrate_meth_expr(mae, verbose = FALSE)

  sig_neg <- sum(result$significant & result$regulation == "Negative")
  expect_gte(sig_neg, 10)
})


test_that("integrate_meth_expr respects cor_cutoff", {
  mae <- .make_test_mae_integrated()

  result_low <- integrate_meth_expr(mae, cor_cutoff = 0.2, verbose = FALSE)
  result_high <- integrate_meth_expr(mae, cor_cutoff = 0.8, verbose = FALSE)

  expect_gte(sum(result_low$significant), sum(result_high$significant))
})


test_that("integrate_meth_expr supports spearman correlation", {
  mae <- .make_test_mae_integrated()
  result <- integrate_meth_expr(mae, correlation_method = "spearman",
                                verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})


test_that("integrate_meth_expr errors on non-MultiAssayExperiment input", {
  expect_error(
    integrate_meth_expr(data.frame(a = 1)),
    "'mae' must be a MultiAssayExperiment object"
  )
})


test_that("integrate_meth_expr errors on missing methylation assay", {
  mae <- .make_test_mae_integrated()
  expect_error(
    integrate_meth_expr(mae, meth_assay = "nonexistent", verbose = FALSE),
    "Methylation assay 'nonexistent' not found"
  )
})


test_that("integrate_meth_expr errors on missing expression assay", {
  mae <- .make_test_mae_integrated()
  expect_error(
    integrate_meth_expr(mae, expr_assay = "nonexistent", verbose = FALSE),
    "Expression assay 'nonexistent' not found"
  )
})


test_that("integrate_meth_expr errors on invalid correlation method", {
  mae <- .make_test_mae_integrated()
  expect_error(
    integrate_meth_expr(mae, correlation_method = "invalid", verbose = FALSE),
    "'correlation_method' must be one of"
  )
})


test_that("integrate_meth_expr errors when not enough common samples", {
  set.seed(42)
  meth <- matrix(runif(50 * 3), nrow = 50, ncol = 3)
  expr <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  rownames(meth) <- paste0("Gene", 1:50)
  rownames(expr) <- paste0("Gene", 1:50)
  colnames(meth) <- paste0("Sample", 1:3)
  colnames(expr) <- paste0("Sample", 1:3)

  se_meth <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = meth))
  se_expr <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = expr))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(methylation = se_meth, expression = se_expr)
  )

  expect_error(
    integrate_meth_expr(mae, min_samples = 5, verbose = FALSE),
    "Not enough common samples"
  )
})
