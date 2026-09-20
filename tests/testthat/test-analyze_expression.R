# Helper to build a test MultiAssayExperiment with expression data
.make_test_mae_expr <- function(n_normal = 10, n_tumor = 10, n_genes = 100) {
  set.seed(123)
  expr <- matrix(rnorm(n_genes * (n_normal + n_tumor), mean = 8, sd = 1.5),
                 nrow = n_genes, ncol = n_normal + n_tumor)
  rownames(expr) <- paste0("Gene", seq_len(n_genes))
  colnames(expr) <- c(paste0("Normal", seq_len(n_normal)),
                      paste0("Tumor", seq_len(n_tumor)))
  # Induce strong upregulation in first 20 genes for tumor samples
  expr[1:20, (n_normal + 1):(n_normal + n_tumor)] <-
    expr[1:20, (n_normal + 1):(n_normal + n_tumor)] + 5

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = expr))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(expression = se)
  )
  clinical <- S4Vectors::DataFrame(
    sample_type = c(rep("Normal", n_normal), rep("Tumor", n_tumor))
  )
  rownames(clinical) <- colnames(expr)
  MultiAssayExperiment::colData(mae) <- clinical
  mae
}


test_that("analyze_expression returns expected columns", {
  mae <- .make_test_mae_expr()
  result <- analyze_expression(mae, log_transform = FALSE, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expected_cols <- c("gene", "logFC", "AveExpr", "t", "P.Value",
                     "adj.P.Val", "direction", "significant")
  expect_true(all(expected_cols %in% colnames(result)))
})


test_that("analyze_expression detects upregulated genes with limma", {
  mae <- .make_test_mae_expr()
  result <- analyze_expression(mae, method = "limma",
                               log_transform = FALSE, verbose = FALSE)

  sig_genes <- result$gene[result$significant]
  expect_gte(length(sig_genes), 5)

  top_sig <- result[result$significant, ]
  expect_gt(sum(top_sig$direction == "Up"), 3)
})


test_that("analyze_expression respects logFC_cutoff", {
  mae <- .make_test_mae_expr()

  result_low <- analyze_expression(mae, logFC_cutoff = 0.3,
                                   log_transform = FALSE, verbose = FALSE)
  result_high <- analyze_expression(mae, logFC_cutoff = 1.5,
                                    log_transform = FALSE, verbose = FALSE)

  expect_gte(sum(result_low$significant), sum(result_high$significant))
})


test_that("analyze_expression works with DESeq2 method", {
  set.seed(123)
  n_normal <- 10
  n_tumor <- 10
  n_genes <- 100
  expr_counts <- matrix(rpois(n_genes * (n_normal + n_tumor), lambda = 100),
                        nrow = n_genes, ncol = n_normal + n_tumor)
  rownames(expr_counts) <- paste0("Gene", seq_len(n_genes))
  colnames(expr_counts) <- c(paste0("Normal", seq_len(n_normal)),
                             paste0("Tumor", seq_len(n_tumor)))
  expr_counts[1:20, (n_normal + 1):(n_normal + n_tumor)] <-
    expr_counts[1:20, (n_normal + 1):(n_normal + n_tumor)] * 3

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = expr_counts)
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(expression = se)
  )
  clinical <- S4Vectors::DataFrame(
    sample_type = c(rep("Normal", n_normal), rep("Tumor", n_tumor))
  )
  rownames(clinical) <- colnames(expr_counts)
  MultiAssayExperiment::colData(mae) <- clinical

  result <- analyze_expression(mae, method = "DESeq2", verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true("gene" %in% colnames(result))
})


test_that("analyze_expression auto-detects group levels", {
  mae <- .make_test_mae_expr()
  result <- analyze_expression(mae, log_transform = FALSE, verbose = FALSE)
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})


test_that("analyze_expression errors on non-MultiAssayExperiment input", {
  expect_error(
    analyze_expression(data.frame(a = 1)),
    "'mae' must be a MultiAssayExperiment object"
  )
})


test_that("analyze_expression errors on missing assay", {
  mae <- .make_test_mae_expr()
  expect_error(
    analyze_expression(mae, assay_name = "nonexistent", verbose = FALSE),
    "Assay 'nonexistent' not found"
  )
})


test_that("analyze_expression errors on invalid method", {
  mae <- .make_test_mae_expr()
  expect_error(
    analyze_expression(mae, method = "invalid", verbose = FALSE),
    "'method' must be one of"
  )
})


test_that("analyze_expression errors on missing group_column", {
  mae <- .make_test_mae_expr()
  expect_error(
    analyze_expression(mae, group_column = "nonexistent", verbose = FALSE),
    "Column 'nonexistent' not found"
  )
})
