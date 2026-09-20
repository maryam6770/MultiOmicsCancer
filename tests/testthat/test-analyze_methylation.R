# Helper to build a test MultiAssayExperiment with methylation data
.make_test_mae_meth <- function(n_normal = 10, n_tumor = 10, n_probes = 100) {
  set.seed(42)
  meth <- matrix(runif(n_probes * (n_normal + n_tumor), 0, 1),
                 nrow = n_probes, ncol = n_normal + n_tumor)
  rownames(meth) <- paste0("CpG", seq_len(n_probes))
  colnames(meth) <- c(paste0("Normal", seq_len(n_normal)),
                      paste0("Tumor", seq_len(n_tumor)))
  # Induce strong hypermethylation in first 20 probes for tumor samples
  meth[1:20, (n_normal + 1):(n_normal + n_tumor)] <-
    pmin(meth[1:20, (n_normal + 1):(n_normal + n_tumor)] + 0.4, 1)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = meth))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(methylation = se)
  )
  clinical <- S4Vectors::DataFrame(
    sample_type = c(rep("Normal", n_normal), rep("Tumor", n_tumor))
  )
  rownames(clinical) <- colnames(meth)
  MultiAssayExperiment::colData(mae) <- clinical
  mae
}


test_that("analyze_methylation returns expected columns", {
  mae <- .make_test_mae_meth()
  result <- analyze_methylation(mae, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expected_cols <- c("probe", "logFC", "AveExpr", "t", "P.Value",
                     "adj.P.Val", "delta_beta", "direction", "significant")
  expect_true(all(expected_cols %in% colnames(result)))
})


test_that("analyze_methylation detects hypermethylated probes", {
  mae <- .make_test_mae_meth()
  result <- analyze_methylation(mae, verbose = FALSE)

  # The 20 probes we deliberately modified should be significant
  sig_probes <- result$probe[result$significant]
  expect_gte(length(sig_probes), 5)

  # Most significant probes should be "Hyper"
  top_sig <- result[result$significant, ]
  expect_gt(sum(top_sig$direction == "Hyper"), 3)
})


test_that("analyze_methylation respects delta_beta cutoff", {
  mae <- .make_test_mae_meth()

  result_low <- analyze_methylation(mae, delta_beta = 0.05, verbose = FALSE)
  result_high <- analyze_methylation(mae, delta_beta = 0.3, verbose = FALSE)

  expect_gte(sum(result_low$significant), sum(result_high$significant))
})


test_that("analyze_methylation auto-detects group levels", {
  mae <- .make_test_mae_meth()
  result <- analyze_methylation(mae, verbose = FALSE)

  # Should detect Normal vs Tumor automatically
  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})


test_that("analyze_methylation accepts explicit group_levels", {
  mae <- .make_test_mae_meth()
  result <- analyze_methylation(
    mae,
    group_levels = c("Normal", "Tumor"),
    verbose = FALSE
  )
  expect_s3_class(result, "data.frame")
})


test_that("analyze_methylation errors on non-MultiAssayExperiment input", {
  expect_error(
    analyze_methylation(data.frame(a = 1)),
    "'mae' must be a MultiAssayExperiment object"
  )
})


test_that("analyze_methylation errors on missing assay", {
  mae <- .make_test_mae_meth()
  expect_error(
    analyze_methylation(mae, assay_name = "nonexistent", verbose = FALSE),
    "Assay 'nonexistent' not found"
  )
})


test_that("analyze_methylation errors on missing group_column", {
  mae <- .make_test_mae_meth()
  expect_error(
    analyze_methylation(mae, group_column = "nonexistent", verbose = FALSE),
    "Column 'nonexistent' not found"
  )
})


test_that("analyze_methylation errors on invalid method", {
  mae <- .make_test_mae_meth()
  expect_error(
    analyze_methylation(mae, method = "invalid", verbose = FALSE),
    "'method' must be 'limma'"
  )
})


test_that("analyze_methylation errors on invalid group_levels length", {
  mae <- .make_test_mae_meth()
  expect_error(
    analyze_methylation(mae, group_levels = c("Normal"), verbose = FALSE),
    "'group_levels' must be a character vector of length 2"
  )
})
