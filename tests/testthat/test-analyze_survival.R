# Helper to build a test MultiAssayExperiment with survival data
.make_test_mae_surv <- function(n_samples = 30, n_genes = 50) {
  set.seed(42)
  expr <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
                 nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", seq_len(n_genes))
  colnames(expr) <- paste0("Sample", seq_len(n_samples))

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = expr))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(expression = se)
  )
  clinical <- S4Vectors::DataFrame(
    OS.time = sample(100:2000, n_samples),
    OS = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.3, 0.7))
  )
  rownames(clinical) <- colnames(expr)
  MultiAssayExperiment::colData(mae) <- clinical
  mae
}


test_that("analyze_survival returns expected structure", {
  mae <- .make_test_mae_surv()
  result <- analyze_survival(mae, feature_name = "Gene1", verbose = FALSE)

  expect_type(result, "list")
  expect_true("km_fit" %in% names(result))
  expect_true("logrank_p" %in% names(result))
  expect_true("cox_model" %in% names(result))
  expect_true("data" %in% names(result))
})


test_that("analyze_survival computes valid log-rank p-value", {
  mae <- .make_test_mae_surv()
  result <- analyze_survival(mae, feature_name = "Gene1", verbose = FALSE)

  expect_true(is.numeric(result$logrank_p))
  expect_gte(result$logrank_p, 0)
  expect_lte(result$logrank_p, 1)
})


test_that("analyze_survival splits samples into High and Low groups", {
  mae <- .make_test_mae_surv()
  result <- analyze_survival(mae, feature_name = "Gene1", verbose = FALSE)

  groups <- result$data$group
  expect_true("High" %in% as.character(groups))
  expect_true("Low" %in% as.character(groups))
})


test_that("analyze_survival works with tertile split", {
  mae <- .make_test_mae_surv()
  result <- analyze_survival(mae, feature_name = "Gene1",
                             split_method = "tertile", verbose = FALSE)
  expect_type(result, "list")
  expect_true("km_fit" %in% names(result))
})


test_that("analyze_survival works with quartile split", {
  mae <- .make_test_mae_surv()
  result <- analyze_survival(mae, feature_name = "Gene1",
                             split_method = "quartile", verbose = FALSE)
  expect_type(result, "list")
})


test_that("analyze_survival returns clinical summary when no feature given", {
  mae <- .make_test_mae_surv()
  result <- analyze_survival(mae, verbose = FALSE)

  expect_type(result, "list")
  expect_true("km_fit" %in% names(result))
  expect_true("message" %in% names(result))
})


test_that("analyze_survival errors on non-MultiAssayExperiment input", {
  expect_error(
    analyze_survival(data.frame(a = 1)),
    "'mae' must be a MultiAssayExperiment object"
  )
})


test_that("analyze_survival errors on invalid split_method", {
  mae <- .make_test_mae_surv()
  expect_error(
    analyze_survival(mae, feature_name = "Gene1",
                     split_method = "invalid", verbose = FALSE),
    "'split_method' must be one of"
  )
})


test_that("analyze_survival errors on missing time column", {
  mae <- .make_test_mae_surv()
  expect_error(
    analyze_survival(mae, time_column = "nonexistent",
                     feature_name = "Gene1", verbose = FALSE),
    "Time column 'nonexistent' not found"
  )
})


test_that("analyze_survival errors on missing event column", {
  mae <- .make_test_mae_surv()
  expect_error(
    analyze_survival(mae, event_column = "nonexistent",
                     feature_name = "Gene1", verbose = FALSE),
    "Event column 'nonexistent' not found"
  )
})


test_that("analyze_survival errors on missing feature", {
  mae <- .make_test_mae_surv()
  expect_error(
    analyze_survival(mae, feature_name = "NonexistentGene", verbose = FALSE),
    "Feature 'NonexistentGene' not found"
  )
})


test_that("analyze_survival errors on missing assay", {
  mae <- .make_test_mae_surv()
  expect_error(
    analyze_survival(mae, feature_name = "Gene1",
                     assay_name = "nonexistent", verbose = FALSE),
    "Assay 'nonexistent' not found"
  )
})
