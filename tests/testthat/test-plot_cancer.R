# Helper: build a data frame mimicking analyze_expression output
.make_test_deg <- function(n = 100) {
  set.seed(42)
  data.frame(
    gene = paste0("Gene", seq_len(n)),
    logFC = c(rnorm(20, 2, 0.5), rnorm(20, -2, 0.5), rnorm(n - 40, 0, 0.5)),
    AveExpr = rnorm(n, 8, 1),
    t = rnorm(n),
    P.Value = runif(n, 0, 0.05),
    adj.P.Val = runif(n, 0, 0.05),
    direction = c(rep("Up", 20), rep("Down", 20), rep("Up", n - 40)),
    significant = c(rep(TRUE, 40), rep(FALSE, n - 40)),
    stringsAsFactors = FALSE
  )
}


# Helper: build a data frame mimicking analyze_methylation output
.make_test_dmr <- function(n = 100) {
  set.seed(42)
  data.frame(
    probe = paste0("CpG", seq_len(n)),
    logFC = c(rnorm(20, 0.3, 0.05), rnorm(20, -0.3, 0.05), rnorm(n - 40, 0, 0.05)),
    AveExpr = runif(n, 0, 1),
    t = rnorm(n),
    P.Value = runif(n, 0, 0.05),
    adj.P.Val = runif(n, 0, 0.05),
    delta_beta = c(rnorm(20, 0.3, 0.05), rnorm(20, -0.3, 0.05), rnorm(n - 40, 0, 0.05)),
    direction = c(rep("Hyper", 20), rep("Hypo", 20), rep("Hyper", n - 40)),
    significant = c(rep(TRUE, 40), rep(FALSE, n - 40)),
    stringsAsFactors = FALSE
  )
}


# Helper: build a data frame mimicking integrate_meth_expr output
.make_test_integrated <- function(n = 100) {
  set.seed(42)
  data.frame(
    gene = paste0("Gene", seq_len(n)),
    probe = paste0("Gene", seq_len(n)),
    correlation = c(rnorm(30, -0.7, 0.1), rnorm(30, 0.7, 0.1), rnorm(n - 60, 0, 0.1)),
    p_value = runif(n, 0, 0.05),
    n_samples = rep(20, n),
    adj_p_value = runif(n, 0, 0.05),
    regulation = c(rep("Negative", 30), rep("Positive", 30), rep("Negative", n - 60)),
    significant = c(rep(TRUE, 60), rep(FALSE, n - 60)),
    stringsAsFactors = FALSE
  )
}


test_that("plot_cancer generates volcano plot for expression data", {
  deg <- .make_test_deg()
  p <- plot_cancer(deg, plot_type = "volcano", fc_cutoff = 1, verbose = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot_cancer generates volcano plot for methylation data", {
  dmr <- .make_test_dmr()
  p <- plot_cancer(dmr, plot_type = "volcano", fc_cutoff = 0.1, verbose = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot_cancer generates correlation plot", {
  integrated <- .make_test_integrated()
  p <- plot_cancer(integrated, plot_type = "correlation", verbose = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot_cancer applies custom title", {
  deg <- .make_test_deg()
  p <- plot_cancer(deg, plot_type = "volcano", title = "My Custom Title",
                   fc_cutoff = 1, verbose = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "My Custom Title")
})


test_that("plot_cancer errors when data is NULL for volcano plot", {
  expect_error(
    plot_cancer(data = NULL, plot_type = "volcano", verbose = FALSE),
    "'data' must be provided"
  )
})


test_that("plot_cancer errors when data is NULL for correlation plot", {
  expect_error(
    plot_cancer(data = NULL, plot_type = "correlation", verbose = FALSE),
    "'data' must be provided"
  )
})


test_that("plot_cancer errors on heatmap plot_type", {
  deg <- .make_test_deg()
  expect_error(
    plot_cancer(deg, plot_type = "heatmap", verbose = FALSE),
    "Heatmap plotting requires a MultiAssayExperiment"
  )
})


test_that("plot_cancer volcano errors when required columns missing", {
  bad_data <- data.frame(a = 1, b = 2)
  expect_error(
    plot_cancer(bad_data, plot_type = "volcano", verbose = FALSE),
    "'data' must contain columns"
  )
})


test_that("plot_cancer correlation errors when required columns missing", {
  bad_data <- data.frame(a = 1, b = 2)
  expect_error(
    plot_cancer(bad_data, plot_type = "correlation", verbose = FALSE),
    "'data' must contain columns 'correlation' and 'adj_p_value'"
  )
})


test_that("plot_cancer survival errors when surv_result is NULL", {
  expect_error(
    plot_cancer(plot_type = "survival", surv_result = NULL),
    "'surv_result' must be provided"
  )
})
