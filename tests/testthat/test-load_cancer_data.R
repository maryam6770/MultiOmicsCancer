test_that("load_cancer_data works with custom data", {
  set.seed(42)
  meth <- matrix(runif(200), nrow = 20, ncol = 10)
  expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(meth) <- paste0("Gene", 1:20)
  rownames(expr) <- paste0("Gene", 1:20)
  colnames(meth) <- paste0("Sample", 1:10)
  colnames(expr) <- paste0("Sample", 1:10)

  mae <- load_cancer_data(
    cancer_type = "test",
    data_types = c("methylation", "expression"),
    source = "custom",
    custom_data = list(methylation = meth, expression = expr),
    verbose = FALSE
  )

  expect_s4_class(mae, "MultiAssayExperiment")
  expect_equal(length(MultiAssayExperiment::experiments(mae)), 2)
  expect_true("methylation" %in% names(MultiAssayExperiment::experiments(mae)))
  expect_true("expression" %in% names(MultiAssayExperiment::experiments(mae)))
})


test_that("load_cancer_data errors on invalid cancer_type", {
  expect_error(
    load_cancer_data(cancer_type = 123, source = "custom"),
    "'cancer_type' must be a single character string"
  )
})


test_that("load_cancer_data errors on invalid data_types", {
  expect_error(
    load_cancer_data(
      cancer_type = "test",
      data_types = c("invalid_type"),
      source = "custom",
      custom_data = list()
    ),
    "'data_types' must be a subset of"
  )
})


test_that("load_cancer_data errors on invalid source", {
  expect_error(
    load_cancer_data(cancer_type = "test", source = "invalid_source"),
    "'source' must be one of"
  )
})


test_that("load_cancer_data errors when custom_data is NULL", {
  expect_error(
    load_cancer_data(cancer_type = "test", source = "custom"),
    "'custom_data' must be provided"
  )
})


test_that("load_cancer_data errors when custom_data is missing requested types", {
  meth <- matrix(runif(100), nrow = 10, ncol = 10)
  expect_error(
    load_cancer_data(
      cancer_type = "test",
      data_types = c("methylation", "expression"),
      source = "custom",
      custom_data = list(methylation = meth)
    ),
    "'custom_data' must contain all requested data types"
  )
})
