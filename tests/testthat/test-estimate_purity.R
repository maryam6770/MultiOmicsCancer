test_that('estimate_purity works', {
  skip_if_not_installed('tidyestimate')
  set.seed(1)
  n <- 30
  genes <- c('CD8A', 'CD4', 'CD19', 'MS4A1', 'CD68', 'PTPRC',
             'COL1A1', 'COL3A1', 'ACTA2', 'VIM')
  expr <- matrix(runif(length(genes) * n, 0, 10),
                 nrow = length(genes), dimnames = list(genes, paste0('S', 1:n)))
  res <- estimate_purity(expr)
  expect_s3_class(res, 'data.frame')
  expect_true('purity' %in% colnames(res))
})

