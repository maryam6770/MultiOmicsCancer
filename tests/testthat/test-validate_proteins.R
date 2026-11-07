test_that('validate_proteins returns data frame', {
  res <- validate_proteins(c('CD8A', 'TNF'), use_uniprot = FALSE)
  expect_s3_class(res, 'data.frame')
  expect_true('protein_exists' %in% colnames(res))
  expect_true(nrow(res) == 2)
})

test_that('validate_proteins errors with empty input', {
  expect_error(validate_proteins(character(0)))
})

