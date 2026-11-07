test_that('build_ims returns a valid IMSResult object', {
    data(example_mae, package = 'MultiOmicsCancer')
    res <- build_ims(
        example_mae,
        immune_genes = c('CD274', 'PDCD1', 'BTLA', 'TNF', 'CD8A'),
        verbose = FALSE
    )
    expect_s3_class(res, 'IMSResult')
    expect_type(res$score, 'double')
    expect_equal(length(res$score), 40)
    expect_true(length(res$selected_cpgs) > 0)
    expect_true(all(res$selected_cpgs %in% rownames(res$gene_map)))
    expect_named(res, c('score', 'score_raw', 'selected_cpgs',
                        'coefficients', 'model', 'gene_map',
                        'n_cpgs_by_gene', 'parameters', 'mae'))
})

test_that('build_ims fails with informative errors', {
    data(example_mae, package = 'MultiOmicsCancer')
    expect_error(
        build_ims(example_mae, immune_genes = 'NONEXISTENT_GENE', verbose = FALSE),
        'Fewer than 3 CpGs'
    )
    expect_error(
        build_ims(example_mae, immune_genes = c('CD274', 'PDCD1'),
                  time_column = 'NONEXISTENT', verbose = FALSE),
        'not found in colData'
    )
})

test_that('build_ims respects seed', {
    data(example_mae, package = 'MultiOmicsCancer')
    res1 <- build_ims(example_mae,
                      immune_genes = c('CD274', 'PDCD1', 'BTLA', 'TNF', 'CD8A'),
                      seed = 42, verbose = FALSE)
    res2 <- build_ims(example_mae,
                      immune_genes = c('CD274', 'PDCD1', 'BTLA', 'TNF', 'CD8A'),
                      seed = 42, verbose = FALSE)
    expect_equal(res1$score, res2$score)
})

test_that('build_ims score is standardized', {
    data(example_mae, package = 'MultiOmicsCancer')
    res <- build_ims(example_mae,
                     immune_genes = c('CD274', 'PDCD1', 'BTLA', 'TNF', 'CD8A'),
                     verbose = FALSE)
    expect_equal(mean(res$score), 0, tolerance = 1e-10)
    expect_equal(sd(res$score), 1, tolerance = 1e-10)
})
