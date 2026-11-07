test_that('bootstrap_validate_ims returns a valid object', {
    data(example_mae, package = 'MultiOmicsCancer')
    set.seed(42)
    bv <- bootstrap_validate_ims(
        example_mae,
        immune_genes = c('CD274', 'PDCD1', 'BTLA', 'TNF', 'CD8A'),
        n_boot = 5,
        verbose = FALSE
    )
    expect_s3_class(bv, 'BootstrapIMSResult')
    expect_true(is.numeric(bv$c_apparent))
    expect_true(is.numeric(bv$optimism))
    expect_true(is.numeric(bv$c_corrected))
    expect_false(is.nan(bv$optimism))
    expect_false(is.nan(bv$c_corrected))
    expect_true(bv$c_apparent >= 0 && bv$c_apparent <= 1)
    expect_equal(bv$n_boot, 5)
    expect_equal(length(bv$c_boot), 5)
})

test_that('bootstrap_validate_ims errors on missing columns', {
    data(example_mae, package = 'MultiOmicsCancer')
    expect_error(
        bootstrap_validate_ims(example_mae,
                               immune_genes = c('CD274', 'PDCD1', 'BTLA'),
                               n_boot = 2,
                               time_column = 'NONEXISTENT',
                               verbose = FALSE),
        'not found in colData'
    )
})
