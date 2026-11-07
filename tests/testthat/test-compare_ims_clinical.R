test_that('compare_ims_clinical returns a valid object', {
    data(example_mae, package = 'MultiOmicsCancer')
    set.seed(42)
    ims <- build_ims(
        example_mae,
        immune_genes = c('CD274', 'PDCD1', 'BTLA', 'TNF', 'CD8A'),
        verbose = FALSE
    )
    cmp <- compare_ims_clinical(
        example_mae, ims,
        clinical_vars = 'sample_type',
        verbose = FALSE
    )
    expect_s3_class(cmp, 'IMSComparison')
    expect_true(is.numeric(cmp$c_clinical))
    expect_true(is.numeric(cmp$c_combined))
    expect_true(is.numeric(cmp$delta_c))
    expect_equal(cmp$delta_c, cmp$c_combined - cmp$c_clinical, tolerance = 1e-10)
    expect_true(cmp$c_clinical >= 0 && cmp$c_clinical <= 1)
    expect_true(cmp$c_combined >= 0 && cmp$c_combined <= 1)
    expect_true(is.numeric(cmp$lrt_stat))
    expect_true(is.numeric(cmp$lrt_p))
    expect_true(cmp$lrt_p >= 0 && cmp$lrt_p <= 1)
    expect_equal(cmp$n_samples, 40)
})

test_that('compare_ims_clinical works without clinical variables', {
    data(example_mae, package = 'MultiOmicsCancer')
    set.seed(42)
    ims <- build_ims(
        example_mae,
        immune_genes = c('CD274', 'PDCD1', 'BTLA', 'TNF', 'CD8A'),
        verbose = FALSE
    )
    cmp <- compare_ims_clinical(
        example_mae, ims,
        clinical_vars = NULL,
        verbose = FALSE
    )
    expect_s3_class(cmp, 'IMSComparison')
    expect_true(is.numeric(cmp$delta_c))
})

test_that('compare_ims_clinical errors on non-IMSResult input', {
    data(example_mae, package = 'MultiOmicsCancer')
    expect_error(
        compare_ims_clinical(example_mae, list(score = 1:40)),
        'ims_result must be an IMSResult'
    )
})
