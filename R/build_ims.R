#' Build an Immune Methylation Score (IMS)
#'
#' Constructs a prognostic score from the methylation levels of CpGs
#' annotated to a panel of immune-related genes. The score is a weighted
#' sum of M-values, with weights estimated by an Elastic Net Cox model.
#' It is a \strong{prognostic risk score}, not a measure of tumour immune
#' infiltration: CpGs are selected by their association with survival,
#' not with immune cell content.
#'
#' @param mae A \code{MultiAssayExperiment} or \code{SummarizedExperiment}
#'   containing a methylation assay and clinical data.
#' @param assay_name Name of the methylation assay.
#' @param immune_genes Character vector of gene symbols. Defaults to a
#'   34-gene immune panel.
#' @param gene_map Optional \code{data.frame} with columns \code{cpg} and
#'   \code{gene}, or a named character vector. If \code{NULL}, the
#'   mapping is taken from \code{rowData(mae)[["gene"]]}.
#' @param region_col Column in \code{rowData} giving the genomic region.
#' @param regions Optional character vector of regions to keep.
#' @param min_sd Minimum standard deviation of M-values.
#' @param alpha Elastic Net mixing parameter (0 = ridge, 1 = lasso).
#' @param nfolds Number of inner cross-validation folds.
#' @param time_column Column with survival time.
#' @param event_column Column with event indicator.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical; print progress.
#'
#' @return An object of class \code{IMSResult} with elements
#'   \code{score}, \code{score_raw}, \code{selected_cpgs},
#'   \code{coefficients}, \code{model}, \code{gene_map},
#'   \code{n_cpgs_by_gene}, \code{parameters}, and \code{mae}.
#'
#' @examples
#' data(example_mae)
#' res <- build_ims(
#'     example_mae,
#'     immune_genes = c("CD274", "PDCD1", "BTLA", "TNF", "CD8A"),
#'     time_column = "OS.time",
#'     event_column = "OS",
#'     verbose = FALSE
#' )
#' res
#' res$selected_cpgs
#' @export
build_ims <- function(
    mae,
    assay_name = "methylation",
    immune_genes = .default_immune_genes(),
    gene_map = NULL,
    region_col = "region",
    regions = NULL,
    min_sd = 0.1,
    alpha = 0.5,
    nfolds = 5L,
    time_column = "OS.time",
    event_column = "OS",
    seed = NULL,
    verbose = TRUE
) {
    se <- .extract_assay(mae, assay_name)
    clin <- as.data.frame(MultiAssayExperiment::colData(mae))

    if (!time_column %in% colnames(clin)) {
        stop("Column '", time_column, "' not found in colData.", call. = FALSE)
    }
    if (!event_column %in% colnames(clin)) {
        stop("Column '", event_column, "' not found in colData.", call. = FALSE)
    }
    time <- as.numeric(clin[[time_column]])
    event <- as.numeric(clin[[event_column]])
    keep_s <- !is.na(time) & !is.na(event) & time > 0
    if (sum(keep_s) < 5) {
        stop("Fewer than 5 samples with valid survival data.", call. = FALSE)
    }

    if (is.null(gene_map)) {
        rd <- SummarizedExperiment::rowData(se)
        if (!"gene" %in% colnames(rd)) {
            stop("No 'gene' column in rowData and no gene_map supplied.",
                 call. = FALSE)
        }
        gene_map <- data.frame(
            cpg = rownames(se),
            gene = as.character(rd[["gene"]]),
            stringsAsFactors = FALSE
        )
    rownames(gene_map) <- gene_map$cpg
    } else if (is.character(gene_map) && !is.null(names(gene_map))) {
        gene_map <- data.frame(cpg = names(gene_map), gene = unname(gene_map),
                               stringsAsFactors = FALSE)
    } else if (is.data.frame(gene_map)) {
        if (!all(c("cpg", "gene") %in% colnames(gene_map))) {
            stop("gene_map must have columns 'cpg' and 'gene'.", call. = FALSE)
        }
    } else {
        stop("gene_map must be NULL, a data.frame, or a named vector.",
             call. = FALSE)
    }

    keep_region <- rep(TRUE, nrow(gene_map))
    if (!is.null(regions)) {
        rd <- SummarizedExperiment::rowData(se)
        if (!region_col %in% colnames(rd)) {
            stop("region_col '", region_col, "' not found in rowData.",
                 call. = FALSE)
        }
        keep_region <- as.character(rd[[region_col]]) %in% regions
    }

    keep_gene <- gene_map$gene %in% immune_genes
    keep <- keep_region & keep_gene
    if (sum(keep) < 3) {
        stop("Fewer than 3 CpGs match the immune gene panel.", call. = FALSE)
    }
    gene_map <- gene_map[keep, , drop = FALSE]

    if (verbose) {
        message(sprintf("CpGs in immune genes: %d", nrow(gene_map)))
    }

    beta <- SummarizedExperiment::assay(se, assay_name)
    beta <- beta[gene_map$cpg, , drop = FALSE]
    beta[beta <= 0] <- 0.001
    beta[beta >= 1] <- 0.999
    if (anyNA(beta)) {
        for (j in seq_len(ncol(beta))) {
            if (anyNA(beta[, j])) {
                beta[is.na(beta[, j]), j] <- stats::median(beta[, j], na.rm = TRUE)
            }
        }
    }
    mval <- log2(beta / (1 - beta))

    sd_m <- apply(mval, 1, stats::sd)
    keep_sd <- !is.na(sd_m) & sd_m > min_sd
    if (sum(keep_sd) < 3) {
        stop("Fewer than 3 CpGs pass the variance filter.", call. = FALSE)
    }
    mval <- mval[keep_sd, , drop = FALSE]
    gene_map <- gene_map[keep_sd, , drop = FALSE]

    if (verbose) {
        message(sprintf("CpGs after variance filter: %d", nrow(mval)))
    }

    X <- t(mval)
    X <- X[keep_s, , drop = FALSE]
    Y <- survival::Surv(time[keep_s], event[keep_s])

    if (!is.null(seed)) set.seed(seed)
    foldid <- sample(rep(seq_len(nfolds), length.out = nrow(X)))
    fit <- .cv_glmnet_safe(x = X, y = Y, family = "cox", cox.ties = "breslow",
                             alpha = alpha, foldid = foldid)

    coefs <- as.numeric(stats::coef(fit, s = "lambda.min"))
    names(coefs) <- colnames(X)
    selected <- names(coefs)[coefs != 0]
    if (length(selected) == 0) {
        stop("Elastic Net selected no CpGs at lambda.min.", call. = FALSE)
    }

    score_raw <- as.numeric(X %*% coefs)
    score <- as.numeric(scale(score_raw))
    names(score) <- rownames(X)
    names(score_raw) <- rownames(X)

    mae_out <- mae
    if (methods::is(mae, "MultiAssayExperiment")) {
        cd <- MultiAssayExperiment::colData(mae)
        cd[["ims_score"]] <- NA_real_
        cd[["ims_score"]][match(rownames(X), rownames(cd))] <- score
        MultiAssayExperiment::colData(mae_out) <- cd
    }

    out <- list(
        score = score,
        score_raw = score_raw,
        selected_cpgs = selected,
        coefficients = coefs[coefs != 0],
        model = fit,
        gene_map = gene_map,
        n_cpgs_by_gene = table(gene_map$gene),
        parameters = list(
            assay_name = assay_name,
            immune_genes = immune_genes,
            regions = regions,
            min_sd = min_sd,
            alpha = alpha,
            nfolds = nfolds,
            n_samples = nrow(X),
            n_events = sum(event[keep_s])
        ),
        mae = mae_out
    )
    class(out) <- c("IMSResult", "list")
    out
}

#' @export
print.IMSResult <- function(x, ...) {
    cat("IMSResult\n")
    cat("  Samples:", x$parameters$n_samples, "\n")
    cat("  Events:", x$parameters$n_events, "\n")
    cat("  CpGs selected:", length(x$selected_cpgs),
        "of", nrow(x$gene_map), "\n")
    cat("  Score range:",
        paste(round(range(x$score), 2), collapse = " to "), "\n")
    invisible(x)
}
